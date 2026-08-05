#!/usr/bin/env python3
"""
Chord diagram (circlize style) in pure matplotlib.
Converts R circlize::chordDiagram → Python matplotlib.

Usage:
    python tax_chord.py [tax_file] [group_file] [--top 10] [--output chord.png]

Input:
    tax_file   : species abundance table (TSV, rows=species, cols=samples, first col=species names)
    group_file : sample grouping (TSV, col1=sample, col2=group)
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patheffects as pe

# ─── Global style ───
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'font.size': 10,
    'axes.unicode_minus': False,
})

# ─── Color palettes ───
# Group colors (left sectors)
GROUP_COLORS = {
    'KC': '#EF9A9A',
    'WB': '#90CAF9',
    'MC': '#F3D32C',
}

# Species colors (right sectors + chords)
SPECIES_COLORS = {
    's__Ruminococcus_bromii':         '#8DD3C7',
    's__Megamonas_unclassified':       '#BEBADA',
    's__Bacteroides_dorei':            '#009329',
    's__Escherichia_coli':             '#FF7F00',
    's__Eubacterium_rectale':          '#CC3366',
    's__Bacteroides_uniformis':        '#c5c5f2',
    's__Bifidobacterium_longum':       '#FB8072',
    's__Roseburia_intestinalis':       '#F781BF',
    's__Bacteroides_ovatus':           '#80B1D3',
    's__Faecalibacterium_prausnitzii': '#A65628',
    'Others':                           '#A4566980',
}

# Fallback palette for species not in the dict
_FALLBACK = ['#8DD3C7','#BEBADA','#009329','#FF7F00','#CC3366',
             '#c5c5f2','#FB8072','#F781BF','#80B1D3','#A65628',
             '#66C2A5','#FC8D62','#E78AC3','#A6D854','#FFD92F']


def get_species_color(name, idx):
    if name in SPECIES_COLORS:
        return SPECIES_COLORS[name]
    if name == 'Others':
        return '#A4566980'
    return _FALLBACK[idx % len(_FALLBACK)]


def get_group_color(name):
    if name in GROUP_COLORS:
        return GROUP_COLORS[name]
    # Auto-generate pastel
    import colorsys
    h = hash(name) % 360 / 360.0
    r, g, b = colorsys.hls_to_rgb(h, 0.75, 0.6)
    return f'#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}'


def process_data(tax_file, group_file, top_n=10):
    """Read data, compute group-mean relative abundance, return top_n + Others."""
    data = pd.read_csv(tax_file, sep='\t', index_col=0)
    group = pd.read_csv(group_file, sep='\t')

    # Standardize group file columns
    if group.shape[1] == 2:
        group.columns = ['sample', 'group']
    else:
        group = group.iloc[:, [0, 1]]
        group.columns = ['sample', 'group']

    # Sort by total abundance
    data['_sum'] = data.sum(axis=1)
    data = data.sort_values('_sum', ascending=False)
    data = data.drop(columns=['_sum'])

    # Group means
    groups = group['group'].unique()
    group_means = pd.DataFrame(index=data.index)
    for g in groups:
        g_samples = group.loc[group['group'] == g, 'sample'].tolist()
        g_samples = [s for s in g_samples if s in data.columns]
        group_means[g] = data[g_samples].mean(axis=1)

    # Top N + Others
    top = group_means.iloc[:top_n]
    others_sum = group_means.iloc[top_n:].sum(axis=0)
    result = top.copy()
    result.loc['Others'] = others_sum

    # Relative abundance per group
    for g in result.columns:
        result[g] = result[g] / result[g].sum()

    return result


def draw_chord(mat, output_path, figsize=(10, 10), dpi=600, gap_deg=2.0, transparency=0.7):
    """
    Draw a chord diagram from a matrix.
    mat: DataFrame, rows=species (including 'Others'), columns=groups.
    """
    groups = list(mat.columns)          # left sectors
    species = list(mat.index)           # right sectors

    # ── Compute sector angular spans ──
    # Each sector's angle is proportional to its total (row sum for species, col sum for groups)
    group_totals = mat.sum(axis=0)      # Series indexed by group
    species_totals = mat.sum(axis=1)    # Series indexed by species

    all_sectors = groups + species      # order: groups on left, species on right
    all_totals = pd.concat([group_totals, species_totals])

    n_sectors = len(all_sectors)
    total_gap = gap_deg * n_sectors
    available_deg = 360.0 - total_gap
    total_value = all_totals.sum()

    sector_angles = {}
    current_angle = 90.0  # start from top

    for name in all_sectors:
        span = (all_totals[name] / total_value) * available_deg
        start = current_angle
        end = current_angle + span
        sector_angles[name] = (start, end)
        current_angle = end + gap_deg

    # ── Setup figure ──
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'aspect': 'equal'})
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.55, 1.55)
    ax.axis('off')

    R_inner = 0.82
    R_outer = 0.90
    R_label = 0.96

    # ── Draw sector arcs ──
    for name in all_sectors:
        start, end = sector_angles[name]
        is_group = name in groups
        color = get_group_color(name) if is_group else get_species_color(name, species.index(name))

        theta1, theta2 = start, end
        wedge = mpatches.Wedge(
            center=(0, 0), r=R_outer, theta1=theta1, theta2=theta2,
            width=R_outer - R_inner, facecolor=color, edgecolor='white', linewidth=0.5
        )
        ax.add_patch(wedge)

        # ── Sector label: horizontal, left/right aligned by position ──
        mid_angle = np.radians((theta1 + theta2) / 2)
        # Point on the outer edge of the sector
        edge_x = R_outer * np.cos(mid_angle)
        edge_y = R_outer * np.sin(mid_angle)

        # Italic for species (Latin names)
        if not is_group and name != 'Others':
            label_text = name.replace('s__', '').replace('_', ' ')
            style = 'italic'
        elif name == 'Others':
            label_text = 'Others'
            style = 'normal'
        else:
            label_text = name
            style = 'normal'

        # Place label horizontally: right side → left-align, left side → right-align
        mid_deg = (theta1 + theta2) / 2 % 360
        if mid_deg <= 180:
            # Right half of circle
            lx = R_label * np.cos(mid_angle) + 0.02
            ha = 'left'
        else:
            # Left half of circle
            lx = R_label * np.cos(mid_angle) - 0.02
            ha = 'right'
        ly = R_label * np.sin(mid_angle)

        # Draw thin leader line from sector edge to label
        ax.plot([edge_x, lx], [edge_y, ly], color='gray', linewidth=0.4, zorder=0)

        ax.text(lx, ly, label_text, fontsize=6.5, fontstyle=style,
                fontweight='bold' if is_group else 'normal',
                ha=ha, va='center')

    # ── Draw chords ──
    # For each (group, species) pair, draw a bezier chord
    for gi, grp in enumerate(groups):
        grp_start, grp_end = sector_angles[grp]
        # Track used angle within this group sector for stacking
        grp_cursor = grp_start

        for si, spc in enumerate(species):
            val = mat.loc[spc, grp]
            if val < 1e-6:
                continue

            spc_start, spc_end = sector_angles[spc]
            spc_span = spc_end - spc_start

            # Proportion of this species in the group
            prop_in_group = val / group_totals[grp]
            chord_span_in_group = prop_in_group * (grp_end - grp_start)

            # Position within group sector
            g_s = grp_cursor
            g_e = grp_cursor + chord_span_in_group
            grp_cursor = g_e

            # Position within species sector
            # Proportion of this group in the species
            prop_in_species = val / species_totals[spc]
            # We need to track per-species cursor too
            # Build a cumulative dict
            if not hasattr(draw_chord, '_spc_cursors'):
                draw_chord._spc_cursors = {}
            if spc not in draw_chord._spc_cursors:
                draw_chord._spc_cursors[spc] = spc_start

            s_s = draw_chord._spc_cursors[spc]
            chord_span_in_species = prop_in_species * spc_span
            s_e = s_s + chord_span_in_species
            draw_chord._spc_cursors[spc] = s_e

            # Draw bezier chord
            chord_color = get_species_color(spc, si)
            _draw_single_chord(ax, g_s, g_e, s_s, s_e, R_inner, chord_color, transparency)

    # Reset cursor cache
    if hasattr(draw_chord, '_spc_cursors'):
        del draw_chord._spc_cursors

    # ── Legend (species only) ──
    legend_handles = []
    for spc in species:
        color = get_species_color(spc, species.index(spc))
        if spc == 'Others':
            label = 'Others'
        else:
            label = spc.replace('s__', '').replace('_', ' ')
        legend_handles.append(mpatches.Patch(facecolor=color, edgecolor='none', label=label))

    leg = ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
                    fontsize=7, frameon=False, title='Species', title_fontsize=8,
                    prop={'style': 'italic'})
    # "Others" and title should not be italic
    for t in leg.get_texts():
        if t.get_text() == 'Others':
            t.set_fontstyle('normal')
    leg.get_title().set_fontstyle('normal')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")


def _draw_single_chord(ax, g_start, g_end, s_start, s_end, R, color, alpha):
    """Draw one chord as a filled bezier polygon from group arc to species arc."""
    n_pts = 30

    # Group arc points
    g_angles = np.linspace(np.radians(g_start), np.radians(g_end), n_pts)
    g_x = R * np.cos(g_angles)
    g_y = R * np.sin(g_angles)

    # Species arc points (reversed so polygon closes properly)
    s_angles = np.linspace(np.radians(s_end), np.radians(s_start), n_pts)
    s_x = R * np.cos(s_angles)
    s_y = R * np.sin(s_angles)

    # Control points for bezier curves (through center)
    g_mid_angle = np.radians((g_start + g_end) / 2)
    s_mid_angle = np.radians((s_start + s_end) / 2)

    # Build path: group arc → curve to species → species arc → curve back to group
    verts = []
    codes = []

    # Start at first group arc point
    verts.append((g_x[0], g_y[0]))
    codes.append(Path.MOVETO)

    # Group arc
    for i in range(1, n_pts):
        verts.append((g_x[i], g_y[i]))
        codes.append(Path.LINETO)

    # Bezier from last group point to first species point
    ctrl_r = R * 0.3  # control point radius (pull toward center)
    ctrl1 = (ctrl_r * np.cos(g_mid_angle), ctrl_r * np.sin(g_mid_angle))
    ctrl2 = (ctrl_r * np.cos(s_mid_angle), ctrl_r * np.sin(s_mid_angle))

    verts.append(ctrl1)
    codes.append(Path.CURVE3)
    verts.append((s_x[0], s_y[0]))
    codes.append(Path.CURVE3)

    # Species arc
    for i in range(1, n_pts):
        verts.append((s_x[i], s_y[i]))
        codes.append(Path.LINETO)

    # Bezier from last species point back to first group point
    verts.append(ctrl2)
    codes.append(Path.CURVE3)
    verts.append((g_x[0], g_y[0]))
    codes.append(Path.CURVE3)

    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor='none',
                                alpha=alpha, linewidth=0)
    ax.add_patch(patch)


def main():
    parser = argparse.ArgumentParser(description='Chord diagram from species abundance data')
    parser.add_argument('tax_file', nargs='?', default='chord_tax.txt',
                        help='Species abundance table (TSV)')
    parser.add_argument('group_file', nargs='?', default='chord_group.txt',
                        help='Group file (TSV: sample, group)')
    parser.add_argument('--top', type=int, default=10,
                        help='Top N species to display (default: 10)')
    parser.add_argument('--output', '-o', default='chord_diagram.png',
                        help='Output image path')
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--gap', type=float, default=2.0,
                        help='Gap between sectors in degrees (default: 2)')
    parser.add_argument('--alpha', type=float, default=0.7,
                        help='Chord transparency (default: 0.7)')
    args = parser.parse_args()

    mat = process_data(args.tax_file, args.group_file, top_n=args.top)
    print(f"Data matrix: {mat.shape[0]} species × {mat.shape[1]} groups")
    print(mat.round(4))

    draw_chord(mat, args.output, dpi=args.dpi, gap_deg=args.gap, transparency=args.alpha)


if __name__ == '__main__':
    main()
