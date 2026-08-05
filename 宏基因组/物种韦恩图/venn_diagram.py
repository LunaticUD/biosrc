#!/usr/bin/env python3
"""
Venn Diagram - Python matplotlib implementation
Converted from R VennDiagram code.
Supports 2-4 groups. Times New Roman font.
Shows raw count + percentage in each region.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from matplotlib import rcParams
import itertools
import sys

# ===== Font Configuration =====
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['mathtext.fontset'] = 'stix'

# ===== Color Palette (RColorBrewer Dark2) =====
DARK2 = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a',
         '#66a61e', '#e6ab02', '#a6761d', '#666666']

# ===== Read Data =====
data = pd.read_csv("tax.txt", sep="\t", index_col=0)
group_info = pd.read_csv("group.txt", sep="\t")

sample_col = group_info.columns[0]  # first column = sample IDs
groups = group_info['group'].unique().tolist()
n_groups = len(groups)

# ===== Compute Presence Sets per Group =====
group_sets = {}
for grp in groups:
    samples = group_info[group_info['group'] == grp][sample_col].tolist()
    common = [s for s in samples if s in data.columns]
    if common:
        row_sum = data[common].sum(axis=1)
        group_sets[grp] = set(row_sum[row_sum > 0].index)
    else:
        group_sets[grp] = set()

total = len(set.union(*group_sets.values())) if group_sets else 0
print(f"Groups: {groups}")
print(f"Total unique taxa: {total}")
for grp in groups:
    print(f"  {grp}: {len(group_sets[grp])} taxa")

# ===== Compute Exclusive Intersection Sizes =====
# Each region = taxa present in ALL listed groups and NONE of the others
regions = {}
for r in range(1, n_groups + 1):
    for combo in itertools.combinations(range(n_groups), r):
        sets_in = [group_sets[groups[i]] for i in combo]
        intersect = set.intersection(*sets_in)
        for j in range(n_groups):
            if j not in combo:
                intersect -= group_sets[groups[j]]
        regions[combo] = len(intersect)

# ===== Format Region Label (raw count + percentage) =====
def fmt_label(count):
    if total == 0:
        return "0"
    pct = count / total * 100
    return f"{count}\n({pct:.1f}%)"

# ===== Draw Venn Diagram =====
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.axis('off')

if n_groups == 2:
    # ---------- 2-set Venn ----------
    R = 0.5
    centers = [(-0.32, 0), (0.32, 0)]

    for i, (cx, cy) in enumerate(centers):
        c = Circle((cx, cy), R, fc=DARK2[i], ec='black', lw=2, alpha=0.5)
        ax.add_patch(c)
        # Group label outside circle
        lx = cx + (-0.35 if i == 0 else 0.35)
        ax.text(lx, cy + R * 0.55, groups[i],
                ha='center', va='center', fontsize=20, fontweight='bold')

    text_pos = {
        (0,):  (-0.5, 0),
        (1,):  (0.5, 0),
        (0, 1): (0, 0),
    }

elif n_groups == 3:
    # ---------- 3-set Venn ----------
    R = 0.5
    centers = [(0, 0.22), (-0.35, -0.18), (0.35, -0.18)]

    for i, (cx, cy) in enumerate(centers):
        c = Circle((cx, cy), R, fc=DARK2[i], ec='black', lw=2, alpha=0.5)
        ax.add_patch(c)

    # Group labels (outside circles, bold)
    ax.text(0,     0.80, groups[0], ha='center', va='center', fontsize=20, fontweight='bold')
    ax.text(-0.72, -0.58, groups[1], ha='center', va='center', fontsize=20, fontweight='bold')
    ax.text(0.72,  -0.58, groups[2], ha='center', va='center', fontsize=20, fontweight='bold')

    text_pos = {
        (0,):      (0,     0.52),
        (1,):      (-0.55, -0.32),
        (2,):      (0.55,  -0.32),
        (0, 1):    (-0.2,  0.08),
        (0, 2):    (0.2,   0.08),
        (1, 2):    (0,    -0.33),
        (0, 1, 2): (0,    -0.02),
    }

elif n_groups == 4:
    # ---------- 4-set Venn (Edwards ellipses) ----------
    ells = [
        (-0.15,  0.15, 0.8, 0.4, 45),
        ( 0.15,  0.15, 0.8, 0.4, 135),
        (-0.15, -0.15, 0.8, 0.4, 315),
        ( 0.15, -0.15, 0.8, 0.4, 225),
    ]
    for i, (cx, cy, w, h, ang) in enumerate(ells):
        e = Ellipse((cx, cy), w, h, angle=ang,
                    fc=DARK2[i], ec='black', lw=2, alpha=0.4)
        ax.add_patch(e)

    # Group labels at corners
    corners = [(-0.62, 0.65), (0.62, 0.65), (-0.62, -0.65), (0.62, -0.65)]
    for i in range(4):
        ax.text(corners[i][0], corners[i][1], groups[i],
                ha='center', va='center', fontsize=20, fontweight='bold')

    text_pos = {
        (0,):         (-0.45,  0.42),
        (1,):         ( 0.45,  0.42),
        (2,):         (-0.45, -0.42),
        (3,):         ( 0.45, -0.42),
        (0, 1):       ( 0,     0.38),
        (0, 2):       (-0.45,  0),
        (0, 3):       (-0.05,  0.18),
        (1, 2):       ( 0.05, -0.18),
        (1, 3):       ( 0.45,  0),
        (2, 3):       ( 0,    -0.38),
        (0, 1, 2):    (-0.18, -0.05),
        (0, 1, 3):    ( 0.18,  0.05),
        (0, 2, 3):    (-0.12, -0.22),
        (1, 2, 3):    ( 0.12, -0.22),
        (0, 1, 2, 3): ( 0,    -0.08),
    }

else:
    print(f"\nWarning: {n_groups} groups detected. Venn diagram supports 2-4 groups.")
    print("For more groups, consider using an UpSet plot.\n")
    print("Intersection summary:")
    for combo in sorted(regions.keys(), key=lambda x: (len(x), x)):
        names = " \u2229 ".join(groups[i] for i in combo)
        cnt = regions[combo]
        pct = cnt / total * 100 if total else 0
        print(f"  {names}: {cnt} ({pct:.1f}%)")
    plt.close()
    sys.exit(0)

# ===== Add Region Text Labels =====
for combo, pos in text_pos.items():
    count = regions.get(combo, 0)
    ax.text(pos[0], pos[1], fmt_label(count),
            ha='center', va='center', fontsize=13)

ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)

plt.tight_layout()
plt.savefig("venn_s.png", dpi=600, bbox_inches='tight')
plt.savefig("venn_s.svg", bbox_inches='tight')
print(f"\nVenn diagram saved: venn_s.png, venn_s.svg")
