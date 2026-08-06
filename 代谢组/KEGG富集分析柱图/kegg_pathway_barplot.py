#!/usr/bin/env python3
"""
KEGG Pathway Enrichment Bar Chart
R ggplot2 → pure matplotlib conversion

Input files:
  - kegg_enrichment.csv : KEGG enrichment result (ID, Description, GeneRatio, BgRatio, pvalue, EnrichmentRatio, geneID, Count, p.adjust)
  - kegg_info.csv       : KEGG pathway background (ID, Pathway.Identifier, Pathway, big.annotion, small.annotion)

Usage:
  python kegg_pathway_barplot.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import textwrap
from matplotlib.patches import Patch

# ═══════════════════════════════════════════════════════
#  Font & Global Style
# ═══════════════════════════════════════════════════════
mpl.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Times New Roman'],
    'axes.unicode_minus': False,
    'pdf.fonttype':       42,   # TrueType → editable in Illustrator
    'ps.fonttype':        42,
})

# ═══════════════════════════════════════════════════════
#  NPG (Nature Publishing Group) Color Palette
# ═══════════════════════════════════════════════════════
NPG_COLORS = [
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#5F5F5F', '#E18727', '#20854E', '#0072B5',
]

# ═══════════════════════════════════════════════════════
#  1. Read Data
# ═══════════════════════════════════════════════════════
kegg = pd.read_csv('kegg_enrichment.csv')
kegg_info = pd.read_csv('kegg_info.csv')

# ═══════════════════════════════════════════════════════
#  2. Filter pvalue < 0.05
# ═══════════════════════════════════════════════════════
kegg_p = kegg[kegg['pvalue'] < 0.05].copy()
print(f'KEGG pathways with pvalue < 0.05: {len(kegg_p)}')

# ═══════════════════════════════════════════════════════
#  3. Clean ID column — strip species prefix (e.g. "ko02010" → "02010")
# ═══════════════════════════════════════════════════════
kegg_p['ID'] = kegg_p['ID'].str.replace(r'^[a-z]+', '', regex=True)
kegg_p['ID'] = kegg_p['ID'].astype(int)

# ═══════════════════════════════════════════════════════
#  4. Extract columns from kegg_info (Class2 = small.annotion)
# ═══════════════════════════════════════════════════════
kegg_info_2class = kegg_info[['ID', 'Pathway', 'small.annotion']].copy()

# ═══════════════════════════════════════════════════════
#  5. Merge enrichment result with pathway background
# ═══════════════════════════════════════════════════════
kegg_p = kegg_p.merge(kegg_info_2class, on='ID', how='inner')

# ═══════════════════════════════════════════════════════
#  6. Extract & rename columns
# ═══════════════════════════════════════════════════════
kegg_2class = kegg_p[['Pathway', 'Count', 'small.annotion']].copy()
kegg_2class = kegg_2class.rename(columns={'small.annotion': 'Class2'})

# ═══════════════════════════════════════════════════════
#  7. Sort: Class2 total desc → Class2 name → Count desc
# ═══════════════════════════════════════════════════════
class2_total = kegg_2class.groupby('Class2')['Count'].transform('sum')
kegg_2class = kegg_2class.assign(Class2_Total=class2_total)
kegg_2class = kegg_2class.sort_values(
    by=['Class2_Total', 'Class2', 'Count'],
    ascending=[False, True, False]
).reset_index(drop=True)

print(f'Pathways to plot: {len(kegg_2class)}')
print(kegg_2class.groupby('Class2')['Count'].agg(['count', 'sum']).to_string())

# ═══════════════════════════════════════════════════════
#  8. Assign colors to Class2 categories
# ═══════════════════════════════════════════════════════
# Preserve the sort-order of categories (largest Class2 first)
class2_categories = list(kegg_2class['Class2'].unique())
color_map = {cat: NPG_COLORS[i % len(NPG_COLORS)] for i, cat in enumerate(class2_categories)}
kegg_2class['color'] = kegg_2class['Class2'].map(color_map)

# ═══════════════════════════════════════════════════════
#  9. Plot
# ═══════════════════════════════════════════════════════
n = len(kegg_2class)
fig_height = max(8, n * 0.42)
fig, ax = plt.subplots(figsize=(12, fig_height))

y_pos = np.arange(n)
counts = kegg_2class['Count'].values
bar_colors = kegg_2class['color'].values

# Horizontal bars
bars = ax.barh(y_pos, counts, height=0.7, color=bar_colors, alpha=0.9, edgecolor='none')

# ── Y-axis: pathway names (wrapped) ──
def wrap_label(text, width=35):
    """Wrap long pathway names to fit in the plot."""
    return '\n'.join(textwrap.wrap(str(text), width=width))

pathway_labels = [wrap_label(p, 35) for p in kegg_2class['Pathway']]
ax.set_yticks(y_pos)
ax.set_yticklabels(pathway_labels, fontsize=10)

# Invert so the top bar = first row (largest Class2)
ax.invert_yaxis()

# ── Count labels at bar ends ──
max_count = int(counts.max())
for i, (cnt, bar) in enumerate(zip(counts, bars)):
    ax.text(cnt + max_count * 0.02, i, str(int(cnt)),
            va='center', ha='left', fontsize=9,
            fontweight='bold', color='#333333')

# ── X-axis ──
ax.set_xlabel('Metabolite Count', fontsize=14, fontweight='bold')
ax.set_ylabel('KEGG Pathway', fontsize=14, fontweight='bold')
ax.set_xlim(0, max_count * 1.25)

# ── Title ──
ax.set_title('KEGG Pathway Classification', fontsize=18,
             fontweight='bold', pad=15)

# ── Legend ──
legend_elements = [
    Patch(facecolor=color_map[cat], alpha=0.9, label=cat)
    for cat in class2_categories
]
legend = ax.legend(
    handles=legend_elements,
    title='Pathway Category',
    loc='lower right',
    fontsize=9, title_fontsize=10,
    framealpha=0.95, edgecolor='#CCCCCC',
    bbox_to_anchor=(0.98, 0.02),
    bbox_transform=ax.transAxes,
)

# ── Spine cleanup ──
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', length=0)   # no tick marks on pathway names

# ── Add subtle horizontal grid lines ──
ax.xaxis.grid(True, linestyle='--', alpha=0.3, color='#888888')
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig('kegg_pathway_barplot.png', dpi=600, bbox_inches='tight')
plt.savefig('kegg_pathway_barplot.pdf', bbox_inches='tight')
plt.savefig('kegg_pathway_barplot.tiff', dpi=600, bbox_inches='tight',
            pil_kwargs={'compression': 'tiff_lzw'})
print('Saved: kegg_pathway_barplot.png / .pdf / .tiff')
