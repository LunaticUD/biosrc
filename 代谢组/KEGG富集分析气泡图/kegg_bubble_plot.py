#!/usr/bin/env python3
"""
KEGG Pathway Enrichment Bubble Plot
R ggplot2 + KeggMetaboliteEnrichment → pure matplotlib conversion

Input file:
  - kegg_enrichment.csv : KEGG enrichment result
    (ID, Description, GeneRatio, BgRatio, pvalue, EnrichmentRatio, geneID, Count, p.adjust)

Usage:
  python kegg_bubble_plot.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import textwrap

# ═══════════════════════════════════════════════════════
#  Font & Global Style
# ═══════════════════════════════════════════════════════
mpl.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Times New Roman'],
    'axes.unicode_minus': False,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

# ═══════════════════════════════════════════════════════
#  1. Read Data
# ═══════════════════════════════════════════════════════
kegg = pd.read_csv('kegg_enrichment.csv')

# ═══════════════════════════════════════════════════════
#  2. Filter pvalue < 0.05, sort by p.adjust, take top N
# ═══════════════════════════════════════════════════════
kegg_sig = kegg[kegg['pvalue'] < 0.05].copy()
kegg_sig = kegg_sig.sort_values('p.adjust').reset_index(drop=True)

# Top 10 pathways (change this number as needed)
TOP_N = 10
top_data = kegg_sig.head(TOP_N).copy()

print(f'Significant pathways (pvalue < 0.05): {len(kegg_sig)}')
print(f'Plotting top {len(top_data)} pathways')

# ═══════════════════════════════════════════════════════
#  3. Clean Description — remove species suffix
# ═══════════════════════════════════════════════════════
top_data['Description'] = top_data['Description'].str.replace(
    r'\s*-\s*.*$', '', regex=True
)

# ═══════════════════════════════════════════════════════
#  4. Compute GeneRatio_value and -log10(p.adjust)
# ═══════════════════════════════════════════════════════
# GeneRatio is like "25/159" → compute 25/159 = 0.1572
top_data['GeneRatio_value'] = top_data['GeneRatio'].apply(
    lambda x: float(x.split('/')[0]) / float(x.split('/')[1])
)

top_data['neg_log_padj'] = -np.log10(top_data['p.adjust'])

# ═══════════════════════════════════════════════════════
#  5. Sort: most significant at top of Y-axis
# ═══════════════════════════════════════════════════════
top_data = top_data.sort_values('p.adjust', ascending=False).reset_index(drop=True)

# Wrap long pathway names
def wrap_label(text, width=35):
    return '\n'.join(textwrap.wrap(str(text), width=width))

top_data['Description_wrap'] = top_data['Description'].apply(lambda x: wrap_label(x, 35))

# ═══════════════════════════════════════════════════════
#  6. Plot Bubble Chart
# ═══════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 8))

x = top_data['GeneRatio_value'].values
y = np.arange(len(top_data))
sizes = top_data['Count'].values
colors = top_data['neg_log_padj'].values

# ── Scatter (bubble) ──
scatter = ax.scatter(
    x, y,
    s=sizes * 30,          # scale count → marker area (adjust multiplier for visual)
    c=colors,               # color = -log10(p.adjust)
    cmap='RdBu_r',          # blue(low) → red(high), matches R's low="#4C72B0" high="#C44E52"
    alpha=0.8,
    edgecolors='white',
    linewidths=0.5,
    vmin=colors.min(),
    vmax=colors.max(),
)

# ── Y-axis: pathway names ──
ax.set_yticks(y)
ax.set_yticklabels(top_data['Description_wrap'].values, fontsize=12)

# ── X-axis ──
ax.set_xlabel('Rich Factor', fontsize=16, fontweight='bold')
ax.set_ylabel('Pathway Name', fontsize=16, fontweight='bold')

# ── Title ──
ax.set_title('KEGG Pathway Enrichment Analysis', fontsize=18,
             fontweight='bold', pad=15)

# ── Colorbar ──
cbar = fig.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('-log10(p.adjust)', fontsize=13, fontweight='bold')
cbar.ax.tick_params(labelsize=11)

# ── Size legend ──
# Create a separate legend for bubble sizes
count_min = int(sizes.min())
count_max = int(sizes.max())
count_mid = (count_min + count_max) // 2
legend_sizes = [count_min, count_mid, count_max]
legend_bubbles = [
    ax.scatter([], [], s=s * 30, c='gray', alpha=0.6, edgecolors='white',
               linewidths=0.5, label=str(s))
    for s in legend_sizes
]
size_legend = ax.legend(
    handles=legend_bubbles,
    title='Metabolite Number',
    loc='lower right',
    fontsize=10, title_fontsize=11,
    framealpha=0.9, edgecolor='#CCCCCC',
    labelspacing=1.2,
    borderpad=0.8,
)

# Add the size legend back (colorbar may have removed it)
ax.add_artist(size_legend)

# ── Spine & grid ──
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', linestyle='--', alpha=0.3, color='#888888')
ax.set_axisbelow(True)
ax.tick_params(axis='x', labelsize=12)

plt.tight_layout()
plt.savefig('kegg_bubble_plot.png', dpi=600, bbox_inches='tight',
            facecolor='white')
plt.savefig('kegg_bubble_plot.pdf', bbox_inches='tight',
            facecolor='white')
plt.savefig('kegg_bubble_plot.tiff', dpi=600, bbox_inches='tight',
            facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
print('Saved: kegg_bubble_plot.png / .pdf / .tiff')
