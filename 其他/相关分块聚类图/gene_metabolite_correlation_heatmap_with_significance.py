# -*- coding: utf-8 -*-
"""
Gene-metabolite correlation clustered heatmap with statistical significance.

Input:
1) gene_expression.xlsx
   - sheet: Expression
   - first column: gene name
   - remaining columns: samples

2) metabolite_abundance.xlsx
   - sheet: Metabolites
   - first column: metabolite name
   - remaining columns: samples

The two files do not need to have exactly the same sample columns.
Only common sample names are used, in the order found in the gene file.

For every gene-metabolite pair, correlation is calculated across matched samples.
Significance is calculated from the correlation test, followed by
Benjamini-Hochberg FDR correction across all gene-metabolite pairs.

Stars shown INSIDE heatmap cells:
    *   q < 0.05
    **  q < 0.01
    *** q < 0.001
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr

# ============================================================
# 0. Settings
# ============================================================

GENE_FILE = r"./gene_expression.xlsx"
GENE_SHEET = "Expression"

METABOLITE_FILE = r"./metabolite_abundance.xlsx"
METABOLITE_SHEET = "Metabolites"

# "pearson" = Pearson correlation
# "spearman" = Spearman rank correlation
CORRELATION_METHOD = "spearman"

# Minimum number of matched, non-missing samples required for a test.
MIN_VALID_SAMPLES = 3

# Number of visible row/column clusters used to create white gaps.
N_ROW_CLUSTERS = 4
N_COL_CLUSTERS = 4

CELL_SIZE = 1.0
GAP_SIZE = 0.65

FONT_SIZE = 13
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["font.size"] = FONT_SIZE
plt.rcParams["font.weight"] = "bold"


# ============================================================
# 1. Read the two data matrices
# ============================================================

gene_df = pd.read_excel(
    GENE_FILE,
    sheet_name=GENE_SHEET,
    index_col=0
)

metabolite_df = pd.read_excel(
    METABOLITE_FILE,
    sheet_name=METABOLITE_SHEET,
    index_col=0
)

gene_df.index = gene_df.index.astype(str)
gene_df.columns = gene_df.columns.astype(str)

metabolite_df.index = metabolite_df.index.astype(str)
metabolite_df.columns = metabolite_df.columns.astype(str)

# Convert all measurements to numeric; invalid entries become NaN.
gene_df = gene_df.apply(pd.to_numeric, errors="coerce")
metabolite_df = metabolite_df.apply(pd.to_numeric, errors="coerce")

# Match samples by column name.
common_samples = [
    sample for sample in gene_df.columns
    if sample in metabolite_df.columns
]

if len(common_samples) < MIN_VALID_SAMPLES:
    raise ValueError(
        f"Only {len(common_samples)} common samples were found. "
        f"At least {MIN_VALID_SAMPLES} matched samples are required."
    )

gene_df = gene_df.loc[:, common_samples]
metabolite_df = metabolite_df.loc[:, common_samples]

gene_labels = gene_df.index.tolist()
metabolite_labels = metabolite_df.index.tolist()

G = gene_df.to_numpy(dtype=float)
M = metabolite_df.to_numpy(dtype=float)

print("Gene file:", GENE_FILE)
print("Metabolite file:", METABOLITE_FILE)
print("Matched samples:", len(common_samples))
print("Sample names:", common_samples)
print("Genes:", len(gene_labels))
print("Metabolites:", len(metabolite_labels))
print("Correlation method:", CORRELATION_METHOD)


# ============================================================
# 2. Correlation and raw p-value for every gene-metabolite pair
# ============================================================

n_genes = G.shape[0]
n_metabolites = M.shape[0]

R = np.full((n_genes, n_metabolites), np.nan, dtype=float)
P = np.full((n_genes, n_metabolites), np.nan, dtype=float)
N = np.zeros((n_genes, n_metabolites), dtype=int)

method = CORRELATION_METHOD.lower()

if method not in {"pearson", "spearman"}:
    raise ValueError(
        'CORRELATION_METHOD must be either "pearson" or "spearman".'
    )

for i in range(n_genes):
    x_all = G[i, :]

    for j in range(n_metabolites):
        y_all = M[j, :]

        valid = np.isfinite(x_all) & np.isfinite(y_all)
        n_valid = int(valid.sum())
        N[i, j] = n_valid

        if n_valid < MIN_VALID_SAMPLES:
            continue

        x = x_all[valid]
        y = y_all[valid]

        # Correlation is undefined when either vector is constant.
        if np.all(x == x[0]) or np.all(y == y[0]):
            continue

        if method == "pearson":
            result = pearsonr(x, y)
        else:
            result = spearmanr(x, y)

        R[i, j] = result.statistic
        P[i, j] = result.pvalue


# ============================================================
# 3. Benjamini-Hochberg FDR correction
# ============================================================

def benjamini_hochberg(p_values):
    """BH-FDR correction while preserving NaN positions."""
    p_values = np.asarray(p_values, dtype=float)
    q_values = np.full(p_values.shape, np.nan, dtype=float)

    valid = np.isfinite(p_values)
    p = p_values[valid]

    if p.size == 0:
        return q_values

    order = np.argsort(p)
    ranked_p = p[order]
    m = ranked_p.size

    adjusted = ranked_p * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    restored = np.empty_like(adjusted)
    restored[order] = adjusted

    q_values[valid] = restored
    return q_values


Q = benjamini_hochberg(P.ravel()).reshape(P.shape)


def q_to_stars(q):
    if not np.isfinite(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


STARS = np.empty(Q.shape, dtype=object)

for i in range(n_genes):
    for j in range(n_metabolites):
        STARS[i, j] = q_to_stars(Q[i, j])


# ============================================================
# 4. Save complete statistical results
# ============================================================

records = []

for i, gene in enumerate(gene_labels):
    for j, metabolite in enumerate(metabolite_labels):
        records.append({
            "Gene": gene,
            "Metabolite": metabolite,
            "N": N[i, j],
            "Correlation_method": method,
            "Correlation_r": R[i, j],
            "Raw_p": P[i, j],
            "FDR_q": Q[i, j],
            "Significance": STARS[i, j]
        })

results_df = pd.DataFrame(records)
results_df.to_csv(
    "gene_metabolite_correlation_statistics.csv",
    index=False
)

# Also save matrices for easy checking.
pd.DataFrame(
    R,
    index=gene_labels,
    columns=metabolite_labels
).to_csv("gene_metabolite_correlation_matrix.csv")

pd.DataFrame(
    Q,
    index=gene_labels,
    columns=metabolite_labels
).to_csv("gene_metabolite_FDR_matrix.csv")


# ============================================================
# 5. Prepare matrix for clustering
# ============================================================

# Clustering cannot handle NaN. For clustering only, replace undefined
# correlations with 0. The original R values remain unchanged for plotting.
R_cluster = np.nan_to_num(R, nan=0.0)

if n_genes < 2 or n_metabolites < 2:
    raise ValueError(
        "At least 2 genes and 2 metabolites are required for clustering."
    )

row_Z = linkage(
    pdist(R_cluster, metric="euclidean"),
    method="average"
)

col_Z = linkage(
    pdist(R_cluster.T, metric="euclidean"),
    method="average"
)

row_order = dendrogram(row_Z, no_plot=True)["leaves"]
col_order = dendrogram(col_Z, no_plot=True)["leaves"]

row_cluster_count = min(N_ROW_CLUSTERS, n_genes)
col_cluster_count = min(N_COL_CLUSTERS, n_metabolites)

row_clusters = fcluster(
    row_Z,
    t=row_cluster_count,
    criterion="maxclust"
)[row_order]

col_clusters = fcluster(
    col_Z,
    t=col_cluster_count,
    criterion="maxclust"
)[col_order]

R_ordered = R[np.ix_(row_order, col_order)]
Q_ordered = Q[np.ix_(row_order, col_order)]
STARS_ordered = STARS[np.ix_(row_order, col_order)]

row_labels_ordered = [gene_labels[i] for i in row_order]
col_labels_ordered = [metabolite_labels[j] for j in col_order]


# ============================================================
# 6. Heatmap positions with cluster gaps
# ============================================================

def positions_with_gaps(cluster_ids):
    positions = []
    current = 0.0

    for k, cluster_id in enumerate(cluster_ids):
        if k > 0 and cluster_id != cluster_ids[k - 1]:
            current += GAP_SIZE

        positions.append(current)
        current += CELL_SIZE

    return np.asarray(positions)


def centers_to_bounds(centers):
    bounds = np.empty(len(centers) + 1, dtype=float)

    bounds[1:-1] = (centers[:-1] + centers[1:]) / 2
    bounds[0] = centers[0] - CELL_SIZE / 2
    bounds[-1] = centers[-1] + CELL_SIZE / 2

    return bounds


row_positions = positions_with_gaps(row_clusters)
col_positions = positions_with_gaps(col_clusters)

x_bounds = centers_to_bounds(col_positions)
y_bounds = centers_to_bounds(row_positions)


# ============================================================
# 7. Figure layout
# ============================================================

# Automatically increase figure size for larger matrices.
fig_width = max(9.6, 4.5 + 0.34 * n_metabolites)
fig_height = max(7.2, 3.5 + 0.28 * n_genes)

fig = plt.figure(
    figsize=(fig_width, fig_height),
    dpi=160
)

gs = gridspec.GridSpec(
    nrows=2,
    ncols=3,
    width_ratios=[1.25, 5.4, 0.75],
    height_ratios=[1.15, 5.2],
    wspace=0.05,
    hspace=0.03
)

ax_col_dendrogram = fig.add_subplot(gs[0, 1])
ax_row_dendrogram = fig.add_subplot(gs[1, 0])
ax_heatmap = fig.add_subplot(gs[1, 1])

ax_colorbar_holder = fig.add_subplot(gs[1, 2])
ax_colorbar_holder.axis("off")

holder_pos = ax_colorbar_holder.get_position()

COLORBAR_RIGHT_SHIFT = 1.20
COLORBAR_WIDTH_FRACTION = 0.22
COLORBAR_HEIGHT_FRACTION = 0.64

cbar_left = holder_pos.x0 + holder_pos.width * COLORBAR_RIGHT_SHIFT
cbar_bottom = holder_pos.y0 + holder_pos.height * 0.18
cbar_width = holder_pos.width * COLORBAR_WIDTH_FRACTION
cbar_height = holder_pos.height * COLORBAR_HEIGHT_FRACTION

ax_colorbar = fig.add_axes([
    cbar_left,
    cbar_bottom,
    cbar_width,
    cbar_height
])


# ============================================================
# 8. Dendrograms aligned to heatmap
# ============================================================

def remap_dendrogram_axis(coord_values, n_leaves, target_positions):
    default_leaf_positions = 5.0 + 10.0 * np.arange(n_leaves)
    return np.interp(
        coord_values,
        default_leaf_positions,
        target_positions
    )


col_ddata = dendrogram(
    col_Z,
    no_plot=True,
    color_threshold=0,
    above_threshold_color="#444444",
    link_color_func=lambda _: "#444444"
)

for icoord, dcoord in zip(
    col_ddata["icoord"],
    col_ddata["dcoord"]
):
    x = remap_dendrogram_axis(
        np.asarray(icoord),
        len(col_order),
        col_positions
    )
    y = np.asarray(dcoord)

    ax_col_dendrogram.plot(
        x,
        y,
        color="#444444",
        linewidth=0.8
    )

ax_col_dendrogram.set_xlim(
    x_bounds[0],
    x_bounds[-1]
)

ax_col_dendrogram.set_ylim(
    0,
    max(max(v) for v in col_ddata["dcoord"]) * 1.03
)

ax_col_dendrogram.set_xticks([])
ax_col_dendrogram.set_yticks([])

for spine in ax_col_dendrogram.spines.values():
    spine.set_visible(False)


row_ddata = dendrogram(
    row_Z,
    no_plot=True,
    color_threshold=0,
    above_threshold_color="#444444",
    link_color_func=lambda _: "#444444"
)

for icoord, dcoord in zip(
    row_ddata["icoord"],
    row_ddata["dcoord"]
):
    y = remap_dendrogram_axis(
        np.asarray(icoord),
        len(row_order),
        row_positions
    )
    x = -np.asarray(dcoord)

    ax_row_dendrogram.plot(
        x,
        y,
        color="#444444",
        linewidth=0.8
    )

row_max_dist = max(
    max(v) for v in row_ddata["dcoord"]
)

ax_row_dendrogram.set_xlim(
    -row_max_dist * 1.03,
    0
)

ax_row_dendrogram.set_ylim(
    y_bounds[-1],
    y_bounds[0]
)

ax_row_dendrogram.set_xticks([])
ax_row_dendrogram.set_yticks([])

for spine in ax_row_dendrogram.spines.values():
    spine.set_visible(False)


# ============================================================
# 9. Correlation heatmap
# ============================================================

# Mask undefined correlations so they appear white.
R_masked = np.ma.masked_invalid(R_ordered)

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("white")

mesh = ax_heatmap.pcolormesh(
    x_bounds,
    y_bounds,
    R_masked,
    cmap=cmap,
    vmin=-1.0,
    vmax=1.0,
    shading="flat",
    edgecolors=(1, 1, 1, 0.20),
    linewidth=0.35
)

ax_heatmap.set_xlim(
    x_bounds[0],
    x_bounds[-1]
)

ax_heatmap.set_ylim(
    y_bounds[-1],
    y_bounds[0]
)

ax_heatmap.set_xticks(col_positions)
ax_heatmap.set_xticklabels(
    col_labels_ordered,
    rotation=55,
    ha="right",
    fontsize=8
)

ax_heatmap.set_yticks(row_positions)
ax_heatmap.set_yticklabels(
    row_labels_ordered,
    fontsize=8
)

ax_heatmap.yaxis.tick_right()

ax_heatmap.tick_params(
    axis="both",
    length=0,
    pad=2
)

for spine in ax_heatmap.spines.values():
    spine.set_visible(False)


# ============================================================
# 10. Put significance stars INSIDE each heatmap cell
# ============================================================

for i, y in enumerate(row_positions):
    for j, x in enumerate(col_positions):
        star = STARS_ordered[i, j]

        if not star:
            continue

        r = R_ordered[i, j]

        # Use white text on stronger colors and black text near the center
        # so significance markers remain readable.
        text_color = "white" if np.isfinite(r) and abs(r) >= 0.55 else "black"

        ax_heatmap.text(
            x,
            y,
            star,
            ha="center",
            va="center",
            fontsize=7.5,
            fontweight="bold",
            color=text_color,
            zorder=10,
            clip_on=True
        )


# ============================================================
# 11. White gaps between clusters
# ============================================================

for k in range(1, len(col_positions)):
    if col_clusters[k] != col_clusters[k - 1]:
        x_mid = (
            col_positions[k - 1]
            + col_positions[k]
        ) / 2

        ax_heatmap.axvspan(
            x_mid - GAP_SIZE / 2,
            x_mid + GAP_SIZE / 2,
            color="white",
            zorder=5,
            linewidth=0
        )


for k in range(1, len(row_positions)):
    if row_clusters[k] != row_clusters[k - 1]:
        y_mid = (
            row_positions[k - 1]
            + row_positions[k]
        ) / 2

        ax_heatmap.axhspan(
            y_mid - GAP_SIZE / 2,
            y_mid + GAP_SIZE / 2,
            color="white",
            zorder=5,
            linewidth=0
        )


# ============================================================
# 12. Colorbar
# ============================================================

colorbar = fig.colorbar(
    mesh,
    cax=ax_colorbar
)

ticks = np.arange(
    -1.0,
    1.01,
    0.2
)

colorbar.set_ticks(ticks)

colorbar.set_ticklabels([
    "0" if abs(v) < 1e-12 else f"{v:.1f}"
    for v in ticks
])

colorbar.ax.tick_params(
    labelsize=8,
    length=0,
    pad=4
)

colorbar.set_label(
    f"{CORRELATION_METHOD.capitalize()} correlation",
    fontsize=9,
    fontweight="bold"
)

colorbar.outline.set_visible(False)


# ============================================================
# 13. Save
# ============================================================

plt.subplots_adjust(
    left=0.04,
    right=0.96,
    top=0.97,
    bottom=0.16
)

output_png = "gene_metabolite_correlation_heatmap.png"
output_svg = "gene_metabolite_correlation_heatmap.svg"

fig.savefig(
    output_png,
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)

fig.savefig(
    output_svg,
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)

plt.close(fig)

print("Saved:", output_png)
print("Saved:", output_svg)
print("Saved: gene_metabolite_correlation_statistics.csv")
print("Saved: gene_metabolite_correlation_matrix.csv")
print("Saved: gene_metabolite_FDR_matrix.csv")
print("Significance: * q<0.05, ** q<0.01, *** q<0.001")
