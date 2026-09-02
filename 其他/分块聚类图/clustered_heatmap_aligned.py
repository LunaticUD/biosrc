# -*- coding: utf-8 -*-
"""
Clustered heatmap similar to the reference image.

Changes in this version:
1. Removed the black rectangle borders around the white cluster gaps.
2. Moved the colorbar farther to the right.
3. Uses the Agg backend, so it works without Qt GUI support.
4. Saves the figure directly as a PNG file.

Dependencies:
    numpy
    scipy
    matplotlib
"""

import matplotlib
matplotlib.use("Agg")  # Avoid Qt backend/plugin errors

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

FONT_SIZE = 13
COLORBAR_FONT_SIZE = 10

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["font.size"] = FONT_SIZE
plt.rcParams["font.weight"] = "bold"


# ============================================================
# 1. Import expression matrix
# ============================================================

# Excel 文件路径
input_file = r"./example_expression_matrix.xlsx"
# 第一列作为基因名
df = pd.read_excel(
    input_file,
    sheet_name="Expression",
    index_col=0
)
# 行名 = gene
row_labels = df.index.astype(str).tolist()
# 列名 = sample
col_labels = df.columns.astype(str).tolist()
# 表达矩阵
X = df.to_numpy(dtype=float)
print("Input file:", input_file)
print("Matrix shape:", X.shape)
print("Genes:", len(row_labels))
print("Samples:", len(col_labels))


# ============================================================
# 2. Hierarchical clustering
# ============================================================
row_Z = linkage(pdist(X, metric="euclidean"), method="average")
col_Z = linkage(pdist(X.T, metric="euclidean"), method="average")

row_order = dendrogram(row_Z, no_plot=True)["leaves"]
col_order = dendrogram(col_Z, no_plot=True)["leaves"]

# Number of clusters used to create visible white gaps.
row_clusters = fcluster(row_Z, t=4, criterion="maxclust")[row_order]
col_clusters = fcluster(col_Z, t=4, criterion="maxclust")[col_order]

X_ordered = X[np.ix_(row_order, col_order)]
row_labels_ordered = [row_labels[i] for i in row_order]
col_labels_ordered = [col_labels[j] for j in col_order]


# ============================================================
# 3. Build cell positions with white gaps between clusters
# ============================================================
CELL_SIZE = 1.0
GAP_SIZE = 0.65


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
# 4. Figure layout
# ============================================================
fig = plt.figure(figsize=(9.6, 7.2), dpi=160)

# The third column is intentionally wider so the colorbar sits farther right.
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

# Use the third GridSpec cell as a holder, then place a narrower colorbar
# farther to the right within that area.
ax_colorbar_holder = fig.add_subplot(gs[1, 2])
ax_colorbar_holder.axis("off")

holder_pos = ax_colorbar_holder.get_position()

# Move the colorbar farther right.
# Increase COLORBAR_RIGHT_SHIFT if you want even more space.
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
# 5. Dendrograms aligned exactly to heatmap cell centers
# ============================================================
# scipy.dendrogram() uses fixed leaf coordinates 5, 15, 25, ... .
# The heatmap below uses custom coordinates with extra cluster gaps, so
# drawing dendrogram() directly on an Axes causes visible misalignment.
# We therefore obtain the dendrogram geometry with no_plot=True and remap
# every dendrogram leaf coordinate onto col_positions / row_positions.

def remap_dendrogram_axis(coord_values, n_leaves, target_positions):
    default_leaf_positions = 5.0 + 10.0 * np.arange(n_leaves)
    return np.interp(coord_values, default_leaf_positions, target_positions)


# ----- Column dendrogram -----
col_ddata = dendrogram(
    col_Z,
    no_plot=True,
    color_threshold=0,
    above_threshold_color="#444444",
    link_color_func=lambda _: "#444444"
)

for icoord, dcoord in zip(col_ddata["icoord"], col_ddata["dcoord"]):
    x = remap_dendrogram_axis(np.asarray(icoord), len(col_order), col_positions)
    y = np.asarray(dcoord)
    ax_col_dendrogram.plot(x, y, color="#444444", linewidth=0.8)

# Use exactly the same x range as the heatmap.
ax_col_dendrogram.set_xlim(x_bounds[0], x_bounds[-1])
ax_col_dendrogram.set_ylim(0, max(max(v) for v in col_ddata["dcoord"]) * 1.03)
ax_col_dendrogram.set_xticks([])
ax_col_dendrogram.set_yticks([])
for spine in ax_col_dendrogram.spines.values():
    spine.set_visible(False)


# ----- Row dendrogram -----
row_ddata = dendrogram(
    row_Z,
    no_plot=True,
    color_threshold=0,
    above_threshold_color="#444444",
    link_color_func=lambda _: "#444444"
)

for icoord, dcoord in zip(row_ddata["icoord"], row_ddata["dcoord"]):
    y = remap_dendrogram_axis(np.asarray(icoord), len(row_order), row_positions)
    x = -np.asarray(dcoord)
    ax_row_dendrogram.plot(x, y, color="#444444", linewidth=0.8)

# Use exactly the same y range and direction as the heatmap.
row_max_dist = max(max(v) for v in row_ddata["dcoord"])
ax_row_dendrogram.set_xlim(-row_max_dist * 1.03, 0)
ax_row_dendrogram.set_ylim(y_bounds[-1], y_bounds[0])
ax_row_dendrogram.set_xticks([])
ax_row_dendrogram.set_yticks([])
for spine in ax_row_dendrogram.spines.values():
    spine.set_visible(False)


# ============================================================
# 7. Heatmap
# ============================================================
mesh = ax_heatmap.pcolormesh(
    x_bounds,
    y_bounds,
    X_ordered,
    cmap="Spectral",
    vmin=-0.8,
    vmax=0.8,
    shading="flat",
    edgecolors=(1, 1, 1, 0.20),
    linewidth=0.35
)

ax_heatmap.set_xlim(x_bounds[0], x_bounds[-1])
ax_heatmap.set_ylim(y_bounds[-1], y_bounds[0])

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
# 8. White cluster gaps
#    IMPORTANT:
#    No black rectangles are drawn around these white gaps.
# ============================================================

# Vertical white gaps
for k in range(1, len(col_positions)):
    if col_clusters[k] != col_clusters[k - 1]:
        x_mid = (col_positions[k - 1] + col_positions[k]) / 2

        ax_heatmap.axvspan(
            x_mid - GAP_SIZE / 2,
            x_mid + GAP_SIZE / 2,
            color="white",
            zorder=5,
            linewidth=0
        )

# Horizontal white gaps
for k in range(1, len(row_positions)):
    if row_clusters[k] != row_clusters[k - 1]:
        y_mid = (row_positions[k - 1] + row_positions[k]) / 2

        ax_heatmap.axhspan(
            y_mid - GAP_SIZE / 2,
            y_mid + GAP_SIZE / 2,
            color="white",
            zorder=5,
            linewidth=0
        )


# ============================================================
# 9. Colorbar
# ============================================================
colorbar = fig.colorbar(mesh, cax=ax_colorbar)

ticks = np.arange(-0.8, 0.81, 0.2)
colorbar.set_ticks(ticks)

# Cleaner labels, e.g. "0" instead of "0.0"
colorbar.set_ticklabels([
    "0" if abs(v) < 1e-12 else f"{v:.1f}"
    for v in ticks
])

colorbar.ax.tick_params(
    labelsize=8,
    length=0,
    pad=4
)

colorbar.outline.set_visible(False)


# ============================================================
# 10. Save
# ============================================================
plt.subplots_adjust(
    left=0.04,
    right=0.96,
    top=0.97,
    bottom=0.13
)

output_file1 = "clustered_heatmap_no_gap_borders.png"
output_file2 = "clustered_heatmap_no_gap_borders.svg"
fig.savefig(
    output_file1,
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)
plt.close(fig)
print(f"Saved: {output_file1}")

fig.savefig(
    output_file2,
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)
plt.close(fig)
print(f"Saved: {output_file2}")
