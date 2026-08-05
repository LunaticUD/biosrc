"""
物种丰度热图 (pheatmap 风格)
R (pheatmap) → Python (scipy + seaborn + matplotlib) 转换
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

# ============================================================
# 全局字体
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",
})

# ============================================================
# 数据读取
# ============================================================

dat = pd.read_csv("tax.txt", sep="\t", index_col=0)
ann_col = pd.read_csv("group.txt", sep="\t", index_col=0)

# 分组因子排序（R: levels=c("KC","WB","MC")）
group_order = {"KC": 0, "WB": 1, "MC": 2}
ann_col["group"] = pd.Categorical(
    ann_col["group"], categories=["KC", "WB", "MC"], ordered=True
)

# 确保样本名对齐
common_samples = dat.columns.intersection(ann_col.index)
dat = dat[common_samples]
ann_col = ann_col.loc[common_samples]

# ============================================================
# 行标准化 (scale="row")
# ============================================================

dat_scaled = dat.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# ============================================================
# 聚类
# ============================================================

# 行聚类
row_linkage = linkage(
    pdist(dat_scaled.values, metric="euclidean"), method="complete"
)
# 列聚类
col_linkage = linkage(
    pdist(dat_scaled.values.T, metric="euclidean"), method="complete"
)

# 聚类排序
row_dendro = dendrogram(row_linkage, no_plot=True)
col_dendro = dendrogram(col_linkage, no_plot=True)
row_order = row_dendro["leaves"]
col_order = col_dendro["leaves"]

# 剪切树 (cutree_rows=2, cutree_cols=3)
row_clusters = fcluster(row_linkage, t=2, criterion="maxclust")
col_clusters = fcluster(col_linkage, t=3, criterion="maxclust")

# 按聚类排序重排数据和注释
dat_sorted = dat_scaled.iloc[row_order, col_order]
ann_sorted = ann_col.iloc[col_order]

# ============================================================
# 颜色配置
# ============================================================

# 热图颜色渐变: navy → white → firebrick3
palette = mcolors.LinearSegmentedColormap.from_list(
    "navy_white_firebrick",
    ["navy", "white", "#CD2626"],
    N=100,
)

# 分组注释颜色
group_colors = {"KC": "#7570b3", "WB": "#1b9e77", "MC": "#d95f02"}

# ============================================================
# 绘图
# ============================================================

n_rows, n_cols = dat_sorted.shape
fig = plt.figure(figsize=(max(8, n_cols * 0.4 + 3), max(6, n_rows * 0.06 + 3)))

# 布局: 左=行树状图, 上=列树状图, 顶=分组注释, 中=热图, 右=颜色条
ax_row_dend = fig.add_axes([0.05, 0.12, 0.12, 0.72])   # 行树状图
ax_col_dend = fig.add_axes([0.22, 0.86, 0.65, 0.10])    # 列树状图
ax_ann      = fig.add_axes([0.22, 0.82, 0.65, 0.03])    # 分组注释条
ax_heatmap  = fig.add_axes([0.22, 0.12, 0.65, 0.72])    # 热图
ax_cbar     = fig.add_axes([0.90, 0.12, 0.02, 0.72])    # 颜色条

# --- 行树状图 ---
dendrogram(
    row_linkage,
    orientation="left",
    ax=ax_row_dend,
    color_threshold=0,
    above_threshold_color="black",
    no_labels=True,
)
ax_row_dend.invert_yaxis()
ax_row_dend.set_xticks([])
ax_row_dend.set_yticks([])
ax_row_dend.spines["top"].set_visible(False)
ax_row_dend.spines["right"].set_visible(False)
ax_row_dend.spines["bottom"].set_visible(False)
ax_row_dend.spines["left"].set_visible(False)

# --- 列树状图 ---
dendrogram(
    col_linkage,
    ax=ax_col_dend,
    color_threshold=0,
    above_threshold_color="black",
    no_labels=True,
)
ax_col_dend.set_xticks([])
ax_col_dend.set_yticks([])
ax_col_dend.spines["top"].set_visible(False)
ax_col_dend.spines["right"].set_visible(False)
ax_col_dend.spines["bottom"].set_visible(False)
ax_col_dend.spines["left"].set_visible(False)

# --- 分组注释条 ---
for i, sample in enumerate(dat_sorted.columns):
    grp = ann_sorted.loc[sample, "group"]
    ax_ann.bar(i, 1, width=1, color=group_colors[grp], edgecolor="none", align="edge")
ax_ann.set_xlim(0, n_cols)
ax_ann.set_ylim(0, 1)
ax_ann.set_xticks([])
ax_ann.set_yticks([])
for spine in ax_ann.spines.values():
    spine.set_visible(False)
ax_ann.patch.set_alpha(0)

# --- 热图 ---
vmax = max(abs(dat_sorted.values.min()), abs(dat_sorted.values.max()))
im = ax_heatmap.imshow(
    dat_sorted.values,
    aspect="auto",
    cmap=palette,
    vmin=-vmax,
    vmax=vmax,
    interpolation="nearest",
)

# 列名（旋转90°）
col_labels = dat_sorted.columns.tolist()
ax_heatmap.set_xticks(range(n_cols))
ax_heatmap.set_xticklabels(col_labels, rotation=90, fontsize=10, ha="center")

# 行名（R: show_rownames=FALSE）
ax_heatmap.set_yticks([])

# 网格线（白色细线分隔单元格）
for i in range(n_rows + 1):
    ax_heatmap.axhline(i - 0.5, color="white", linewidth=0.3)
for j in range(n_cols + 1):
    ax_heatmap.axvline(j - 0.5, color="white", linewidth=0.3)

# 行切割线 (cutree_rows=2)
row_cluster_boundaries = []
prev_c = row_clusters[row_order[0]]
for idx, i in enumerate(row_order):
    c = row_clusters[i]
    if c != prev_c:
        row_cluster_boundaries.append(idx)
    prev_c = c
for b in row_cluster_boundaries:
    ax_heatmap.axhline(b - 0.5, color="black", linewidth=1.5)

# 列切割线 (cutree_cols=3)
col_cluster_boundaries = []
prev_c = col_clusters[col_order[0]]
for idx, j in enumerate(col_order):
    c = col_clusters[j]
    if c != prev_c:
        col_cluster_boundaries.append(idx)
    prev_c = c
for b in col_cluster_boundaries:
    ax_heatmap.axvline(b - 0.5, color="black", linewidth=1.5)

# 去掉边框刻度
ax_heatmap.tick_params(axis="y", length=0)

# --- 颜色条 ---
cbar = fig.colorbar(im, cax=ax_cbar)
cbar.ax.tick_params(labelsize=9)

# --- 分组图例 ---
legend_patches = [Patch(facecolor=group_colors[g], label=g) for g in ["KC", "WB", "MC"]]
ax_ann.legend(
    handles=legend_patches,
    loc="upper left",
    bbox_to_anchor=(1.01, 1.2),
    frameon=False,
    fontsize=11,
    title="Group",
    title_fontsize=12,
)

# ============================================================
# 保存
# ============================================================

fig.savefig("tax_heatmap.png", dpi=600, bbox_inches="tight")
print("热图已保存: tax_heatmap.png")
fig.savefig("tax_heatmap.svg", dpi=600, bbox_inches="tight")
print("热图已保存: tax_heatmap.svg")
plt.show()
