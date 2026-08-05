"""
物种丰度堆叠柱状图 & 冲积图
R (ggplot2 + ggalluvial) → Python (matplotlib) 转换
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# 全局字体设为 Times New Roman
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",  # 数学公式也用类似 Times 的字体
})

# ============================================================
# 数据读取与处理
# ============================================================

data = pd.read_csv("tax.txt", sep="\t", index_col=0)

# 按行和降序排列
data["sum"] = data.sum(axis=1)
data1 = data.sort_values("sum", ascending=False)
data1 = data1.drop(columns=["sum"])

# 分组计算平均值
# R 1-indexed: Group1=cols[1:8], Group2=cols[9:16], Group3=cols[17:23]
# Python 0-indexed:
data1["Group1"] = data1.iloc[:, 0:8].mean(axis=1)
data1["Group2"] = data1.iloc[:, 8:16].mean(axis=1)
data1["Group3"] = data1.iloc[:, 16:23].mean(axis=1)

# 提取平均值列
data2 = data1[["Group1", "Group2", "Group3"]]

# 取丰度前10的物种，计算相对丰度
top10 = data2.iloc[:10] / data2.sum(axis=0)
others = 1 - top10.sum(axis=0)
top10.loc["Others"] = others

# 检查数据
print("行名:", top10.index.tolist())
print("列名:", top10.columns.tolist())
print("各列总和:", top10.sum(axis=0).values)

# ============================================================
# 配色方案 (RColorBrewer Set3)
# ============================================================

set3_colors = [
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#CCEBC5", "#FFED6F",
]
n_taxa = len(top10)
colors = set3_colors[:n_taxa]
# 显示名：下划线→空格，物种学名斜体，Others不斜体
display_names = {}
for name in top10.index.tolist():
    if name == "Others":
        display_names[name] = "Others"
    else:
        display_names[name] = name.replace("_", " ")

color_map = dict(zip(top10.index.tolist(), colors))

# R: factor(levels = rev(rownames)) → Others在底部，丰度最高的在顶部
stack_order = list(reversed(top10.index.tolist()))

groups = top10.columns.tolist()
x_pos = np.arange(len(groups))


def apply_theme(ax):
    """统一绘图样式（对应R代码的theme设置）"""
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlabel("Samples", fontsize=13)
    ax.set_ylabel("Relative Abundance (%)", fontsize=13)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=12)
    ax.grid(False)


# ============================================================
# 图1: 堆叠柱状图
# ============================================================

fig1, ax1 = plt.subplots(figsize=(6, 4))

bottom = np.zeros(len(groups))
for taxon in stack_order:
    values = top10.loc[taxon].values * 100
    label = display_names[taxon]
    ax1.bar(
        x_pos, values, bottom=bottom,
        label=label if taxon == "Others" else rf"$\mathit{{{label}}}$",
        width=0.6, color=color_map[taxon], edgecolor="white", linewidth=0.3,
    )
    bottom += values

apply_theme(ax1)

# 图例：反转使丰度最高的在顶部
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(
    handles[::-1], labels[::-1],
    loc="center left", bbox_to_anchor=(1, 0.5),
    frameon=False, fontsize=11,
)

fig1.tight_layout()
fig1.savefig("stacked_barplot.tif", dpi=600, bbox_inches="tight")
print("图1已保存: stacked_barplot.tif")

# ============================================================
# 图2: 堆叠柱状图 + 冲积连线 (alluvial)
# ============================================================

fig2, ax2 = plt.subplots(figsize=(6, 4))

bar_width = 0.6
half_w = bar_width / 2

# 计算每个分类在每个组中的位置 (bottom, top)
positions = {}
for i, group in enumerate(groups):
    cumsum = 0
    for taxon in stack_order:
        val = top10.loc[taxon, group] * 100
        positions[(taxon, i)] = (cumsum, cumsum + val)
        cumsum += val

# 画柱子 (stratum)
bottom = np.zeros(len(groups))
for taxon in stack_order:
    values = top10.loc[taxon].values * 100
    label = display_names[taxon]
    ax2.bar(
        x_pos, values, bottom=bottom,
        label=label if taxon == "Others" else rf"$\mathit{{{label}}}$",
        width=bar_width, color=color_map[taxon],
        edgecolor="white", linewidth=0.5,
    )
    bottom += values

# 画冲积连线 (alluvium) — curve_type="linear"
for taxon in stack_order:
    for i in range(len(groups) - 1):
        left_right_edge = x_pos[i] + half_w      # 左柱右边缘
        right_left_edge = x_pos[i + 1] - half_w   # 右柱左边缘

        l_bottom, l_top = positions[(taxon, i)]
        r_bottom, r_top = positions[(taxon, i + 1)]

        verts = [
            (left_right_edge, l_bottom),
            (right_left_edge, r_bottom),
            (right_left_edge, r_top),
            (left_right_edge, l_top),
        ]
        poly = Polygon(
            verts, alpha=0.4, facecolor=color_map[taxon],
            edgecolor="white", linewidth=0.5,
        )
        ax2.add_patch(poly)

apply_theme(ax2)
ax2.set_xlim(x_pos[0] - 0.5, x_pos[-1] + 0.5)

# 图例
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(
    handles[::-1], labels[::-1],
    loc="center left", bbox_to_anchor=(1, 0.5),
    frameon=False, fontsize=11,
)

fig2.tight_layout()
fig2.savefig("stacked_barplot_2.tif", dpi=600, bbox_inches="tight")
print("图2已保存: stacked_barplot_2.tif")

plt.show()
