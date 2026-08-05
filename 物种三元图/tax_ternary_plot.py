"""
物种丰度三元图 (Ternary Plot)
R (ggtern) → Python (matplotlib) 纯手绘转换
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# ============================================================
# 全局字体
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",
})

# ============================================================
# 三元坐标 → 笛卡尔坐标 转换
# ============================================================
# 顶点定义: 左下=KCG, 右下=WBG, 顶=MCG
SQRT3_2 = np.sqrt(3) / 2

# 顶点笛卡尔坐标
V_KCG = np.array([0.0, 0.0])         # 左下 (1,0,0)
V_WBG  = np.array([1.0, 0.0])         # 右下 (0,1,0)
V_MCG  = np.array([0.5, SQRT3_2])     # 顶部 (0,0,1)


def ternary_to_cart(kcg, wbg, mcg):
    """三元坐标 (a,b,c) → 笛卡尔 (x,y)"""
    x = wbg + 0.5 * mcg
    y = SQRT3_2 * mcg
    return x, y


# ============================================================
# 数据读取与处理
# ============================================================

data = pd.read_csv("tax.txt", sep="\t", index_col=0)

# 按行和降序排列
data["sum"] = data.sum(axis=1)
data1 = data.sort_values("sum", ascending=False)
data1 = data1.drop(columns=["sum"])

# 分组计算平均值
# 分组列数自动检测：根据实际列数分配，也可手动指定
n_cols = len(data1.columns)
n_g1 = 8  # Group1 列数
n_g2 = 8  # Group2 列数
# Group3 为剩余列
data1["KCG"] = data1.iloc[:, 0:n_g1].mean(axis=1)
data1["WBG"] = data1.iloc[:, n_g1:n_g1+n_g2].mean(axis=1)
data1["MCG"] = data1.iloc[:, n_g1+n_g2:n_cols].mean(axis=1)

# 提取平均值列
data2 = data1[["KCG", "WBG", "MCG"]]

# 取丰度前10的物种，计算相对丰度
top10_raw = data2.iloc[:10]
top10 = top10_raw / top10_raw.sum(axis=0)

# 计算三组平均丰度（用于点大小排序）
top10["Abundance"] = top10[["KCG", "WBG", "MCG"]].mean(axis=1)

# 按平均丰度降序排列
top10 = top10.sort_values("Abundance", ascending=False)

# 显示名：下划线→空格，物种斜体
display_names = []
for name in top10.index.tolist():
    if name == "Others":
        display_names.append("Others")
    else:
        display_names.append(name.replace("_", " "))

# ============================================================
# 配色方案
# ============================================================

palette = [
    "#9999CC", "#CC6666", "#CC79A7", "#D55E00", "#0072B2",
    "#F0E442", "#009E73", "#56B4E9", "#E69F00", "#B2182B",
]

n_taxa = len(top10)
colors = palette[:n_taxa]

# ============================================================
# 绘图辅助函数
# ============================================================


def draw_ternary_frame(ax, n_grid=5):
    """画三元图边框和网格"""
    # 外边框
    triangle = MplPolygon(
        [V_KCG, V_WBG, V_MCG],
        closed=True, fill=False, edgecolor="black", linewidth=1.2,
    )
    ax.add_patch(triangle)

    # 网格线
    for i in range(1, n_grid):
        frac = i / n_grid

        # MCG = const 的平行线（平行于底边）
        p1 = ternary_to_cart(1 - frac, 0, frac)
        p2 = ternary_to_cart(0, 1 - frac, frac)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="grey", linewidth=0.4, alpha=0.5)

        # KCG = const 的平行线（平行于右边）
        p1 = ternary_to_cart(frac, 1 - frac, 0)
        p2 = ternary_to_cart(frac, 0, 1 - frac)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="grey", linewidth=0.4, alpha=0.5)

        # WBG = const 的平行线（平行于左边）
        p1 = ternary_to_cart(1 - frac, frac, 0)
        p2 = ternary_to_cart(0, frac, 1 - frac)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="grey", linewidth=0.4, alpha=0.5)

    # 刻度标签
    for i in range(0, n_grid + 1):
        frac = i / n_grid

        # KCG 轴（左边，从下往上）
        px, py = ternary_to_cart(1 - frac, 0, frac)
        ax.text(px - 0.04, py, f"{frac:.1f}", ha="right", va="center", fontsize=11, color="grey")

        # WBG 轴（底边，从左往右）
        px, py = ternary_to_cart(1 - frac, frac, 0)
        ax.text(px, py - 0.04, f"{frac:.1f}", ha="center", va="top", fontsize=11, color="grey")

        # MCG 轴（右边，从上往下）
        px, py = ternary_to_cart(0, 1 - frac, frac)
        ax.text(px + 0.04, py, f"{frac:.1f}", ha="left", va="center", fontsize=11, color="grey")


def draw_ternary_labels(ax):
    """画顶点标签"""
    # 顶点名称
    ax.text(V_KCG[0] - 0.06, V_KCG[1] - 0.04, "KCG", fontsize=12,
            ha="right", va="top", color="black", fontweight="bold")
    ax.text(V_WBG[0] + 0.06, V_WBG[1] - 0.04, "WBG", fontsize=12,
            ha="left", va="top", color="black", fontweight="bold")
    ax.text(V_MCG[0], V_MCG[1] + 0.04, "MCG", fontsize=12,
            ha="center", va="bottom", color="black", fontweight="bold")


# ============================================================
# 图2: 三元图 + 右侧图例
# ============================================================

fig2, ax2 = plt.subplots(figsize=(8.4, 7.4))

draw_ternary_frame(ax2)
draw_ternary_labels(ax2)

# 画散点 + 收集图例句柄
legend_handles = []
for idx, (taxon, row) in enumerate(top10.iterrows()):
    x, y = ternary_to_cart(row["KCG"], row["WBG"], row["MCG"])
    size = row["Abundance"] * 3000
    ax2.scatter(x, y, s=size, c=colors[idx], alpha=0.8,
                edgecolors="white", linewidths=0.5, zorder=5)

    # 图例：统一大小
    label = display_names[idx]
    label_fmt = label if taxon == "Others" else rf"$\mathit{{{label}}}$"
    h = ax2.scatter([], [], s=80, c=colors[idx], label=label_fmt,
                    edgecolors="white", linewidths=0.5)
    legend_handles.append(h)

ax2.set_xlim(-0.15, 1.15)
ax2.set_ylim(-0.15, SQRT3_2 + 0.1)
ax2.set_aspect("equal")
ax2.axis("off")

# 图例放在图内右侧
leg = ax2.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(0.78, 0.68),
    ncol=1,
    frameon=False,
    fontsize=13,
    labelspacing=0.8,
    borderpad=1,
    handletextpad=0.5,
    title="Taxonomy",
    title_fontsize=23,
)
# 图例标题颜色
leg.get_title().set_color("black")

fig2.tight_layout()
fig2.savefig("sanyuan.tif", dpi=600, bbox_inches="tight")
fig2.savefig("sanyuan.pdf", bbox_inches="tight")
print("图2已保存: sanyuan.tif, sanyuan.pdf")

plt.show()
