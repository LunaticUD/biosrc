# ============================================================
# 02_Differential_Barplot.py
# Top Differential Metabolites
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib
matplotlib.use("Agg")

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 9,
    "axes.linewidth": 1,
    "pdf.fonttype": 42
})


# ============================================================
# 1. 数据
# ============================================================

df = pd.read_csv(
    "All_Metabolite_Statistics.csv",
    index_col=0
)


sig = df[
    (df["VIP"] > 1) &
    (df["log2FC"].abs() >= 1) &
    (df["P_value"] < 0.05)
].copy()


if sig.empty:
    raise ValueError(
        "No differential metabolites were found."
    )


# ============================================================
# 2. Top 15
# ============================================================

top_n = 15

sig["abs_log2FC"] = sig[
    "log2FC"
].abs()

top_df = (
    sig
    .sort_values(
        "abs_log2FC",
        ascending=False
    )
    .head(top_n)
    .sort_values("log2FC")
)


# ============================================================
# 3. 绘图
# ============================================================

fig, ax = plt.subplots(
    figsize=(6.5, 6.5)
)

colors = [
    "#D73027" if x > 0
    else "#1A9850"
    for x in top_df["log2FC"]
]


y = np.arange(
    len(top_df)
)


bars = ax.barh(
    y,
    top_df["log2FC"],
    color=colors,
    height=0.65,
    edgecolor="none"
)


# ============================================================
# 4. 数值
# ============================================================

max_abs = top_df[
    "log2FC"
].abs().max()

offset = max_abs * 0.035


for i, value in enumerate(
    top_df["log2FC"]
):

    if value >= 0:

        ax.text(
            value + offset,
            i,
            f"{value:.2f}",
            va="center",
            ha="left",
            fontsize=9
        )

    else:

        ax.text(
            value - offset,
            i,
            f"{value:.2f}",
            va="center",
            ha="right",
            fontsize=9
        )


# ============================================================
# 5. 格式
# ============================================================

ax.axvline(
    0,
    linewidth=0.8,
    color="black"
)

ax.set_yticks(y)

ax.set_yticklabels(
    top_df.index
)

ax.set_xlabel(
    r"$\log_{2}$(Fold Change)"
)

ax.set_ylabel("")

ax.set_title(
    "Top Differential Metabolites",
    fontweight="bold"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.30
)

ax.set_axisbelow(True)


# ============================================================
# 6. 图例
# ============================================================

legend_handles = [
    Patch(
        facecolor="#D73027",
        label="Up"
    ),
    Patch(
        facecolor="#1A9850",
        label="Down"
    )
]

ax.legend(
    handles=legend_handles,
    frameon=False,
    loc="lower right"
)


ax.set_xlim(
    -max_abs * 1.30,
    max_abs * 1.30
)


plt.tight_layout()

plt.savefig(
    "Top_Differential_Metabolites.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "Top_Differential_Metabolites.pdf",
    bbox_inches="tight"
)

plt.close()

print("Differential metabolite bar plot saved.")
