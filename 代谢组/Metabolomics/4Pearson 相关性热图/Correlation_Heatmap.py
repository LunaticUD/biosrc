# ============================================================
# 03_Correlation_Heatmap.py
# Pearson Correlation Heatmap
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use("Agg")

plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.titlesize": 13,
    "pdf.fonttype": 42
})


# ============================================================
# 1. 数据
# ============================================================

sig = pd.read_csv(
    "Differential_Metabolites.csv",
    index_col=0
)

expr_df = pd.read_csv(
    "Expression_Matrix.csv",
    index_col=0
)


common_ids = sig.index.intersection(
    expr_df.index
)


if len(common_ids) < 2:
    raise ValueError(
        "At least two metabolites are required."
    )


sig_expr = expr_df.loc[
    common_ids
]


# ============================================================
# 2. Pearson相关矩阵
# ============================================================

corr = sig_expr.T.corr(
    method="pearson"
)


# ============================================================
# 3. 下三角Mask
# ============================================================

mask = np.triu(
    np.ones_like(
        corr,
        dtype=bool
    ),
    k=1
)


# ============================================================
# 4. 绘图
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 7)
)


sns.heatmap(
    corr,
    mask=mask,
    cmap="RdBu_r",
    center=0,
    vmin=-1,
    vmax=1,
    square=True,
    linewidths=0.35,
    linecolor="white",
    annot=False,
    cbar_kws={
        "label": "Pearson correlation",
        "shrink": 0.82,
        "pad": 0.02
    },
    ax=ax
)


ax.set_xlabel("")
ax.set_ylabel("")

ax.set_title(
    "Pearson Correlation of Differential Metabolites",
    fontsize=13,
    fontweight="bold",
    pad=15
)

ax.tick_params(
    axis="x",
    rotation=90,
    labelsize=8
)

ax.tick_params(
    axis="y",
    rotation=0,
    labelsize=8
)


plt.tight_layout()

plt.savefig(
    "Correlation_Heatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "Correlation_Heatmap.pdf",
    bbox_inches="tight"
)

plt.close()

print(
    f"Correlation heatmap saved: "
    f"{len(common_ids)} metabolites."
)
