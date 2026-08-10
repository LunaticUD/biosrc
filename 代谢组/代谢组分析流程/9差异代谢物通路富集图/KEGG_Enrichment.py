# ============================================================
# 09_KEGG_Enrichment.py
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")

# ============================================================
# 1. SCI style
# ============================================================

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10


# ============================================================
# 2. Read differential metabolites
# ============================================================

sig = pd.read_csv(
    "Differential_Metabolites.csv",
    index_col=0
)


# ============================================================
# 3. Read KEGG annotation
# ============================================================

annotation_df = pd.read_csv(
    "Metabolite_KEGG_Annotation.csv"
)


# 只保留差异代谢物
sig_annotation = annotation_df[
    annotation_df["Metabolite"].isin(
        sig.index
    )
].copy()


if sig_annotation.empty:

    raise ValueError(
        "No differential metabolites matched KEGG annotation."
    )


# ============================================================
# 4. Pathway differential counts
# ============================================================

pathway_hits = (
    sig_annotation
    .groupby("Pathway")
    ["Metabolite"]
    .nunique()
)


# ============================================================
# 5. Background file
# ============================================================

background_file = (
    "KEGG_Pathway_Background.csv"
)


if not os.path.exists(
    background_file
):

    print(
        "Generating example KEGG background..."
    )

    np.random.seed(42)

    pathways = (
        annotation_df["Pathway"]
        .drop_duplicates()
        .tolist()
    )

    background_df = pd.DataFrame({

        "Pathway":
            pathways,

        "Background_Count":
            np.random.randint(
                30,
                120,
                size=len(pathways)
            )
    })


    # 确保背景数 > 差异命中数
    for pathway in pathways:

        hit = pathway_hits.get(
            pathway,
            0
        )

        idx = (
            background_df["Pathway"]
            == pathway
        )

        background_df.loc[
            idx,
            "Background_Count"
        ] = np.maximum(
            background_df.loc[
                idx,
                "Background_Count"
            ],
            hit + 5
        )


    background_df.to_csv(
        background_file,
        index=False,
        encoding="utf-8-sig"
    )

else:

    background_df = pd.read_csv(
        background_file
    )


# ============================================================
# 6. Global background parameters
#
# 示例版：
# M = 所有 pathway background count 总和
# N = 差异代谢物数量
#
# 正式分析时更推荐使用真实可注释代谢物全集
# ============================================================

M = int(
    background_df[
        "Background_Count"
    ].sum()
)

N = len(
    sig.index.unique()
)


# ============================================================
# 7. Enrichment
# ============================================================

results = []


for _, row in background_df.iterrows():

    pathway = row["Pathway"]

    K = int(
        row["Background_Count"]
    )

    k = int(
        pathway_hits.get(
            pathway,
            0
        )
    )


    if k == 0:
        continue


    # Rich factor
    rich_factor = (
        k / K
    )


    # Hypergeometric upper-tail P
    p_value = hypergeom.sf(
        k - 1,
        M,
        K,
        N
    )


    results.append({

        "Pathway":
            pathway,

        "Count":
            k,

        "Background_Count":
            K,

        "RichFactor":
            rich_factor,

        "P_value":
            p_value
    })


kegg_df = pd.DataFrame(
    results
)


if kegg_df.empty:

    raise ValueError(
        "No KEGG pathways available for enrichment."
    )


# ============================================================
# 8. FDR correction
# ============================================================

kegg_df["FDR"] = multipletests(
    kegg_df["P_value"],
    method="fdr_bh"
)[1]


kegg_df["neg_log10_FDR"] = (
    -np.log10(
        kegg_df["FDR"].clip(
            lower=1e-300
        )
    )
)


# ============================================================
# 9. Save enrichment table
# ============================================================

kegg_df.to_csv(
    "KEGG_Enrichment_Result.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 10. Select top pathways
# ============================================================

top_n = 15


plot_data = (
    kegg_df
    .sort_values(
        [
            "FDR",
            "RichFactor"
        ],
        ascending=[
            True,
            False
        ]
    )
    .head(top_n)
    .sort_values(
        "RichFactor",
        ascending=True
    )
)


# ============================================================
# 11. Bubble sizes
# ============================================================

count_min = (
    plot_data["Count"].min()
)

count_max = (
    plot_data["Count"].max()
)


if count_max == count_min:

    point_sizes = np.repeat(
        180,
        len(plot_data)
    )

else:

    point_sizes = (
        80
        +
        (
            plot_data["Count"]
            - count_min
        )
        /
        (
            count_max
            - count_min
        )
        * 420
    )


# ============================================================
# 12. Plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 6)
)


y_pos = np.arange(
    len(plot_data)
)


scatter = ax.scatter(

    plot_data["RichFactor"],

    y_pos,

    s=point_sizes,

    c=plot_data[
        "neg_log10_FDR"
    ],

    cmap="coolwarm",

    alpha=0.88,

    edgecolors="black",

    linewidths=0.4
)


# ============================================================
# 13. Y axis
# ============================================================

ax.set_yticks(
    y_pos
)

ax.set_yticklabels(
    plot_data["Pathway"],
    fontsize=9
)


# ============================================================
# 14. Axis labels
# ============================================================

ax.set_xlabel(
    "Rich factor",
    fontsize=12
)

ax.set_ylabel(
    "KEGG pathway",
    fontsize=12
)

ax.set_title(
    "KEGG Pathway Enrichment",
    fontsize=14,
    pad=12
)


# ============================================================
# 15. SCI style
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

ax.tick_params(
    axis="x",
    direction="out",
    width=1,
    length=4,
    labelsize=10
)

ax.tick_params(
    axis="y",
    length=0
)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.25
)

ax.set_axisbelow(True)


# ============================================================
# 16. Colorbar
# ============================================================

cbar = fig.colorbar(
    scatter,
    ax=ax,
    fraction=0.035,
    pad=0.04
)

cbar.set_label(
    r"$-\log_{10}(\mathrm{FDR})$",
    fontsize=10
)

cbar.ax.tick_params(
    labelsize=9
)


# ============================================================
# 17. Count legend
# ============================================================

size_values = np.linspace(
    count_min,
    count_max,
    3
)


size_handles = []
size_labels = []


for value in size_values:

    if count_max == count_min:

        size = 180

    else:

        size = (
            80
            +
            (
                value
                - count_min
            )
            /
            (
                count_max
                - count_min
            )
            * 420
        )


    handle = ax.scatter(
        [],
        [],
        s=size,
        facecolor="#BDBDBD",
        edgecolor="black",
        linewidth=0.4
    )


    size_handles.append(
        handle
    )

    size_labels.append(
        str(
            int(
                round(value)
            )
        )
    )


ax.legend(
    size_handles,
    size_labels,
    title="Metabolite count",
    loc="lower right",
    frameon=True,
    edgecolor="#CCCCCC",
    fontsize=9,
    title_fontsize=10
)


# ============================================================
# 18. Save
# ============================================================

plt.tight_layout()


fig.savefig(
    "09_KEGG_Enrichment.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)


fig.savefig(
    "09_KEGG_Enrichment.pdf",
    format="pdf",
    bbox_inches="tight",
    facecolor="white"
)


plt.close(fig)


print(
    "09_KEGG_Enrichment completed."
)
