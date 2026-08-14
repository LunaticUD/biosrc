# ============================================================
# 08_KEGG_Classification.py
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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

print(
    f"Differential metabolites loaded: {len(sig)}"
)


# ============================================================
# 3. KEGG annotation
#
# 如果没有真实注释文件，则自动生成示例数据
# ============================================================

annotation_file = "Metabolite_KEGG_Annotation.csv"


if not os.path.exists(annotation_file):

    print(
        "Metabolite_KEGG_Annotation.csv not found."
    )

    print(
        "Generating example KEGG annotations..."
    )

    np.random.seed(42)

    pathway_class_map = {

        "ABC transporters":
            "Environmental Information Processing",

        "Biosynthesis of amino acids":
            "Metabolism",

        "Purine metabolism":
            "Metabolism",

        "Pyrimidine metabolism":
            "Metabolism",

        "Glycerophospholipid metabolism":
            "Metabolism",

        "Fatty acid metabolism":
            "Metabolism",

        "Citrate cycle (TCA cycle)":
            "Metabolism",

        "Glycolysis / Gluconeogenesis":
            "Metabolism",

        "Phenylalanine metabolism":
            "Metabolism",

        "Glutathione metabolism":
            "Metabolism",

        "Protein digestion and absorption":
            "Organismal Systems",

        "Apoptosis":
            "Cellular Processes"
    }

    pathways = list(
        pathway_class_map.keys()
    )

    annotation_df = pd.DataFrame({
        "Metabolite": sig.index,
        "Pathway": np.random.choice(
            pathways,
            size=len(sig),
            replace=True
        )
    })

    annotation_df["Class"] = (
        annotation_df["Pathway"]
        .map(pathway_class_map)
    )

    annotation_df.to_csv(
        annotation_file,
        index=False,
        encoding="utf-8-sig"
    )

    print(
        f"Example annotation saved to: "
        f"{annotation_file}"
    )

else:

    annotation_df = pd.read_csv(
        annotation_file
    )


# ============================================================
# 4. Keep only differential metabolites
# ============================================================

annotation_df = annotation_df[
    annotation_df["Metabolite"].isin(
        sig.index
    )
].copy()


if annotation_df.empty:

    raise ValueError(
        "No differential metabolites matched KEGG annotation."
    )


# ============================================================
# 5. Count pathway hits
# ============================================================

kegg_class_df = (
    annotation_df
    .groupby(
        [
            "Pathway",
            "Class"
        ]
    )
    .size()
    .reset_index(
        name="Count"
    )
)


# ============================================================
# 6. Percentage
# ============================================================

total_hits = (
    kegg_class_df["Count"].sum()
)

kegg_class_df["Percentage"] = (
    kegg_class_df["Count"]
    / total_hits
    * 100
).round(1)


# 保存用于检查
kegg_class_df.to_csv(
    "KEGG_Classification_Result.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 7. Sort
# ============================================================

plot_data = kegg_class_df.sort_values(
    "Count",
    ascending=True
).copy()


# ============================================================
# 8. Class colors
# ============================================================

class_palette = {

    "Metabolism":
        "#66C2A5",

    "Environmental Information Processing":
        "#FC8D62",

    "Genetic Information Processing":
        "#8DA0CB",

    "Cellular Processes":
        "#E78AC3",

    "Organismal Systems":
        "#A6D854",

    "Human Diseases":
        "#FFD92F"
}


bar_colors = [

    class_palette.get(
        cls,
        "#BDBDBD"
    )

    for cls in plot_data["Class"]
]


# ============================================================
# 9. Plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 6)
)


bars = ax.barh(
    plot_data["Pathway"],
    plot_data["Count"],
    color=bar_colors,
    height=0.72,
    alpha=0.88,
    edgecolor="none"
)


# ============================================================
# 10. Labels
# ============================================================

max_count = (
    plot_data["Count"].max()
)

offset = max_count * 0.025


for bar, (_, row) in zip(
    bars,
    plot_data.iterrows()
):

    ax.text(
        row["Count"] + offset,
        bar.get_y()
        + bar.get_height() / 2,

        f'{int(row["Count"])} '
        f'({row["Percentage"]:.1f}%)',

        va="center",
        ha="left",
        fontsize=9
    )


# ============================================================
# 11. Axis
# ============================================================

ax.set_xlabel(
    "Number of differential metabolites",
    fontsize=12
)

ax.set_ylabel("")

ax.set_title(
    "KEGG Pathway Classification",
    fontsize=14,
    pad=12
)


ax.tick_params(
    axis="x",
    direction="out",
    width=1,
    length=4,
    labelsize=10
)

ax.tick_params(
    axis="y",
    length=0,
    labelsize=9
)


# ============================================================
# 12. SCI style
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.25
)

ax.set_axisbelow(True)


# ============================================================
# 13. Legend
# ============================================================

used_classes = list(
    dict.fromkeys(
        plot_data["Class"]
    )
)


legend_handles = [

    Patch(
        facecolor=class_palette.get(
            cls,
            "#BDBDBD"
        ),
        edgecolor="none",
        label=cls
    )

    for cls in used_classes
]


ax.legend(
    handles=legend_handles,
    title="KEGG Class",
    loc="lower right",
    frameon=True,
    edgecolor="#CCCCCC",
    fontsize=9,
    title_fontsize=10
)


# ============================================================
# 14. X limits
# ============================================================

ax.set_xlim(
    0,
    max_count * 1.25
)


# ============================================================
# 15. Save
# ============================================================

plt.tight_layout()


fig.savefig(
    "08_KEGG_Classification.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)


fig.savefig(
    "08_KEGG_Classification.pdf",
    format="pdf",
    bbox_inches="tight",
    facecolor="white"
)


plt.close(fig)


print(
    "08_KEGG_Classification completed."
)
