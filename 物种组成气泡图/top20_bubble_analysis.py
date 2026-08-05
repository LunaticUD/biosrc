# -*- coding: utf-8 -*-
"""
Top-20 species relative-abundance bubble plot.

Converted from the supplied R/ggplot2 workflow.

Input files
-----------
tax_bubble.txt:
    Rows = species
    Columns = samples
    First column = species name

group_bubble.txt:
    Required columns = sample, group

Outputs
-------
top20_relative_abundance.csv
top20_bubble_long.csv
top20_bubble.png
"""

from pathlib import Path
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# User settings
# ============================================================
WORK_DIR = Path(".")
# For your original folder:
# WORK_DIR = Path(r"C:/Users/14266/Desktop/meta")

TAX_FILE = WORK_DIR / "tax_bubble.txt"
GROUP_FILE = WORK_DIR / "group_bubble.txt"

TOP_N = 20
GROUP_ORDER = ["KC", "WB", "MC"]

GROUP_COLORS = {
    "KC": "#edae11",
    "WB": "#0f6657",
    "MC": "#c74732",
}

# Equivalent intent to:
# scale_size(range=c(0, 15), breaks=c(0,10,20,40), limits=c(0,400))
SIZE_BREAKS = [0, 10, 20, 40]
SIZE_LIMIT = 400
MAX_MARKER_DIAMETER_PT = 42

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# Times New Roman is normally installed on Windows.
# DejaVu Serif is only a fallback for environments where it is unavailable.
mpl.rcParams.update({
    "font.family": ["Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 9,
    "axes.unicode_minus": False,
})


# ============================================================
# Data loading and validation
# ============================================================
def load_data(tax_file: Path, group_file: Path):
    tax = pd.read_csv(tax_file, sep="\t", index_col=0)
    group = pd.read_csv(group_file, sep="\t", dtype=str)

    required = {"sample", "group"}
    missing = required.difference(group.columns)
    if missing:
        raise ValueError(
            f"group file is missing required column(s): {sorted(missing)}"
        )

    if group["sample"].duplicated().any():
        duplicates = group.loc[
            group["sample"].duplicated(),
            "sample",
        ].tolist()
        raise ValueError(f"Duplicate sample names: {duplicates}")

    tax = tax.apply(pd.to_numeric, errors="raise")

    if (tax < 0).any().any():
        raise ValueError("The abundance table contains negative values.")

    zero_total = tax.columns[tax.sum(axis=0) == 0].tolist()
    if zero_total:
        raise ValueError(
            "Samples with zero total abundance: " + ", ".join(zero_total)
        )

    tax_samples = set(tax.columns)
    group_samples = set(group["sample"])

    missing_group = sorted(tax_samples - group_samples)
    missing_tax = sorted(group_samples - tax_samples)

    if missing_group:
        raise ValueError(
            "Samples in the abundance table but not in group file: "
            + ", ".join(missing_group)
        )
    if missing_tax:
        raise ValueError(
            "Samples in group file but not in abundance table: "
            + ", ".join(missing_tax)
        )

    # Align by sample name rather than relying on row/column position.
    group = group.set_index("sample").reindex(tax.columns)
    group.index.name = "sample"
    group = group.reset_index()

    return tax, group


# ============================================================
# Data preparation
# ============================================================
def prepare_top_species(tax: pd.DataFrame, group: pd.DataFrame):
    """
    Match the R workflow:
      1. Rank species by total abundance.
      2. Keep the top 20.
      3. Divide by each sample's total abundance and multiply by 100.
      4. Transpose and combine with sample grouping.
      5. Convert wide data to long data.
    """
    row_sums = tax.sum(axis=1)
    ranked_tax = tax.loc[row_sums.sort_values(ascending=False).index]

    top_n = min(TOP_N, len(ranked_tax))

    # Important: denominator is the total abundance of all species,
    # matching data1[1:20, ] / apply(data1, 2, sum) * 100.
    sample_totals = ranked_tax.sum(axis=0)
    top_relative = (
        ranked_tax.iloc[:top_n, :]
        .div(sample_totals, axis=1)
        .mul(100)
    )

    # Samples x top species.
    sample_wide = top_relative.T
    sample_wide.index.name = "sample"
    sample_wide = sample_wide.reset_index()

    data_wide = group[["sample", "group"]].merge(
        sample_wide,
        on="sample",
        how="left",
        validate="one_to_one",
    )

    # Explicit equivalent of reshape2::melt(..., id.vars=c("sample","group")).
    data_long = data_wide.melt(
        id_vars=["sample", "group"],
        var_name="species",
        value_name="relative_abundance",
    )

    data_long["group"] = pd.Categorical(
        data_long["group"],
        categories=GROUP_ORDER,
        ordered=True,
    )

    return top_relative, data_long


# ============================================================
# Plotting
# ============================================================
def bubble_area(values):
    """
    Map abundance to marker area.

    ggplot2's size scale uses perceptual area scaling. Matplotlib's scatter
    parameter `s` is already marker area in points squared.
    """
    values = np.asarray(values, dtype=float)
    normalized = np.clip(values / SIZE_LIMIT, 0, 1)
    return normalized * (MAX_MARKER_DIAMETER_PT ** 2)


def plot_bubble(data_long: pd.DataFrame, species_order, output_file: Path):
    present_groups = [
        group for group in GROUP_ORDER
        if group in set(data_long["group"].dropna().astype(str))
    ]

    if not present_groups:
        present_groups = list(
            dict.fromkeys(data_long["group"].dropna().astype(str))
        )

    sample_counts = [
        data_long.loc[
            data_long["group"].astype(str) == group,
            "sample",
        ].nunique()
        for group in present_groups
    ]

    width_ratios = [
        max(count, 1) for count in sample_counts
    ]

    fig, axes = plt.subplots(
        1,
        len(present_groups),
        figsize=(10, 6.2),
        sharey=True,
        gridspec_kw={
            "width_ratios": width_ratios,
            "wspace": 0.05,
        },
    )

    if len(present_groups) == 1:
        axes = [axes]

    # Highest-ranked species appears at the top.
    y_positions = {
        species: len(species_order) - 1 - i
        for i, species in enumerate(species_order)
    }

    for panel_index, (ax, group_name) in enumerate(
        zip(axes, present_groups)
    ):
        panel = data_long.loc[
            data_long["group"].astype(str) == group_name
        ].copy()

        sample_order = list(dict.fromkeys(panel["sample"]))
        x_positions = {
            sample: i for i, sample in enumerate(sample_order)
        }

        x = panel["sample"].map(x_positions).to_numpy()
        y = panel["species"].map(y_positions).to_numpy()
        sizes = bubble_area(panel["relative_abundance"].to_numpy())

        ax.scatter(
            x,
            y,
            s=sizes,
            color=GROUP_COLORS.get(group_name, "#666666"),
            alpha=0.90,
            edgecolors="none",
        )

        ax.set_title(group_name, pad=7)
        ax.set_xticks(range(len(sample_order)))
        ax.set_xticklabels(
            sample_order,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax.set_xlim(-0.6, len(sample_order) - 0.4)
        ax.set_ylim(-0.7, len(species_order) - 0.3)
        ax.grid(False)

        # White panel background with a black border, similar to theme_bw().
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(0.8)

        if panel_index == 0:
            ax.set_yticks(range(len(species_order)))
            ax.set_yticklabels(list(reversed(species_order)))
            ax.set_ylabel("Species")
        else:
            ax.tick_params(axis="y", left=False, labelleft=False)

        ax.set_xlabel("Samples")

    # Group color legend.
    color_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=8,
            markerfacecolor=GROUP_COLORS.get(group, "#666666"),
            markeredgecolor=GROUP_COLORS.get(group, "#666666"),
            label=group,
        )
        for group in present_groups
    ]

    color_legend = fig.legend(
        handles=color_handles,
        title="Group",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.855, 0.91),
    )

    # Bubble-size legend.
    size_handles = [
        plt.scatter(
            [],
            [],
            s=bubble_area([value])[0],
            color="#808080",
            alpha=0.85,
            edgecolors="none",
            label=str(value),
        )
        for value in SIZE_BREAKS
    ]

    size_legend = fig.legend(
        handles=size_handles,
        title="Relative Abundance (%)",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.855, 0.53),
        labelspacing=1.25,
    )

    fig.add_artist(color_legend)
    fig.add_artist(size_legend)

    fig.subplots_adjust(
        left=0.28,
        right=0.84,
        bottom=0.23,
        top=0.91,
    )
    fig.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.10,
    )
    plt.close(fig)


# ============================================================
# Main workflow
# ============================================================
def main():
    tax, group = load_data(TAX_FILE, GROUP_FILE)
    top_relative, data_long = prepare_top_species(tax, group)

    top_relative.to_csv(
        WORK_DIR / "top20_relative_abundance.csv",
        encoding="utf-8-sig",
    )
    data_long.to_csv(
        WORK_DIR / "top20_bubble_long.csv",
        index=False,
        encoding="utf-8-sig",
    )

    species_order = list(top_relative.index)

    plot_bubble(
        data_long,
        species_order,
        WORK_DIR / "top20_bubble.png",
    )

    print("Analysis completed.")
    print(f"Samples: {tax.shape[1]}")
    print(f"Species in input: {tax.shape[0]}")
    print(f"Species plotted: {len(species_order)}")
    print(
        "Groups: "
        + ", ".join(
            group["group"].drop_duplicates().astype(str)
        )
    )
    print("Generated:")
    print("  top20_relative_abundance.csv")
    print("  top20_bubble_long.csv")
    print("  top20_bubble.png")


if __name__ == "__main__":
    main()
