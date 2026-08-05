# -*- coding: utf-8 -*-
"""
Alpha diversity analysis converted from the supplied R script.

Input files
-----------
tax.txt:
    Rows = taxa/species
    Columns = samples
    First column = taxon name

group.txt:
    Required columns = sample, group

Outputs
-------
alpha_indices.csv
alpha_indices_long.csv
spe_alpha_1.png
spe_alpha.png
"""

from itertools import combinations
from pathlib import Path
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


# ============================================================
# User settings
# ============================================================
WORK_DIR = Path(".")  # Example: Path(r"C:/Users/14266/Desktop/meta")
TAX_FILE = WORK_DIR / "tax.txt"
GROUP_FILE = WORK_DIR / "group.txt"
RANDOM_SEED = 20260805

# Use Times New Roman for all figure text.
# On Windows this font is normally already installed.
mpl.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 13,
    "axes.unicode_minus": False,
})


# ============================================================
# Data loading and validation
# ============================================================
def load_data(tax_file: Path, group_file: Path):
    """Load the abundance table and sample-group table."""
    tax = pd.read_csv(tax_file, sep="\t", index_col=0)
    group = pd.read_csv(group_file, sep="\t", dtype=str)

    required_columns = {"sample", "group"}
    missing = required_columns.difference(group.columns)
    if missing:
        raise ValueError(
            f"group.txt is missing required column(s): {sorted(missing)}"
        )

    if group["sample"].duplicated().any():
        duplicates = group.loc[group["sample"].duplicated(), "sample"].tolist()
        raise ValueError(f"Duplicate sample names in group.txt: {duplicates}")

    # Convert all abundance values to numeric.
    tax = tax.apply(pd.to_numeric, errors="raise")

    if (tax < 0).any().any():
        raise ValueError("The abundance table contains negative values.")

    zero_sum_samples = tax.columns[tax.sum(axis=0) == 0].tolist()
    if zero_sum_samples:
        raise ValueError(
            f"These samples have a total abundance of zero: {zero_sum_samples}"
        )

    tax_samples = set(tax.columns)
    group_samples = set(group["sample"])

    missing_group = sorted(tax_samples - group_samples)
    missing_tax = sorted(group_samples - tax_samples)

    if missing_group:
        raise ValueError(
            "Samples present in tax.txt but absent from group.txt: "
            + ", ".join(missing_group)
        )
    if missing_tax:
        raise ValueError(
            "Samples present in group.txt but absent from tax.txt: "
            + ", ".join(missing_tax)
        )

    # Preserve the sample order from tax.txt.
    group = group.set_index("sample").reindex(tax.columns)
    group.index.name = "sample"
    group = group.reset_index()

    return tax, group


# ============================================================
# Alpha-diversity calculations
# ============================================================
def calculate_alpha_indices(tax: pd.DataFrame, group: pd.DataFrame):
    """
    Reproduce vegan::diversity calculations.

    Shannon:
        -sum(p * ln(p))

    Simpson in vegan:
        1 - sum(p^2)

    Inverse Simpson:
        1 / sum(p^2)
    """
    # R equivalent:
    # data <- data / apply(data, 2, sum)
    relative_abundance = tax.div(tax.sum(axis=0), axis=1)

    # R equivalent:
    # a <- as.data.frame(t(data))
    sample_by_taxon = relative_abundance.T

    p = sample_by_taxon.to_numpy(dtype=float)

    # Avoid log(0): terms with p == 0 contribute zero.
    log_p = np.zeros_like(p)
    positive = p > 0
    log_p[positive] = np.log(p[positive])

    shannon = -(p * log_p).sum(axis=1)
    concentration = np.square(p).sum(axis=1)
    simpson = 1.0 - concentration
    invsimpson = 1.0 / concentration

    alpha = pd.DataFrame({
        "sample": sample_by_taxon.index,
        "shannon": shannon,
        "simpson": simpson,
        "invsimpson": invsimpson,
    })

    alpha = alpha.merge(group[["sample", "group"]], on="sample", how="left")
    return alpha


# ============================================================
# Statistical comparisons
# ============================================================
def pairwise_mannwhitney(
    data: pd.DataFrame,
    metric: str,
    comparisons,
):
    """
    Pairwise two-sided Mann-Whitney U tests.

    This corresponds closely to ggpubr::stat_compare_means() with
    pairwise Wilcoxon rank-sum comparisons.
    """
    results = []

    for group_1, group_2 in comparisons:
        values_1 = data.loc[data["group"] == group_1, metric].dropna()
        values_2 = data.loc[data["group"] == group_2, metric].dropna()

        if len(values_1) == 0 or len(values_2) == 0:
            p_value = np.nan
        else:
            _, p_value = mannwhitneyu(
                values_1,
                values_2,
                alternative="two-sided",
                method="auto",
            )

        results.append({
            "group_1": group_1,
            "group_2": group_2,
            "p_value": p_value,
        })

    return pd.DataFrame(results)


def format_p_value(p_value: float):
    if np.isnan(p_value):
        return "p = NA"
    if p_value < 0.0001:
        return "p < 0.0001"
    return f"p = {p_value:.4f}"


# ============================================================
# Plotting functions
# ============================================================
def add_p_value_brackets(
    ax,
    data: pd.DataFrame,
    metric: str,
    group_order,
    comparisons,
):
    """Add pairwise p-value brackets above the boxplots."""
    test_results = pairwise_mannwhitney(data, metric, comparisons)

    all_values = data[metric].dropna().to_numpy()
    data_min = float(np.min(all_values))
    data_max = float(np.max(all_values))
    data_range = data_max - data_min

    if data_range == 0:
        data_range = max(abs(data_max), 1.0)

    step = data_range * 0.14
    bracket_height = data_range * 0.025
    first_y = data_max + data_range * 0.10

    position_map = {
        group_name: position
        for position, group_name in enumerate(group_order, start=1)
    }

    for index, row in test_results.iterrows():
        x1 = position_map[row["group_1"]]
        x2 = position_map[row["group_2"]]
        y = first_y + index * step

        ax.plot(
            [x1, x1, x2, x2],
            [y, y + bracket_height, y + bracket_height, y],
            linewidth=1.0,
        )
        ax.text(
            (x1 + x2) / 2,
            y + bracket_height,
            format_p_value(row["p_value"]),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    upper_limit = first_y + max(len(test_results), 1) * step + data_range * 0.08
    lower_limit = data_min - data_range * 0.08
    ax.set_ylim(lower_limit, upper_limit)


def draw_alpha_panel(
    ax,
    data: pd.DataFrame,
    metric: str,
    ylabel: str,
    group_order,
    comparisons,
    rng,
    panel_title=None,
):
    """Draw one boxplot panel with jittered sample points."""
    values_by_group = [
        data.loc[data["group"] == group_name, metric].dropna().to_numpy()
        for group_name in group_order
    ]

    positions = np.arange(1, len(group_order) + 1)
    default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    group_colors = [
        default_colors[i % len(default_colors)]
        for i in range(len(group_order))
    ]

    boxplot = ax.boxplot(
        values_by_group,
        positions=positions,
        widths=0.60,
        patch_artist=True,
        showfliers=False,
        medianprops={"linewidth": 1.2},
    )

    for patch, color in zip(boxplot["boxes"], group_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)

    for position, values, color in zip(
        positions,
        values_by_group,
        group_colors,
    ):
        jitter = rng.normal(loc=0.0, scale=0.055, size=len(values))
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=18,
            alpha=0.75,
            color=color,
            zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(group_order, rotation=30, ha="right")
    ax.set_xlabel("Group")
    ax.set_ylabel(ylabel)
    ax.grid(False)

    if panel_title:
        ax.set_title(panel_title)

    add_p_value_brackets(
        ax=ax,
        data=data,
        metric=metric,
        group_order=group_order,
        comparisons=comparisons,
    )

    return [
        Patch(facecolor=color, alpha=0.65, label=group_name)
        for color, group_name in zip(group_colors, group_order)
    ]


def save_three_panel_plot(alpha: pd.DataFrame, output_file: Path):
    """Equivalent to plot_grid(p1, p2, p3, ncol = 3)."""
    group_order = list(dict.fromkeys(alpha["group"]))
    comparisons = list(combinations(group_order, 2))
    rng = np.random.default_rng(RANDOM_SEED)

    fig, axes = plt.subplots(1, 3, figsize=(8.9, 5.9))

    metrics = [
        ("shannon", "Shannon index"),
        ("simpson", "Simpson index"),
        ("invsimpson", "Invsimpson index"),
    ]

    legend_handles = None
    for ax, (metric, ylabel) in zip(axes, metrics):
        legend_handles = draw_alpha_panel(
            ax=ax,
            data=alpha,
            metric=metric,
            ylabel=ylabel,
            group_order=group_order,
            comparisons=comparisons,
            rng=rng,
        )

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=len(group_order),
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


def save_faceted_plot(alpha_long: pd.DataFrame, output_file: Path):
    """
    Equivalent to facet_wrap(~ variable, scales = "free_y").

    Each panel has an independent y-axis scale.
    """
    metric_labels = {
        "shannon": "Shannon",
        "simpson": "Simpson",
        "invsimpson": "Invsimpson",
    }

    group_order = list(dict.fromkeys(alpha_long["group"]))
    comparisons = list(combinations(group_order, 2))
    rng = np.random.default_rng(RANDOM_SEED)

    fig, axes = plt.subplots(1, 3, figsize=(8.9, 5.9))

    legend_handles = None
    for ax, metric in zip(
        axes,
        ["shannon", "simpson", "invsimpson"],
    ):
        panel_data = alpha_long.loc[
            alpha_long["variable"] == metric,
            ["sample", "group", "value"],
        ].rename(columns={"value": metric})

        legend_handles = draw_alpha_panel(
            ax=ax,
            data=panel_data,
            metric=metric,
            ylabel="Alpha index",
            group_order=group_order,
            comparisons=comparisons,
            rng=rng,
            panel_title=metric_labels[metric],
        )

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=len(group_order),
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main workflow
# ============================================================
def main():
    tax, group = load_data(TAX_FILE, GROUP_FILE)
    alpha = calculate_alpha_indices(tax, group)

    # Wide-format result.
    alpha.to_csv(
        WORK_DIR / "alpha_indices.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # R equivalent:
    # reshape2::melt(spe_alpha, id = c("sample", "group"))
    alpha_long = alpha.melt(
        id_vars=["sample", "group"],
        value_vars=["shannon", "simpson", "invsimpson"],
        var_name="variable",
        value_name="value",
    )

    alpha_long.to_csv(
        WORK_DIR / "alpha_indices_long.csv",
        index=False,
        encoding="utf-8-sig",
    )

    save_three_panel_plot(
        alpha,
        WORK_DIR / "spe_alpha_1.png",
    )
    save_faceted_plot(
        alpha_long,
        WORK_DIR / "spe_alpha.png",
    )

    print("Analysis completed.")
    print(f"Samples: {len(alpha)}")
    print(f"Groups: {', '.join(alpha['group'].drop_duplicates())}")
    print("Generated:")
    print("  alpha_indices.csv")
    print("  alpha_indices_long.csv")
    print("  spe_alpha_1.png")
    print("  spe_alpha.png")


if __name__ == "__main__":
    main()
