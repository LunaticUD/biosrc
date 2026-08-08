# -*- coding: utf-8 -*-
"""
Bray-Curtis PCoA + PERMANOVA + confidence ellipses + marginal densities.

Converted from the supplied R workflow.

Required input files
--------------------
tax.txt
    Rows = taxa/species
    Columns = samples
    First column = taxon name

group.txt
    Required columns = sample, group

Generated outputs
-----------------
bray-crutis.txt
pcoa_coordinates.csv
permanova_results.csv
bra_pcoa1.png
bray_pcoa.png
"""

from pathlib import Path
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2, gaussian_kde

# ============================================================
# User settings
# ============================================================
WORK_DIR = Path(".")
# Example:
# WORK_DIR = Path(r"C:/Users/14266/Desktop/meta")

TAX_FILE = WORK_DIR / "tax.txt"
GROUP_FILE = WORK_DIR / "group.txt"

PERMUTATIONS = 999
ELLIPSE_LEVEL = 0.90
RANDOM_SEED = 20260805

POINT_COLORS = ["#2f2f2f", "#972b1d", "#e6a84b"]
FILL_COLORS = ["#f2ccac", "#a1d5b9", "#e1abbc"]
MARKERS = ["o", "s", "^", "D", "P", "X"]

# Configure all chart text as Times New Roman.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
mpl.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.unicode_minus": False,
})

# ============================================================
# Data loading and validation
# ============================================================
def load_data(tax_file: Path, group_file: Path):
    """Load and align the species abundance and grouping tables."""
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

    tax = tax.apply(pd.to_numeric, errors="raise")

    if (tax < 0).any().any():
        raise ValueError("The abundance table contains negative values.")

    zero_sum_samples = tax.columns[tax.sum(axis=0) == 0].tolist()
    if zero_sum_samples:
        raise ValueError(
            "Samples with zero total abundance: " + ", ".join(zero_sum_samples)
        )

    tax_samples = set(tax.columns)
    group_samples = set(group["sample"])

    only_in_tax = sorted(tax_samples - group_samples)
    only_in_group = sorted(group_samples - tax_samples)

    if only_in_tax:
        raise ValueError(
            "Samples present in tax.txt but absent from group.txt: "
            + ", ".join(only_in_tax)
        )
    if only_in_group:
        raise ValueError(
            "Samples present in group.txt but absent from tax.txt: "
            + ", ".join(only_in_group)
        )

    # Preserve sample order from tax.txt.
    group = group.set_index("sample").reindex(tax.columns)
    group.index.name = "sample"
    group = group.reset_index()

    return tax, group


# ============================================================
# Bray-Curtis distance and PCoA
# ============================================================
def relative_abundance_by_sample(tax: pd.DataFrame):
    """
    Match the R preprocessing:
        data <- data / apply(data, 2, sum)
        data <- t(data)
    """
    relative = tax.div(tax.sum(axis=0), axis=1)
    return relative.T


def bray_curtis_matrix(sample_by_taxon: pd.DataFrame):
    """Calculate a square Bray-Curtis distance matrix."""
    condensed = pdist(sample_by_taxon.to_numpy(dtype=float), metric="braycurtis")
    matrix = squareform(condensed)

    return pd.DataFrame(
        matrix,
        index=sample_by_taxon.index,
        columns=sample_by_taxon.index,
    )


def classical_pcoa(distance_matrix: pd.DataFrame, dimensions: int = 3):
    """
    Classical principal coordinates analysis.

    For Bray-Curtis distances, small negative eigenvalues may occur because
    the distance is not always perfectly Euclidean. Coordinates are built
    from positive eigenvalues only.
    """
    distances = distance_matrix.to_numpy(dtype=float)
    n = distances.shape[0]

    centering = np.eye(n) - np.ones((n, n)) / n
    b_matrix = -0.5 * centering @ np.square(distances) @ centering

    eigenvalues, eigenvectors = np.linalg.eigh(b_matrix)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive_mask = eigenvalues > np.finfo(float).eps
    positive_eigenvalues = eigenvalues[positive_mask]
    positive_eigenvectors = eigenvectors[:, positive_mask]

    if len(positive_eigenvalues) < dimensions:
        raise ValueError(
            f"Only {len(positive_eigenvalues)} positive PCoA axes are available."
        )

    coordinates = (
        positive_eigenvectors[:, :dimensions]
        * np.sqrt(positive_eigenvalues[:dimensions])
    )

    # Corrected version of the R expression:
    # eig_percent <- eig / sum(positive eig) * 100
    explained_percent = (
        positive_eigenvalues / positive_eigenvalues.sum() * 100.0
    )

    coordinate_columns = [
        f"PCoA{i}" for i in range(1, dimensions + 1)
    ]
    coordinates_df = pd.DataFrame(
        coordinates,
        index=distance_matrix.index,
        columns=coordinate_columns,
    )
    coordinates_df.index.name = "sample"

    return coordinates_df, eigenvalues, explained_percent


# ============================================================
# One-factor PERMANOVA
# ============================================================
def permanova_one_factor(
    distance_matrix: pd.DataFrame,
    groups,
    permutations: int = 999,
    seed: int = 20260805,
):
    """
    One-factor PERMANOVA based on a distance matrix.

    This follows the projection-matrix formulation underlying adonis2:
      R2 = model SS / total SS
      pseudo-F = (model SS / df_model) / (residual SS / df_residual)

    Group labels are permuted to obtain the p-value.
    """
    distances = distance_matrix.to_numpy(dtype=float)
    labels = np.asarray(groups, dtype=object)
    n = len(labels)

    if distances.shape != (n, n):
        raise ValueError("Distance matrix size and group vector length differ.")

    unique_groups = list(dict.fromkeys(labels))
    number_of_groups = len(unique_groups)

    if number_of_groups < 2:
        raise ValueError("PERMANOVA requires at least two groups.")

    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ np.square(distances) @ centering
    total_ss = float(np.trace(gower))

    def model_statistics(current_labels):
        dummy = pd.get_dummies(
            pd.Categorical(current_labels, categories=unique_groups),
            drop_first=False,
            dtype=float,
        ).to_numpy()

        design = np.column_stack([np.ones(n), dummy[:, 1:]])
        hat = design @ np.linalg.pinv(design.T @ design) @ design.T

        model_ss = float(np.trace(hat @ gower))
        residual_ss = total_ss - model_ss

        df_model = number_of_groups - 1
        df_residual = n - number_of_groups

        pseudo_f = (
            (model_ss / df_model)
            / (residual_ss / df_residual)
        )
        r_squared = model_ss / total_ss

        return pseudo_f, r_squared, model_ss, residual_ss

    observed_f, r_squared, model_ss, residual_ss = model_statistics(labels)

    rng = np.random.default_rng(seed)
    permuted_f = np.empty(permutations, dtype=float)

    for i in range(permutations):
        permuted_labels = rng.permutation(labels)
        permuted_f[i] = model_statistics(permuted_labels)[0]

    p_value = (
        1 + np.count_nonzero(permuted_f >= observed_f)
    ) / (permutations + 1)

    return {
        "pseudo_F": observed_f,
        "R2": r_squared,
        "p_value": p_value,
        "permutations": permutations,
        "model_SS": model_ss,
        "residual_SS": residual_ss,
        "df_model": number_of_groups - 1,
        "df_residual": n - number_of_groups,
    }


# ============================================================
# Plot helpers
# ============================================================
def build_group_styles(group_order):
    """Assign colors and marker shapes in group order."""
    styles = {}
    for i, group_name in enumerate(group_order):
        styles[group_name] = {
            "point_color": POINT_COLORS[i % len(POINT_COLORS)],
            "fill_color": FILL_COLORS[i % len(FILL_COLORS)],
            "marker": MARKERS[i % len(MARKERS)],
        }
    return styles


def style_main_axis(ax, x_label, y_label):
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.axhline(0, color="#BEBEBE", linestyle="--", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#BEBEBE", linestyle="--", linewidth=0.8, zorder=0)
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color("black")

    ax.tick_params(colors="black")


def scatter_groups(ax, pcoa_result, group_order, styles):
    for group_name in group_order:
        subset = pcoa_result.loc[pcoa_result["group"] == group_name]
        style = styles[group_name]

        ax.scatter(
            subset["PCoA1"],
            subset["PCoA2"],
            s=65,
            c=style["point_color"],
            marker=style["marker"],
            label=group_name,
            edgecolors="none",
            zorder=3,
        )


def add_confidence_ellipses(
    ax,
    pcoa_result,
    group_order,
    styles,
    level=0.90,
):
    """
    Add covariance-based confidence ellipses, comparable to stat_ellipse().
    """
    scale = np.sqrt(chi2.ppf(level, df=2))

    for group_name in group_order:
        subset = pcoa_result.loc[
            pcoa_result["group"] == group_name,
            ["PCoA1", "PCoA2"],
        ].dropna()

        if len(subset) < 3:
            continue

        values = subset.to_numpy(dtype=float)
        covariance = np.cov(values, rowvar=False)

        if not np.all(np.isfinite(covariance)):
            continue

        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        if np.any(eigenvalues <= 0):
            continue

        angle = np.degrees(
            np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
        )
        width, height = 2 * scale * np.sqrt(eigenvalues)

        style = styles[group_name]
        ellipse = Ellipse(
            xy=values.mean(axis=0),
            width=width,
            height=height,
            angle=angle,
            facecolor=style["fill_color"],
            edgecolor=style["point_color"],
            linestyle="--",
            linewidth=0.8,
            alpha=0.30,
            zorder=1,
        )
        ax.add_patch(ellipse)


def density_values(values, grid):
    """Return KDE values, handling nearly constant data safely."""
    values = np.asarray(values, dtype=float)

    if len(values) < 2 or np.isclose(np.std(values, ddof=1), 0):
        return np.zeros_like(grid)

    kde = gaussian_kde(values)
    return kde(grid)


def make_legend_handles(group_order, styles):
    handles = []
    for group_name in group_order:
        style = styles[group_name]
        handles.append(
            Line2D(
                [0],
                [0],
                marker=style["marker"],
                linestyle="None",
                markerfacecolor=style["point_color"],
                markeredgecolor=style["point_color"],
                markersize=7,
                label=group_name,
            )
        )
    return handles


def save_simple_pcoa(
    pcoa_result,
    explained_percent,
    permanova_result,
    output_file,
):
    """Create the simple PCoA plot."""
    group_order = list(dict.fromkeys(pcoa_result["group"]))
    styles = build_group_styles(group_order)

    fig, ax = plt.subplots(figsize=(5, 5))

    scatter_groups(ax, pcoa_result, group_order, styles)

    x_label = f"PCoA 1 ({explained_percent[0]:.1f}%)"
    y_label = f"PCoA 2 ({explained_percent[1]:.1f}%)"
    style_main_axis(ax, x_label, y_label)

    legend = ax.legend(
        handles=make_legend_handles(group_order, styles),
        title=None,
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.08),
    )
    ax.add_artist(legend)

    caption = (
        f"PERMANOVA R²: {permanova_result['R2']:.2f}; "
        f"P-value: {permanova_result['p_value']:.3f}"
    )
    fig.text(0.5, 0.015, caption, ha="center", va="bottom", fontsize=9)

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


def save_enhanced_pcoa(
    pcoa_result,
    explained_percent,
    permanova_result,
    output_file,
):
    """
    Create the enhanced PCoA plot with:
      - 90% group confidence ellipses
      - top and right marginal density plots
    """
    group_order = list(dict.fromkeys(pcoa_result["group"]))
    styles = build_group_styles(group_order)

    fig = plt.figure(figsize=(5, 5))
    grid = GridSpec(
        2,
        2,
        width_ratios=[4.0, 1.15],
        height_ratios=[1.15, 4.0],
        hspace=0.04,
        wspace=0.04,
    )

    ax_top = fig.add_subplot(grid[0, 0])
    ax_main = fig.add_subplot(grid[1, 0])
    ax_right = fig.add_subplot(grid[1, 1])
    ax_empty = fig.add_subplot(grid[0, 1])
    ax_empty.axis("off")

    scatter_groups(ax_main, pcoa_result, group_order, styles)
    add_confidence_ellipses(
        ax_main,
        pcoa_result,
        group_order,
        styles,
        level=ELLIPSE_LEVEL,
    )

    x_label = f"PCoA 1 ({explained_percent[0]:.1f}%)"
    y_label = f"PCoA 2 ({explained_percent[1]:.1f}%)"
    style_main_axis(ax_main, x_label, y_label)

    legend = ax_main.legend(
        handles=make_legend_handles(group_order, styles),
        title=None,
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.08),
    )
    ax_main.add_artist(legend)

    x_min, x_max = ax_main.get_xlim()
    y_min, y_max = ax_main.get_ylim()
    x_grid = np.linspace(x_min, x_max, 300)
    y_grid = np.linspace(y_min, y_max, 300)

    for group_name in group_order:
        subset = pcoa_result.loc[pcoa_result["group"] == group_name]
        style = styles[group_name]

        x_density = density_values(subset["PCoA1"], x_grid)
        y_density = density_values(subset["PCoA2"], y_grid)

        ax_top.fill_between(
            x_grid,
            0,
            x_density,
            color=style["fill_color"],
            alpha=0.55,
            linewidth=0,
        )
        ax_right.fill_betweenx(
            y_grid,
            0,
            y_density,
            color=style["fill_color"],
            alpha=0.55,
            linewidth=0,
        )

    ax_top.set_xlim(ax_main.get_xlim())
    ax_right.set_ylim(ax_main.get_ylim())

    ax_top.axis("off")
    ax_right.axis("off")

    caption = (
        f"PERMANOVA R²: {permanova_result['R2']:.2f}; "
        f"P-value: {permanova_result['p_value']:.3f}"
    )
    fig.text(0.45, 0.008, caption, ha="center", va="bottom", fontsize=9)

    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.13, top=0.97)
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main workflow
# ============================================================
def main():
    tax, group = load_data(TAX_FILE, GROUP_FILE)

    sample_by_taxon = relative_abundance_by_sample(tax)
    bray = bray_curtis_matrix(sample_by_taxon)

    # The filename spelling is preserved from the supplied R code.
    bray.to_csv(
        WORK_DIR / "bray-crutis.txt",
        sep="\t",
        float_format="%.10f",
    )

    coordinates, eigenvalues, explained_percent = classical_pcoa(
        bray,
        dimensions=3,
    )

    pcoa_result = (
        coordinates.reset_index()
        .merge(group[["sample", "group"]], on="sample", how="left")
    )
    pcoa_result.to_csv(
        WORK_DIR / "pcoa_coordinates.csv",
        index=False,
        encoding="utf-8-sig",
    )

    permanova_result = permanova_one_factor(
        distance_matrix=bray,
        groups=group["group"].to_numpy(),
        permutations=PERMUTATIONS,
        seed=RANDOM_SEED,
    )

    pd.DataFrame([permanova_result]).to_csv(
        WORK_DIR / "permanova_results.csv",
        index=False,
        encoding="utf-8-sig",
    )

    save_simple_pcoa(
        pcoa_result,
        explained_percent,
        permanova_result,
        WORK_DIR / "bra_pcoa1.png",
    )

    save_enhanced_pcoa(
        pcoa_result,
        explained_percent,
        permanova_result,
        WORK_DIR / "bray_pcoa.png",
    )

    positive_count = int(np.count_nonzero(eigenvalues > 0))
    negative_count = int(np.count_nonzero(eigenvalues < 0))

    print("Analysis completed.")
    print(f"Samples: {len(pcoa_result)}")
    print(f"Groups: {', '.join(pcoa_result['group'].drop_duplicates())}")
    print(
        "PCoA explained variance: "
        f"Axis 1 = {explained_percent[0]:.1f}%, "
        f"Axis 2 = {explained_percent[1]:.1f}%, "
        f"Axis 3 = {explained_percent[2]:.1f}%"
    )
    print(
        f"Eigenvalues: {positive_count} positive, "
        f"{negative_count} negative"
    )
    print(
        "PERMANOVA: "
        f"R2 = {permanova_result['R2']:.4f}, "
        f"pseudo-F = {permanova_result['pseudo_F']:.4f}, "
        f"p = {permanova_result['p_value']:.4f}"
    )
    print("Generated:")
    print("  bray-crutis.txt")
    print("  pcoa_coordinates.csv")
    print("  permanova_results.csv")
    print("  bra_pcoa1.png")
    print("  bray_pcoa.png")


if __name__ == "__main__":
    main()
