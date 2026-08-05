# -*- coding: utf-8 -*-
"""
Bray-Curtis NMDS analysis converted from the supplied R code.

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
bray-crutis.txt
pcoa_coordinates_nmds_workflow.csv
nmds_coordinates.csv
nmds_summary.csv
NMDS_stressplot.png
NMDS1.png
NMDS2.png
"""

from pathlib import Path
import inspect
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
from sklearn.isotonic import IsotonicRegression
from sklearn.manifold import MDS


# ============================================================
# User settings
# ============================================================
WORK_DIR = Path(".")
# Example:
# WORK_DIR = Path(r"C:/Users/14266/Desktop/meta")

TAX_FILE = WORK_DIR / "tax.txt"
GROUP_FILE = WORK_DIR / "group.txt"

NMDS_DIMENSIONS = 3
NMDS_RESTARTS = 20
NMDS_MAX_ITERATIONS = 1000
ELLIPSE_LEVEL = 0.90
RANDOM_SEED = 20260805

POINT_COLORS = ["#314826", "#ca7442", "#5a558b"]
FILL_COLORS = ["#f2ccac", "#a1d5b9", "#e1abbc"]
MARKERS = ["o", "s", "^", "D", "P", "X"]

# Configure all plot text as Times New Roman.
# On Windows, Times New Roman is normally installed by default.
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
# Data loading and preprocessing
# ============================================================
def load_data(tax_file: Path, group_file: Path):
    """Load and align abundance and grouping tables."""
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

    # Preserve the sample order in tax.txt.
    group = group.set_index("sample").reindex(tax.columns)
    group.index.name = "sample"
    group = group.reset_index()

    return tax, group


def relative_abundance_by_sample(tax: pd.DataFrame):
    """
    Match the R preprocessing:
        data <- data / apply(data, 2, sum)
        data <- t(data)
    """
    relative = tax.div(tax.sum(axis=0), axis=1)
    return relative.T


def bray_curtis_matrix(sample_by_taxon: pd.DataFrame):
    """Calculate the square Bray-Curtis distance matrix."""
    condensed = pdist(
        sample_by_taxon.to_numpy(dtype=float),
        metric="braycurtis",
    )
    matrix = squareform(condensed)

    return pd.DataFrame(
        matrix,
        index=sample_by_taxon.index,
        columns=sample_by_taxon.index,
    )


# ============================================================
# Classical PCoA retained from the original workflow
# ============================================================
def classical_pcoa(distance_matrix: pd.DataFrame, dimensions: int = 3):
    """Classical multidimensional scaling, equivalent to cmdscale()."""
    distances = distance_matrix.to_numpy(dtype=float)
    n = distances.shape[0]

    centering = np.eye(n) - np.ones((n, n)) / n
    b_matrix = -0.5 * centering @ np.square(distances) @ centering

    eigenvalues, eigenvectors = np.linalg.eigh(b_matrix)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = eigenvalues > np.finfo(float).eps
    positive_eigenvalues = eigenvalues[positive]
    positive_eigenvectors = eigenvectors[:, positive]

    if len(positive_eigenvalues) < dimensions:
        raise ValueError(
            f"Only {len(positive_eigenvalues)} positive PCoA axes are available."
        )

    coordinates = (
        positive_eigenvectors[:, :dimensions]
        * np.sqrt(positive_eigenvalues[:dimensions])
    )

    columns = [f"PCoA{i}" for i in range(1, dimensions + 1)]
    result = pd.DataFrame(
        coordinates,
        index=distance_matrix.index,
        columns=columns,
    )
    result.index.name = "sample"
    return result


# ============================================================
# Non-metric multidimensional scaling
# ============================================================
def fit_nmds(distance_matrix: pd.DataFrame):
    """
    Fit non-metric MDS using scikit-learn's SMACOF implementation.

    The helper supports both newer and older scikit-learn APIs.
    Stress is requested as normalized Kruskal Stress-1 where supported.
    """
    signature = inspect.signature(MDS)
    parameters = signature.parameters

    common_arguments = {
        "n_components": NMDS_DIMENSIONS,
        "n_init": NMDS_RESTARTS,
        "max_iter": NMDS_MAX_ITERATIONS,
        "eps": 1e-6,
        "random_state": RANDOM_SEED,
        "n_jobs": -1,
    }

    if "metric_mds" in parameters:
        common_arguments.update({
            "metric_mds": False,
            "metric": "precomputed",
            # Random multi-start initialization is closer to vegan::metaMDS()
            # and preserves the requested number of dimensions.
            "init": "random",
        })
    else:
        common_arguments.update({
            "metric": False,
            "dissimilarity": "precomputed",
        })

    if "normalized_stress" in parameters:
        common_arguments["normalized_stress"] = True

    model = MDS(**common_arguments)
    coordinates = model.fit_transform(distance_matrix.to_numpy(dtype=float))

    coordinate_columns = [
        f"MDS{i}" for i in range(1, NMDS_DIMENSIONS + 1)
    ]
    coordinates_df = pd.DataFrame(
        coordinates,
        index=distance_matrix.index,
        columns=coordinate_columns,
    )
    coordinates_df.index.name = "sample"

    stress = float(model.stress_)

    # For old versions returning raw stress, calculate Stress-1 ourselves.
    if "normalized_stress" not in parameters:
        ordination_distances = pdist(coordinates)
        observed = squareform(
            distance_matrix.to_numpy(dtype=float),
            checks=False,
        )
        isotonic = IsotonicRegression(increasing=True, out_of_bounds="clip")
        fitted_disparities = isotonic.fit_transform(
            observed,
            ordination_distances,
        )
        numerator = np.sum(
            np.square(ordination_distances - fitted_disparities)
        )
        denominator = np.sum(np.square(ordination_distances))
        stress = float(np.sqrt(numerator / denominator))

    return coordinates_df, stress


# ============================================================
# Shepard/stress plot
# ============================================================
def stressplot_data(distance_matrix, nmds_coordinates):
    """
    Prepare observed dissimilarities, ordination distances and monotonic fit.
    """
    observed = squareform(
        distance_matrix.to_numpy(dtype=float),
        checks=False,
    )
    ordination_distances = pdist(
        nmds_coordinates.to_numpy(dtype=float),
        metric="euclidean",
    )

    sorter = np.argsort(observed)
    observed_sorted = observed[sorter]
    ordination_sorted = ordination_distances[sorter]

    isotonic = IsotonicRegression(increasing=True, out_of_bounds="clip")
    fitted_sorted = isotonic.fit_transform(
        observed_sorted,
        ordination_sorted,
    )

    residual_ss = np.sum(
        np.square(ordination_sorted - fitted_sorted)
    )
    total_ss = np.sum(
        np.square(ordination_sorted - ordination_sorted.mean())
    )
    nonmetric_r2 = (
        1.0 - residual_ss / total_ss
        if total_ss > 0
        else np.nan
    )

    linear_correlation = np.corrcoef(
        observed_sorted,
        ordination_sorted,
    )[0, 1]
    linear_r2 = float(linear_correlation ** 2)

    return (
        observed_sorted,
        ordination_sorted,
        fitted_sorted,
        float(nonmetric_r2),
        linear_r2,
    )


def save_stressplot(
    distance_matrix,
    nmds_coordinates,
    stress,
    output_file,
):
    """Create a Shepard diagram comparable to vegan::stressplot()."""
    (
        observed,
        ordination,
        fitted,
        nonmetric_r2,
        linear_r2,
    ) = stressplot_data(distance_matrix, nmds_coordinates)

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        observed,
        ordination,
        s=18,
        alpha=0.60,
        color="#4c4c4c",
        edgecolors="none",
        label="Sample pairs",
    )
    ax.plot(
        observed,
        fitted,
        linewidth=1.5,
        color="#972b1d",
        label="Monotonic fit",
    )

    ax.set_xlabel("Observed Bray-Curtis dissimilarity")
    ax.set_ylabel("Ordination distance")
    ax.set_title("NMDS Shepard Plot")
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.9)

    annotation = (
        f"Stress-1 = {stress:.3f}\n"
        f"Non-metric R² = {nonmetric_r2:.3f}\n"
        f"Linear R² = {linear_r2:.3f}"
    )
    ax.text(
        0.03,
        0.97,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    ax.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)

    return nonmetric_r2, linear_r2


# ============================================================
# Plot helpers
# ============================================================
def build_group_styles(group_order):
    styles = {}

    for i, group_name in enumerate(group_order):
        styles[group_name] = {
            "point_color": POINT_COLORS[i % len(POINT_COLORS)],
            "fill_color": FILL_COLORS[i % len(FILL_COLORS)],
            "marker": MARKERS[i % len(MARKERS)],
        }

    return styles


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


def style_main_axis(ax):
    ax.set_xlabel("NMDS 1")
    ax.set_ylabel("NMDS 2")
    ax.axhline(
        0,
        color="#BEBEBE",
        linestyle="--",
        linewidth=0.8,
        zorder=0,
    )
    ax.axvline(
        0,
        color="#BEBEBE",
        linestyle="--",
        linewidth=0.8,
        zorder=0,
    )
    ax.grid(False)
    ax.tick_params(colors="black")

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.9)


def scatter_groups(ax, nmds_result, group_order, styles):
    for group_name in group_order:
        subset = nmds_result.loc[nmds_result["group"] == group_name]
        style = styles[group_name]

        ax.scatter(
            subset["MDS1"],
            subset["MDS2"],
            s=65,
            color=style["point_color"],
            marker=style["marker"],
            edgecolors="none",
            zorder=3,
        )


def add_confidence_ellipses(
    ax,
    nmds_result,
    group_order,
    styles,
    level=0.90,
):
    """
    Add covariance-based confidence ellipses, comparable to stat_ellipse().
    """
    scale = np.sqrt(chi2.ppf(level, df=2))

    for group_name in group_order:
        subset = nmds_result.loc[
            nmds_result["group"] == group_name,
            ["MDS1", "MDS2"],
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


def save_basic_nmds(nmds_result, stress, output_file):
    """Save the basic two-dimensional NMDS scatter plot."""
    group_order = list(dict.fromkeys(nmds_result["group"]))
    styles = build_group_styles(group_order)

    fig, ax = plt.subplots(figsize=(5, 5))

    scatter_groups(ax, nmds_result, group_order, styles)
    style_main_axis(ax)

    ax.legend(
        handles=make_legend_handles(group_order, styles),
        title=None,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
    )

    fig.text(
        0.5,
        0.015,
        f"Stress = {stress:.3f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


def colored_boxplot(ax, values_by_group, positions, fill_colors, vertical=True):
    """Draw colored grouped marginal boxplots."""
    boxplot = ax.boxplot(
        values_by_group,
        positions=positions,
        widths=0.55,
        vert=vertical,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.0},
        whiskerprops={"color": "black", "linewidth": 0.8},
        capprops={"color": "black", "linewidth": 0.8},
    )

    for patch, color in zip(boxplot["boxes"], fill_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.7)


def save_enhanced_nmds(nmds_result, stress, output_file):
    """
    Save NMDS with:
      - 90% confidence ellipses
      - grouped marginal boxplots on both axes
    """
    group_order = list(dict.fromkeys(nmds_result["group"]))
    styles = build_group_styles(group_order)

    fig = plt.figure(figsize=(5, 5))
    grid = GridSpec(
        2,
        2,
        width_ratios=[4.0, 1.25],
        height_ratios=[1.25, 4.0],
        hspace=0.05,
        wspace=0.05,
    )

    ax_top = fig.add_subplot(grid[0, 0])
    ax_main = fig.add_subplot(grid[1, 0])
    ax_right = fig.add_subplot(grid[1, 1])
    ax_empty = fig.add_subplot(grid[0, 1])
    ax_empty.axis("off")

    scatter_groups(ax_main, nmds_result, group_order, styles)
    add_confidence_ellipses(
        ax_main,
        nmds_result,
        group_order,
        styles,
        level=ELLIPSE_LEVEL,
    )
    style_main_axis(ax_main)

    ax_main.legend(
        handles=make_legend_handles(group_order, styles),
        title=None,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
    )

    x_values = [
        nmds_result.loc[
            nmds_result["group"] == group_name,
            "MDS1",
        ].dropna().to_numpy()
        for group_name in group_order
    ]
    y_values = [
        nmds_result.loc[
            nmds_result["group"] == group_name,
            "MDS2",
        ].dropna().to_numpy()
        for group_name in group_order
    ]

    fill_colors = [
        styles[group_name]["fill_color"]
        for group_name in group_order
    ]
    positions = np.arange(1, len(group_order) + 1)

    colored_boxplot(
        ax_top,
        x_values,
        positions,
        fill_colors,
        vertical=False,
    )
    colored_boxplot(
        ax_right,
        y_values,
        positions,
        fill_colors,
        vertical=True,
    )

    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_right.set_xticks([])
    ax_right.set_yticks([])

    for axis in (ax_top, ax_right):
        axis.grid(False)
        for spine in axis.spines.values():
            spine.set_visible(False)

    fig.text(
        0.44,
        0.008,
        f"Stress = {stress:.3f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )

    fig.subplots_adjust(
        left=0.14,
        right=0.97,
        bottom=0.13,
        top=0.97,
    )
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main workflow
# ============================================================
def main():
    tax, group = load_data(TAX_FILE, GROUP_FILE)

    sample_by_taxon = relative_abundance_by_sample(tax)
    bray = bray_curtis_matrix(sample_by_taxon)

    # Preserve the filename spelling used by the R code.
    bray.to_csv(
        WORK_DIR / "bray-crutis.txt",
        sep="\t",
        float_format="%.10f",
    )

    # Retain the cmdscale step from the original R workflow.
    pcoa_coordinates = classical_pcoa(bray, dimensions=3)
    pcoa_coordinates.reset_index().to_csv(
        WORK_DIR / "pcoa_coordinates_nmds_workflow.csv",
        index=False,
        encoding="utf-8-sig",
    )

    nmds_coordinates, stress = fit_nmds(bray)

    nmds_result = (
        nmds_coordinates.reset_index()
        .merge(group[["sample", "group"]], on="sample", how="left")
    )
    nmds_result.to_csv(
        WORK_DIR / "nmds_coordinates.csv",
        index=False,
        encoding="utf-8-sig",
    )

    nonmetric_r2, linear_r2 = save_stressplot(
        bray,
        nmds_coordinates,
        stress,
        WORK_DIR / "NMDS_stressplot.png",
    )

    summary = pd.DataFrame([{
        "dimensions": NMDS_DIMENSIONS,
        "stress_1": stress,
        "nonmetric_R2": nonmetric_r2,
        "linear_R2": linear_r2,
        "restarts": NMDS_RESTARTS,
        "max_iterations": NMDS_MAX_ITERATIONS,
        "random_seed": RANDOM_SEED,
    }])
    summary.to_csv(
        WORK_DIR / "nmds_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    save_basic_nmds(
        nmds_result,
        stress,
        WORK_DIR / "NMDS1.png",
    )
    save_enhanced_nmds(
        nmds_result,
        stress,
        WORK_DIR / "NMDS2.png",
    )

    interpretation = (
        "good"
        if stress <= 0.10
        else "acceptable"
        if stress <= 0.20
        else "weak"
    )

    print("Analysis completed.")
    print(f"Samples: {len(nmds_result)}")
    print(f"Groups: {', '.join(nmds_result['group'].drop_duplicates())}")
    print(f"NMDS dimensions: {NMDS_DIMENSIONS}")
    print(f"Stress-1: {stress:.4f} ({interpretation})")
    print(f"Non-metric R2: {nonmetric_r2:.4f}")
    print(f"Linear R2: {linear_r2:.4f}")
    print("Generated:")
    print("  bray-crutis.txt")
    print("  pcoa_coordinates_nmds_workflow.csv")
    print("  nmds_coordinates.csv")
    print("  nmds_summary.csv")
    print("  NMDS_stressplot.png")
    print("  NMDS1.png")
    print("  NMDS2.png")


if __name__ == "__main__":
    main()
