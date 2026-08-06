# -*- coding: utf-8 -*-
"""
OTU differential-abundance analysis and Manhattan plots.

This script follows the supplied R/DESeq2 workflow. It uses PyDESeq2 when that
package is installed. If PyDESeq2 is unavailable, it falls back to a documented
negative-binomial GLM approximation so the complete example remains runnable.

Input
-----
OTU_example.txt:
    Rows = OTUs
    Sample columns = raw integer counts
    Last column = semicolon-separated taxonomy string

metadata_example.txt:
    Required columns = sample, Group

Outputs
-------
control_treat.DESeq2.txt
control_treat.txt
control_treat.DESeq2.select.txt
control_treat.DESeq2.up.txt
control_treat.DESeq2.down.txt
DEOTU_CON_Day_28_vs_CON_Day_0.csv
OTU_CON_Day_28_vs_CON_Day_0.csv
man_otu.png
man_phylum.png
analysis_summary.txt
"""

from pathlib import Path
import importlib.util
import logging
import re
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# ============================================================
# User settings
# ============================================================
WORK_DIR = Path("./outdir")
# For your original directory:
# WORK_DIR = Path(r"C:/Users/14266/Desktop/OTU")

OTU_FILE = "./OTU_example.txt"
METADATA_FILE = "./metadata_example.txt"

SAMPLE_COLUMN = "sample"
GROUP_COLUMN = "Group"

# This matches the supplied R contrast:
# results(dds, contrast=c("Group", "CON Day 28", "CON Day 0"))
NUMERATOR_GROUP = "CON Day 28"
DENOMINATOR_GROUP = "CON Day 0"

RELATIVE_ABUNDANCE_FILTER = 0.0005
OUTPUT_ABUNDANCE_FILTER = 0.0001
PADJ_CUTOFF = 0.01
LOG2FC_CUTOFF = 1.0
PLOT_PVALUE_CUTOFF = 0.05
NEGLOG10_CAP = 30.0
RANDOM_SEED = 20260805

# auto: use PyDESeq2 if available, otherwise use the fallback NB-GLM.
# pydeseq2: require PyDESeq2 and raise an error if it is absent.
# fallback: always use the pure-Python NB-GLM approximation.
ANALYSIS_BACKEND = "auto"

# "ratio" matches the default DESeq2 normalization more closely.
# For extremely sparse microbiome tables, "poscounts" can be substituted.
PYDESEQ2_SIZE_FACTOR_METHOD = "ratio"

DARK2_COLORS = [
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#666666",
    "#1F78B4",
    "#B15928",
]

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=RuntimeWarning)

mpl.rcParams.update({
    "font.family": ["Times New Roman", "DejaVu Serif"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 8,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "axes.unicode_minus": False,
})


# ============================================================
# Input handling
# ============================================================
def load_inputs(otu_file: Path, metadata_file: Path):
    otu_raw = pd.read_csv(otu_file, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_file, sep="\t", dtype=str)

    required_metadata = {SAMPLE_COLUMN, GROUP_COLUMN}
    missing = required_metadata.difference(metadata.columns)
    if missing:
        raise ValueError(
            f"metadata is missing required column(s): {sorted(missing)}"
        )

    if metadata[SAMPLE_COLUMN].duplicated().any():
        duplicates = metadata.loc[
            metadata[SAMPLE_COLUMN].duplicated(),
            SAMPLE_COLUMN,
        ].tolist()
        raise ValueError(f"Duplicate metadata sample names: {duplicates}")

    if otu_raw.shape[1] < 2:
        raise ValueError(
            "OTU table must contain at least one sample column and one taxonomy column."
        )

    taxonomy_column = otu_raw.columns[-1]
    taxonomy = otu_raw[[taxonomy_column]].copy()
    taxonomy.columns = ["tax"]

    counts = otu_raw.iloc[:, :-1].apply(pd.to_numeric, errors="raise")

    if (counts < 0).any().any():
        raise ValueError("OTU counts contain negative values.")

    if not np.allclose(counts.to_numpy(), np.round(counts.to_numpy())):
        raise ValueError("DESeq-style analysis requires raw integer counts.")

    counts = counts.astype(int)

    count_samples = set(counts.columns)
    metadata_samples = set(metadata[SAMPLE_COLUMN])

    only_counts = sorted(count_samples - metadata_samples)
    only_metadata = sorted(metadata_samples - count_samples)

    if only_counts:
        raise ValueError(
            "Samples present in OTU table but absent from metadata: "
            + ", ".join(only_counts)
        )
    if only_metadata:
        raise ValueError(
            "Samples present in metadata but absent from OTU table: "
            + ", ".join(only_metadata)
        )

    metadata = metadata.set_index(SAMPLE_COLUMN).reindex(counts.columns)
    metadata.index.name = SAMPLE_COLUMN

    group_values = set(metadata[GROUP_COLUMN])
    required_groups = {NUMERATOR_GROUP, DENOMINATOR_GROUP}
    if not required_groups.issubset(group_values):
        raise ValueError(
            f"Required groups {sorted(required_groups)} were not both found. "
            f"Available groups: {sorted(group_values)}"
        )

    taxonomy["OTU"] = taxonomy.index.astype(str)
    taxonomy["phylum"] = taxonomy["tax"].map(extract_phylum)

    return counts, metadata, taxonomy


def extract_phylum(value):
    """Extract the p__ taxonomic level from a semicolon-delimited string."""
    if pd.isna(value):
        return "Unclassified"

    match = re.search(r"(?:^|;\s*)p__([^;]+)", str(value))
    if match:
        phylum = match.group(1).strip()
        return phylum if phylum else "Unclassified"

    return "Unclassified"


def filter_otus(counts: pd.DataFrame):
    """
    Keep an OTU when its relative abundance exceeds the threshold in at least
    one sample, matching:
        rowSums(otu_relative > threshold) >= 1
    """
    sample_totals = counts.sum(axis=0)
    if (sample_totals == 0).any():
        zero_samples = sample_totals.index[sample_totals == 0].tolist()
        raise ValueError(
            "Samples with zero total counts: " + ", ".join(zero_samples)
        )

    relative = counts.div(sample_totals, axis=1)
    keep = (relative > RELATIVE_ABUNDANCE_FILTER).sum(axis=1) >= 1

    filtered_counts = counts.loc[keep].copy()
    filtered_relative = relative.loc[keep].copy()

    if filtered_counts.empty:
        raise ValueError("No OTUs passed the relative-abundance filter.")

    return filtered_counts, filtered_relative


# ============================================================
# Differential abundance backends
# ============================================================
def run_pydeseq2(counts: pd.DataFrame, metadata: pd.DataFrame):
    """Run the current PyDESeq2 single-factor Wald-test workflow."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    counts_samples_by_otu = counts.T.astype(int)

    # Use a simple internal factor name so formulas and contrasts remain robust
    # even when the original metadata column or group levels contain spaces.
    metadata_for_model = metadata[[GROUP_COLUMN]].rename(
        columns={GROUP_COLUMN: "condition"}
    )
    metadata_for_model["condition"] = pd.Categorical(
        metadata_for_model["condition"],
        categories=[DENOMINATOR_GROUP, NUMERATOR_GROUP],
        ordered=True,
    )

    inference = DefaultInference(n_cpus=1)

    dds = DeseqDataSet(
        counts=counts_samples_by_otu,
        metadata=metadata_for_model,
        design="~condition",
        refit_cooks=True,
        size_factors_fit_type=PYDESEQ2_SIZE_FACTOR_METHOD,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["condition", NUMERATOR_GROUP, DENOMINATOR_GROUP],
        inference=inference,
        alpha=PADJ_CUTOFF,
        cooks_filter=True,
        independent_filter=True,
        quiet=True,
    )
    stats.summary()

    result = stats.results_df.copy()
    result.index = result.index.astype(str)
    result.index.name = "OTU"

    required = {
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    }
    missing = required.difference(result.columns)
    if missing:
        raise RuntimeError(
            f"PyDESeq2 result is missing columns: {sorted(missing)}"
        )

    return result, "PyDESeq2"


def estimate_poscounts_size_factors(counts_samples_by_otu: pd.DataFrame):
    """
    Estimate DESeq2-style poscounts size factors.

    Geometric means are calculated from positive counts only. Each sample's
    factor is the median ratio over OTUs with positive counts in that sample.
    """
    values = counts_samples_by_otu.to_numpy(dtype=float)

    geometric_means = np.empty(values.shape[1], dtype=float)
    for j in range(values.shape[1]):
        positive = values[:, j] > 0
        if positive.any():
            geometric_means[j] = np.exp(np.mean(np.log(values[positive, j])))
        else:
            geometric_means[j] = np.nan

    size_factors = []
    for i in range(values.shape[0]):
        valid = (
            (values[i, :] > 0)
            & np.isfinite(geometric_means)
            & (geometric_means > 0)
        )
        ratios = values[i, valid] / geometric_means[valid]

        if len(ratios) == 0:
            size_factors.append(np.nan)
        else:
            size_factors.append(np.median(ratios))

    size_factors = np.asarray(size_factors, dtype=float)

    if np.any(~np.isfinite(size_factors)) or np.any(size_factors <= 0):
        library_sizes = values.sum(axis=1)
        positive_library = library_sizes[library_sizes > 0]
        geometric_library = np.exp(np.mean(np.log(positive_library)))
        size_factors = library_sizes / geometric_library

    # Center factors to geometric mean 1.
    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    return size_factors


def benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full_like(pvalues, np.nan)

    valid = np.isfinite(pvalues)
    if valid.any():
        adjusted[valid] = multipletests(
            pvalues[valid],
            method="fdr_bh",
        )[1]

    return adjusted


def run_fallback_nb_glm(counts: pd.DataFrame, metadata: pd.DataFrame):
    """
    Pure-Python negative-binomial GLM approximation.

    It uses:
      - poscounts median-ratio size factors,
      - one NB GLM per OTU,
      - a method-of-moments dispersion estimate shrunk toward the global median,
      - Wald tests,
      - Benjamini-Hochberg adjustment.

    It is useful for testing and plotting but is not numerically identical to
    R DESeq2 because it does not reproduce DESeq2's full dispersion-MAP and
    independent-filtering implementation.
    """
    count_matrix = counts.T.astype(float)
    size_factors = estimate_poscounts_size_factors(count_matrix)

    group_indicator = (
        metadata[GROUP_COLUMN].to_numpy() == NUMERATOR_GROUP
    ).astype(float)
    design = sm.add_constant(group_indicator)
    offset = np.log(size_factors)

    normalized = count_matrix.div(size_factors, axis=0)
    norm_means = normalized.mean(axis=0)
    norm_vars = normalized.var(axis=0, ddof=1)

    raw_alpha = (
        (norm_vars - norm_means)
        / np.maximum(np.square(norm_means), 1e-12)
    )
    raw_alpha = raw_alpha.clip(lower=1e-8, upper=10.0)

    finite_positive = raw_alpha[
        np.isfinite(raw_alpha) & (raw_alpha > 0)
    ]
    global_alpha = (
        float(np.median(finite_positive))
        if len(finite_positive)
        else 0.1
    )

    # Mild geometric shrinkage toward a shared dispersion.
    shrunk_alpha = np.sqrt(raw_alpha * global_alpha).clip(1e-8, 10.0)

    records = []

    for otu in count_matrix.columns:
        y = count_matrix[otu].to_numpy(dtype=float)
        alpha = float(shrunk_alpha.loc[otu])

        record = {
            "OTU": str(otu),
            "baseMean": float(norm_means.loc[otu]),
            "log2FoldChange": np.nan,
            "lfcSE": np.nan,
            "stat": np.nan,
            "pvalue": np.nan,
            "dispersion": alpha,
        }

        try:
            model = sm.GLM(
                y,
                design,
                family=sm.families.NegativeBinomial(alpha=alpha),
                offset=offset,
            )
            fit = model.fit(maxiter=200, disp=0)

            coefficient = float(fit.params[1])
            standard_error = float(fit.bse[1])

            record["log2FoldChange"] = coefficient / np.log(2)
            record["lfcSE"] = standard_error / np.log(2)
            record["stat"] = coefficient / standard_error
            record["pvalue"] = 2 * norm.sf(abs(record["stat"]))
        except Exception:
            numerator_values = normalized.loc[
                metadata[GROUP_COLUMN] == NUMERATOR_GROUP,
                otu,
            ]
            denominator_values = normalized.loc[
                metadata[GROUP_COLUMN] == DENOMINATOR_GROUP,
                otu,
            ]

            numerator_mean = float(numerator_values.mean())
            denominator_mean = float(denominator_values.mean())
            record["log2FoldChange"] = np.log2(
                (numerator_mean + 0.5) / (denominator_mean + 0.5)
            )

        records.append(record)

    result = pd.DataFrame(records).set_index("OTU")
    result["padj"] = benjamini_hochberg(result["pvalue"])
    return result, "Fallback NB-GLM approximation"


def differential_analysis(counts: pd.DataFrame, metadata: pd.DataFrame):
    has_pydeseq2 = importlib.util.find_spec("pydeseq2") is not None

    if ANALYSIS_BACKEND == "pydeseq2" and not has_pydeseq2:
        raise ModuleNotFoundError(
            "PyDESeq2 is not installed. Run: pip install pydeseq2"
        )

    if ANALYSIS_BACKEND in {"auto", "pydeseq2"} and has_pydeseq2:
        try:
            return run_pydeseq2(counts, metadata)
        except Exception as error:
            if ANALYSIS_BACKEND == "pydeseq2":
                raise
            print(
                "PyDESeq2 could not run; using fallback backend.\n"
                f"Reason: {error}"
            )

    return run_fallback_nb_glm(counts, metadata)


# ============================================================
# Result processing
# ============================================================
def group_mean_relative_abundance(
    relative: pd.DataFrame,
    metadata: pd.DataFrame,
):
    group_means = {}

    for group_name in [DENOMINATOR_GROUP, NUMERATOR_GROUP]:
        samples = metadata.index[
            metadata[GROUP_COLUMN] == group_name
        ]
        group_means[group_name] = relative[samples].mean(axis=1)

    abundance = pd.DataFrame(group_means)
    abundance["Abundance"] = abundance.mean(axis=1)
    abundance["Abundance_percent"] = abundance["Abundance"] * 100
    return abundance


def classify_significance(result: pd.DataFrame):
    result = result.copy()

    result["sig"] = "none"
    result.loc[
        (result["log2FoldChange"] >= LOG2FC_CUTOFF)
        & (result["padj"] < PADJ_CUTOFF),
        "sig",
    ] = "up"
    result.loc[
        (result["log2FoldChange"] <= -LOG2FC_CUTOFF)
        & (result["padj"] < PADJ_CUTOFF),
        "sig",
    ] = "down"

    return result


def prepare_manhattan_data(
    result,
    abundance,
    taxonomy,
):
    data = (
        result
        .join(abundance[["Abundance", "Abundance_percent"]], how="left")
        .join(
            taxonomy.set_index("OTU")[["phylum"]],
            how="left",
        )
    )

    data.index.name = "otu"
    data = data.reset_index()

    data = data.loc[
        data["Abundance"] > OUTPUT_ABUNDANCE_FILTER
    ].copy()

    data["phylum"] = data["phylum"].fillna("Unclassified")
    safe_p = data["pvalue"].clip(lower=np.nextafter(0, 1))
    data["neglogp"] = -np.log10(safe_p)
    data["neglogp_plot"] = data["neglogp"].clip(upper=NEGLOG10_CAP)

    # Preserve the supplied R meaning:
    # negative LFC -> "enriched"; positive LFC -> "depleted".
    data["level"] = np.select(
        [
            data["pvalue"] >= PLOT_PVALUE_CUTOFF,
            (
                (data["pvalue"] < PLOT_PVALUE_CUTOFF)
                & (data["log2FoldChange"] < 0)
            ),
        ],
        ["nosig", "enriched"],
        default="depleted",
    )

    return data


def safe_filename_label(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def write_outputs(result, manhattan_data):
    result_sorted_p = result.sort_values("pvalue", na_position="last")
    result_sorted_p.to_csv(
        WORK_DIR / "control_treat.DESeq2.txt",
        sep="\t",
        na_rep="NA",
    )

    result_sig = classify_significance(result)
    result_sig = result_sig.sort_values(
        ["padj", "log2FoldChange"],
        ascending=[True, False],
        na_position="last",
    )

    result_sig.to_csv(
        WORK_DIR / "control_treat.txt",
        sep="\t",
        na_rep="NA",
    )

    selected = result_sig.loc[result_sig["sig"].isin(["up", "down"])]
    selected.to_csv(
        WORK_DIR / "control_treat.DESeq2.select.txt",
        sep="\t",
        na_rep="NA",
    )
    result_sig.loc[result_sig["sig"] == "up"].to_csv(
        WORK_DIR / "control_treat.DESeq2.up.txt",
        sep="\t",
        na_rep="NA",
    )
    result_sig.loc[result_sig["sig"] == "down"].to_csv(
        WORK_DIR / "control_treat.DESeq2.down.txt",
        sep="\t",
        na_rep="NA",
    )

    diff = result_sig.loc[
        (result_sig["padj"] < PADJ_CUTOFF)
        & (result_sig["log2FoldChange"].abs() > LOG2FC_CUTOFF)
    ]

    comparison = (
        f"{safe_filename_label(NUMERATOR_GROUP)}"
        f"_vs_"
        f"{safe_filename_label(DENOMINATOR_GROUP)}"
    )

    diff.to_csv(
        WORK_DIR / f"DEOTU_{comparison}.csv",
        encoding="utf-8-sig",
    )
    manhattan_data.to_csv(
        WORK_DIR / f"OTU_{comparison}.csv",
        index=False,
        encoding="utf-8-sig",
    )

    return result_sig, diff


# ============================================================
# Manhattan plotting
# ============================================================
def point_sizes(abundance_percent):
    values = np.asarray(abundance_percent, dtype=float)
    max_value = max(np.nanpercentile(values, 95), 0.01)
    normalized = np.clip(values / max_value, 0, 1)
    return 18 + normalized * 155


def phylum_colors(phylum_order):
    return {
        phylum: DARK2_COLORS[i % len(DARK2_COLORS)]
        for i, phylum in enumerate(phylum_order)
    }


def add_legends(fig, phylum_order, color_map):
    phylum_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            markerfacecolor=color_map[p],
            markeredgecolor=color_map[p],
            markersize=6,
            label=p,
        )
        for p in phylum_order
    ]

    level_markers = {
        "enriched": "^",
        "depleted": "v",
        "nosig": "o",
    }
    level_handles = [
        Line2D(
            [0], [0],
            marker=marker,
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#666666",
            markersize=6,
            label=level,
        )
        for level, marker in level_markers.items()
    ]

    abundance_values = [0.1, 0.5, 1.0]
    abundance_handles = [
        plt.scatter(
            [],
            [],
            s=point_sizes([value])[0],
            color="#777777",
            alpha=0.6,
            edgecolors="none",
            label=f"{value:g}",
        )
        for value in abundance_values
    ]

    legend1 = fig.legend(
        handles=phylum_handles,
        title="Phylum",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.99),
        ncol=min(4, max(1, len(phylum_handles))),
    )
    legend2 = fig.legend(
        handles=level_handles,
        title="Level",
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.69, 0.99),
        ncol=3,
    )
    legend3 = fig.legend(
        handles=abundance_handles,
        title="Mean abundance (%)",
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        ncol=3,
    )

    fig.add_artist(legend1)
    fig.add_artist(legend2)
    fig.add_artist(legend3)


def plot_manhattan_otu(data: pd.DataFrame, output_file: Path):
    rng = np.random.default_rng(RANDOM_SEED)

    phylum_order = [
        value for value in pd.unique(data["phylum"])
        if value != "Low Abundance"
    ]
    if "Low Abundance" in set(data["phylum"]):
        phylum_order.append("Low Abundance")

    color_map = phylum_colors(phylum_order)
    marker_map = {
        "enriched": "^",
        "depleted": "v",
        "nosig": "o",
    }

    plot_data = data.reset_index(drop=True).copy()
    plot_data["x"] = np.arange(len(plot_data), dtype=float)
    plot_data["x_jitter"] = plot_data["x"] + rng.uniform(
        -0.22, 0.22, len(plot_data)
    )

    fig, ax = plt.subplots(figsize=(9, 4))

    ax.axhline(
        -np.log10(PLOT_PVALUE_CUTOFF),
        linestyle="--",
        linewidth=1,
        color="lightgrey",
        zorder=0,
    )

    for level, marker in marker_map.items():
        subset = plot_data.loc[plot_data["level"] == level]

        ax.scatter(
            subset["x_jitter"],
            subset["neglogp_plot"],
            s=point_sizes(subset["Abundance_percent"]),
            c=subset["phylum"].map(color_map),
            marker=marker,
            alpha=0.60,
            linewidths=0.7,
            edgecolors="white",
        )

    title = (
        f"Differential OTU in {NUMERATOR_GROUP} "
        f"vs {DENOMINATOR_GROUP}"
    )
    ax.set_title(title)
    ax.set_xlabel("OTU")
    ax.set_ylabel("-log10(P)")
    ax.set_xticks([])
    ax.tick_params(axis="x", bottom=False)
    ax.grid(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    add_legends(fig, phylum_order, color_map)

    fig.tight_layout(rect=(0, 0, 1, 0.73))
    fig.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.10,
    )
    plt.close(fig)


def plot_manhattan_phylum(data: pd.DataFrame, output_file: Path):
    rng = np.random.default_rng(RANDOM_SEED)

    phylum_order = [
        value for value in pd.unique(data["phylum"])
        if value != "Low Abundance"
    ]
    if "Low Abundance" in set(data["phylum"]):
        phylum_order.append("Low Abundance")

    color_map = phylum_colors(phylum_order)
    marker_map = {
        "enriched": "^",
        "depleted": "v",
        "nosig": "o",
    }
    x_map = {phylum: i for i, phylum in enumerate(phylum_order)}

    plot_data = data.copy()
    plot_data["x"] = plot_data["phylum"].map(x_map).astype(float)
    plot_data["x_jitter"] = plot_data["x"] + rng.uniform(
        -0.25, 0.25, len(plot_data)
    )

    fig, ax = plt.subplots(figsize=(9, 4))

    ax.axhline(
        -np.log10(PLOT_PVALUE_CUTOFF),
        linestyle="--",
        linewidth=1,
        color="lightgrey",
        zorder=0,
    )

    for level, marker in marker_map.items():
        subset = plot_data.loc[plot_data["level"] == level]

        ax.scatter(
            subset["x_jitter"],
            subset["neglogp_plot"],
            s=point_sizes(subset["Abundance_percent"]),
            c=subset["phylum"].map(color_map),
            marker=marker,
            alpha=0.60,
            linewidths=0.7,
            edgecolors="white",
        )

    title = (
        f"Differential OTU in {NUMERATOR_GROUP} "
        f"vs {DENOMINATOR_GROUP}"
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("-log10(P)")
    ax.set_xticks(range(len(phylum_order)))
    ax.set_xticklabels(
        phylum_order,
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.grid(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    add_legends(fig, phylum_order, color_map)

    fig.tight_layout(rect=(0, 0, 1, 0.73))
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
    counts, metadata, taxonomy = load_inputs(
        OTU_FILE,
        METADATA_FILE,
    )

    filtered_counts, filtered_relative = filter_otus(counts)

    result, backend_name = differential_analysis(
        filtered_counts,
        metadata,
    )

    abundance = group_mean_relative_abundance(
        filtered_relative,
        metadata,
    )

    manhattan_data = prepare_manhattan_data(
        result,
        abundance,
        taxonomy,
    )

    result_sig, diff = write_outputs(
        result,
        manhattan_data,
    )

    plot_manhattan_otu(
        manhattan_data,
        WORK_DIR / "man_otu.png",
    )
    plot_manhattan_phylum(
        manhattan_data,
        WORK_DIR / "man_phylum.png",
    )

    up_count = int((result_sig["sig"] == "up").sum())
    down_count = int((result_sig["sig"] == "down").sum())

    summary_lines = [
        f"Backend: {backend_name}",
        f"Input OTUs: {len(counts)}",
        f"Filtered OTUs: {len(filtered_counts)}",
        f"Samples: {counts.shape[1]}",
        f"Comparison: {NUMERATOR_GROUP} vs {DENOMINATOR_GROUP}",
        f"Significant up: {up_count}",
        f"Significant down: {down_count}",
        f"Selected differential OTUs: {len(diff)}",
        "",
        "Interpretation:",
        (
            f"Positive log2FoldChange = higher in {NUMERATOR_GROUP}; "
            f"negative log2FoldChange = higher in {DENOMINATOR_GROUP}."
        ),
        (
            "For continuity with the supplied R plot, 'enriched' denotes "
            "negative LFC and 'depleted' denotes positive LFC."
        ),
    ]

    (WORK_DIR / "analysis_summary.txt").write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print("\n".join(summary_lines))
    print("\nGenerated man_otu.png and man_phylum.png.")


if __name__ == "__main__":
    main()
