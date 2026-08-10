# ============================================================
# 06B_KMeans_Trend.py
# SCI-style K-Means Metabolite Trend Analysis
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
import matplotlib
matplotlib.use("Agg")

plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 1,
    "pdf.fonttype": 42
})


# ============================================================
# 1. 数据
# ============================================================

expr_df = pd.read_csv(
    "KTrend_Expression_Matrix.csv",
    index_col=0
)

group_df = pd.read_csv(
    "KTrend_Sample_Group.csv",
    index_col=0
)

group_series = group_df.iloc[:, 0]


# ============================================================
# 2. 分组顺序
# ============================================================

group_order = [
    "T0",
    "T1",
    "T2",
    "T3",
    "T4",
    "T5",
    "T6"
]


# ============================================================
# 3. 每组计算平均表达量
#
# 行 = metabolite
# 列 = T0 ... T6
# ============================================================

group_mean = pd.DataFrame(
    index=expr_df.index
)


for group in group_order:

    samples = group_series[
        group_series == group
    ].index

    group_mean[group] = (
        expr_df[
            samples
        ].mean(axis=1)
    )


# ============================================================
# 4. 每个代谢物做Z-score
# ============================================================

row_mean = group_mean.mean(
    axis=1
)

row_std = group_mean.std(
    axis=1,
    ddof=0
)

row_std = row_std.replace(
    0,
    np.nan
)


z_data = (
    group_mean
    .sub(
        row_mean,
        axis=0
    )
    .div(
        row_std,
        axis=0
    )
    .dropna()
)


# ============================================================
# 5. K-Means
# ============================================================

n_clusters = 4


model = KMeans(
    n_clusters=n_clusters,
    random_state=42,
    n_init=50
)


cluster_labels = model.fit_predict(
    z_data.values
)


cluster_series = pd.Series(
    cluster_labels,
    index=z_data.index
)


# ============================================================
# 6. 根据峰值位置排序cluster
#
# 只是为了让Panel顺序比较稳定
# ============================================================

cluster_info = []


for cluster_id in range(
    n_clusters
):

    ids = cluster_series[
        cluster_series == cluster_id
    ].index

    mean_profile = (
        z_data.loc[
            ids
        ]
        .mean(axis=0)
    )

    peak_position = np.argmax(
        mean_profile.values
    )

    cluster_info.append(
        (
            cluster_id,
            peak_position
        )
    )


sorted_clusters = [
    x[0]
    for x in sorted(
        cluster_info,
        key=lambda x: x[1]
    )
]


# ============================================================
# 7. 图形
# ============================================================

fig, axes = plt.subplots(
    1,
    n_clusters,
    figsize=(12, 4.8),
    sharey=True
)


cluster_colors = [
    "#D73027",
    "#F46D43",
    "#FDAE61",
    "#4575B4"
]


x = np.arange(
    len(group_order)
)


result = []


# ============================================================
# 8. 分Cluster绘制
# ============================================================

for panel_i, cluster_id in enumerate(
    sorted_clusters
):

    ax = axes[
        panel_i
    ]

    metabolite_ids = cluster_series[
        cluster_series == cluster_id
    ].index

    cluster_data = z_data.loc[
        metabolite_ids
    ]

    color = cluster_colors[
        panel_i
    ]


    # --------------------------------------------------------
    # Individual metabolites
    # --------------------------------------------------------

    for metabolite in metabolite_ids:

        ax.plot(
            x,
            cluster_data.loc[
                metabolite
            ].values,
            color=color,
            linewidth=0.55,
            alpha=0.40
        )


    # --------------------------------------------------------
    # Cluster平均趋势
    # --------------------------------------------------------

    mean_profile = (
        cluster_data
        .mean(axis=0)
    )


    ax.plot(
        x,
        mean_profile.values,
        color="black",
        linewidth=2.6,
        zorder=20
    )


    # --------------------------------------------------------
    # 保存cluster结果
    # --------------------------------------------------------

    for metabolite in metabolite_ids:

        result.append({
            "Metabolite":
                metabolite,

            "Sub_class":
                panel_i + 1
        })


    # --------------------------------------------------------
    # Panel标题
    # --------------------------------------------------------

    ax.set_title(
        f"Sub class {panel_i + 1}, "
        f"total: {len(metabolite_ids)}",
        fontweight="bold",
        pad=12
    )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        group_order,
        rotation=90
    )


    ax.axhline(
        0,
        color="gray",
        linewidth=0.5,
        linestyle="--",
        alpha=0.5
    )


    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)


    ax.tick_params(
        direction="out"
    )


    if panel_i == 0:

        ax.set_ylabel(
            "Standardized value"
        )


# ============================================================
# 9. Y轴范围
# ============================================================

max_abs = np.nanmax(
    np.abs(
        z_data.values
    )
)

ylim = np.ceil(
    max_abs * 1.10
)


for ax in axes:

    ax.set_ylim(
        -ylim,
        ylim
    )


# ============================================================
# 10. 保存
# ============================================================

plt.tight_layout()

plt.savefig(
    "KMeans_Trend.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "KMeans_Trend.pdf",
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 11. 保存Cluster分类结果
# ============================================================

result_df = pd.DataFrame(
    result
)

result_df = result_df.sort_values(
    [
        "Sub_class",
        "Metabolite"
    ]
)

result_df.to_csv(
    "KMeans_Cluster_Result.csv",
    index=False
)


print(
    "K-Means trend plot saved."
)

print(
    result_df["Sub_class"]
    .value_counts()
    .sort_index()
)
