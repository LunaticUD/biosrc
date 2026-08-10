# ============================================================
# 05_Clustered_Heatmap.py
# SCI-style Bidirectional Clustered Heatmap
#
# 输入：
#   Differential_Metabolites.csv
#   Expression_Matrix.csv
#   Sample_Group.csv
#
# 输出：
#   05_Clustered_Heatmap.png
#   05_Clustered_Heatmap.pdf
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.patches import Patch
import matplotlib
matplotlib.use("Agg")

# ============================================================
# 1. SCI 全局绘图参数
# ============================================================

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10


# ============================================================
# 2. 双向聚类热图函数
# ============================================================

def plot_clustered_heatmap(
    sig_df,
    expr_df,
    group_series,
    figsize=(10, 9),
    save_path="05_Clustered_Heatmap.png"
):
    """
    绘制差异代谢物双向层次聚类热图。

    Parameters
    ----------
    sig_df : pd.DataFrame
        差异代谢物统计结果。
        index = 代谢物名称。

    expr_df : pd.DataFrame
        表达矩阵。
        行 = 代谢物。
        列 = 样本。

    group_series : pd.Series
        样本分组。
        index = 样本名称。
        value = 分组名称。

    figsize : tuple
        图片大小。

    save_path : str
        PNG 保存路径。
    """

    # ========================================================
    # 1. 匹配差异代谢物
    # ========================================================

    common_ids = sig_df.index.intersection(
        expr_df.index
    )

    if len(common_ids) < 2:
        raise ValueError(
            "可用于聚类分析的差异代谢物少于 2 个。"
        )

    data_to_plot = expr_df.loc[
        common_ids
    ].copy()


    # ========================================================
    # 2. 匹配样本
    # ========================================================

    common_samples = data_to_plot.columns.intersection(
        group_series.index
    )

    if len(common_samples) < 2:
        raise ValueError(
            "表达矩阵与分组信息匹配到的样本少于 2 个。"
        )

    data_to_plot = data_to_plot.loc[
        :,
        common_samples
    ]

    groups = group_series.loc[
        common_samples
    ]


    print(
        f"Heatmap data: "
        f"{data_to_plot.shape[0]} metabolites × "
        f"{data_to_plot.shape[1]} samples"
    )


    # ========================================================
    # 3. 每个代谢物做 Z-score
    #
    # 对每一行：
    # z = (x - mean) / std
    # ========================================================

    row_mean = data_to_plot.mean(
        axis=1
    )

    row_std = data_to_plot.std(
        axis=1,
        ddof=0
    )

    # 防止标准差为 0
    row_std = row_std.replace(
        0,
        np.nan
    )

    zscore_data = data_to_plot.sub(
        row_mean,
        axis=0
    ).div(
        row_std,
        axis=0
    )

    # 删除无法计算 Z-score 的代谢物
    zscore_data = zscore_data.dropna(
        axis=0,
        how="any"
    )


    if zscore_data.shape[0] < 2:
        raise ValueError(
            "Z-score 标准化后剩余代谢物少于 2 个。"
        )


    # ========================================================
    # 4. Group 颜色
    # ========================================================

    unique_groups = list(
        pd.unique(groups)
    )

    # 这里显式设置常用组颜色
    # 如果组数超过这里定义的数量，会自动从 Set2 补充
    base_group_colors = [
        "#1F77B4",
        "#E67E22",
        "#2CA02C",
        "#D62728",
        "#9467BD",
        "#8C564B",
        "#E377C2",
        "#7F7F7F"
    ]

    if len(unique_groups) <= len(base_group_colors):

        group_palette = base_group_colors[
            :len(unique_groups)
        ]

    else:

        group_palette = sns.color_palette(
            "Set2",
            n_colors=len(unique_groups)
        )


    group_colors = {
        group: color
        for group, color in zip(
            unique_groups,
            group_palette
        )
    }


    # 每个样本对应的颜色
    col_colors = groups.map(
        group_colors
    )

    col_colors.name = "Group"


    # ========================================================
    # 5. 绘制 clustermap
    # ========================================================

    g = sns.clustermap(

        zscore_data,

        # -------------------------
        # 热图颜色
        # -------------------------
        cmap="RdYlGn_r",

        center=0,

        vmin=-2.5,
        vmax=2.5,

        # -------------------------
        # 双向聚类
        # -------------------------
        row_cluster=True,
        col_cluster=True,

        method="average",
        metric="euclidean",

        # -------------------------
        # 顶部 Group 色条
        # -------------------------
        col_colors=col_colors,

        # -------------------------
        # 热图样式
        # -------------------------
        linewidths=0,

        figsize=figsize,

        dendrogram_ratio=(
            0.12,
            0.12
        ),

        colors_ratio=0.025,

        cbar_kws={
            "label": "Z-score"
        }
    )


    # ========================================================
    # 6. 热图坐标轴设置
    # ========================================================

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")


    # X 轴样本名称
    g.ax_heatmap.tick_params(
        axis="x",
        labelsize=8,
        rotation=45,
        length=0
    )


    # Y 轴代谢物名称
    g.ax_heatmap.tick_params(
        axis="y",
        labelsize=7,
        rotation=0,
        length=0
    )


    # 设置刻度字体
    for label in g.ax_heatmap.get_xticklabels():

        label.set_fontfamily(
            "Times New Roman"
        )


    for label in g.ax_heatmap.get_yticklabels():

        label.set_fontfamily(
            "Times New Roman"
        )


    # ========================================================
    # 7. 顶部 Group 注释标题
    # ========================================================

    if hasattr(g, "ax_col_colors"):

        g.ax_col_colors.set_ylabel(
            "Group",
            rotation=0,
            ha="left",
            va="center",
            fontsize=9,
            fontfamily="Times New Roman"
        )

        g.ax_col_colors.yaxis.set_label_coords(
            1.01,
            0.5
        )


    # ========================================================
    # 8. 调整 Z-score colorbar
    #
    # 右侧同一列的上半部分
    #
    # set_position:
    # [left, bottom, width, height]
    # ========================================================

    g.cax.set_position(
        [
            0.86,   # left
            0.60,   # bottom
            0.025,  # width
            0.20    # height
        ]
    )


    g.cax.tick_params(
        labelsize=9,
        width=0.8,
        length=3
    )


    g.cax.set_ylabel(
        "Z-score",
        fontsize=10,
        fontfamily="Times New Roman",
        labelpad=8
    )


    for label in g.cax.get_yticklabels():

        label.set_fontfamily(
            "Times New Roman"
        )


    # ========================================================
    # 9. 创建 Group 图例
    # ========================================================

    legend_handles = [

        Patch(
            facecolor=group_colors[group],
            edgecolor="none",
            label=str(group)
        )

        for group in unique_groups
    ]


    # ========================================================
    # 10. Group 图例放在 colorbar 下方
    #
    # 与 colorbar 使用相同 left = 0.86
    # 因而视觉上处于同一列
    # ========================================================

    group_legend = g.fig.legend(

        handles=legend_handles,

        title="Group",

        loc="upper left",

        bbox_to_anchor=(
            0.86,
            0.54
        ),

        frameon=False,

        fontsize=9,

        title_fontsize=10,

        handlelength=1.0,

        handleheight=1.0,

        handletextpad=0.5,

        labelspacing=0.5,

        borderaxespad=0
    )


    # 图例字体统一为 Times New Roman
    plt.setp(
        group_legend.get_texts(),
        fontfamily="Times New Roman"
    )

    plt.setp(
        group_legend.get_title(),
        fontfamily="Times New Roman"
    )


    # ========================================================
    # 11. 标题
    # ========================================================

    g.fig.suptitle(
        " ",
        fontsize=14,
        fontfamily="Times New Roman",
        y=0.98
    )


    # ========================================================
    # 12. 调整热图区域
    #
    # 给右侧 colorbar + legend 留空间
    # ========================================================

    heatmap_pos = g.ax_heatmap.get_position()

    g.ax_heatmap.set_position(
        [
            heatmap_pos.x0,
            heatmap_pos.y0,
            0.70 - heatmap_pos.x0,
            heatmap_pos.height
        ]
    )


    # 行聚类树
    row_dendrogram_pos = (
        g.ax_row_dendrogram.get_position()
    )

    g.ax_row_dendrogram.set_position(
        [
            row_dendrogram_pos.x0,
            row_dendrogram_pos.y0,
            row_dendrogram_pos.width,
            heatmap_pos.height
        ]
    )


    # 列聚类树
    col_dendrogram_pos = (
        g.ax_col_dendrogram.get_position()
    )

    g.ax_col_dendrogram.set_position(
        [
            heatmap_pos.x0,
            col_dendrogram_pos.y0,
            0.70 - heatmap_pos.x0,
            col_dendrogram_pos.height
        ]
    )


    # 顶部 Group 色条同步宽度
    if hasattr(g, "ax_col_colors"):

        col_colors_pos = (
            g.ax_col_colors.get_position()
        )

        g.ax_col_colors.set_position(
            [
                heatmap_pos.x0,
                col_colors_pos.y0,
                0.70 - heatmap_pos.x0,
                col_colors_pos.height
            ]
        )


    # ========================================================
    # 13. 保存 PNG
    # ========================================================

    g.fig.savefig(
        save_path,
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )


    # ========================================================
    # 14. 同时保存 PDF
    # ========================================================

    if save_path.lower().endswith(
        ".png"
    ):

        pdf_path = (
            save_path[:-4]
            + ".pdf"
        )

        g.fig.savefig(
            pdf_path,
            format="pdf",
            bbox_inches="tight",
            facecolor="white"
        )


    plt.close(
        g.fig
    )


    print(
        f"Clustered heatmap saved to: "
        f"{save_path}"
    )


# ============================================================
# 3. 读取前面差异筛选脚本生成的数据
# ============================================================

sig = pd.read_csv(
    "Differential_Metabolites.csv",
    index_col=0
)


expr_df = pd.read_csv(
    "Expression_Matrix.csv",
    index_col=0
)


group_df = pd.read_csv(
    "Sample_Group.csv",
    index_col=0
)


# Sample_Group.csv 通常只有一列 Group
group_series = group_df.iloc[
    :,
    0
]

group_series.name = "Group"


# ============================================================
# 4. 绘图
# ============================================================

plot_clustered_heatmap(
    sig_df=sig,
    expr_df=expr_df,
    group_series=group_series,
    figsize=(10, 9),
    save_path="05_Clustered_Heatmap.png"
)
