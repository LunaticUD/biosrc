import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

# ============================================================
# SCI-style Volcano Plot
# ============================================================

def plot_volcano_sci(
    stats_df,
    fc_thresh=1.0,
    p_thresh=0.05,
    label_top_n=3,
    save_path="Volcano_Plot.png"
):

    # --------------------------------------------------------
    # 1. 数据准备
    # --------------------------------------------------------

    df = stats_df.copy()

    if "Metabolite" not in df.columns:
        df["Metabolite"] = df.index.astype(str)

    df["-Log10P"] = -np.log10(
        df["P_value"].clip(lower=1e-300)
    )


    # --------------------------------------------------------
    # 2. 分类
    # --------------------------------------------------------

    def classify(row):

        if (
            row["P_value"] < p_thresh
            and row["log2FC"] >= fc_thresh
        ):
            return "Up"

        elif (
            row["P_value"] < p_thresh
            and row["log2FC"] <= -fc_thresh
        ):
            return "Down"

        else:
            return "Insignificant"


    df["Group"] = df.apply(
        classify,
        axis=1
    )


    up_n = (df["Group"] == "Up").sum()
    down_n = (df["Group"] == "Down").sum()
    insig_n = (df["Group"] == "Insignificant").sum()


    # --------------------------------------------------------
    # 3. SCI字体设置
    # --------------------------------------------------------

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 11


    # --------------------------------------------------------
    # 4. 建图
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7.5, 6)
    )


    colors = {
        "Up": "#D73027",
        "Down": "#4575B4",
        "Insignificant": "#BDBDBD"
    }


    # --------------------------------------------------------
    # 5. 绘制散点
    # --------------------------------------------------------

    # 先画不显著点，避免覆盖显著点
    for group_name in [
        "Insignificant",
        "Down",
        "Up"
    ]:

        temp = df[
            df["Group"] == group_name
        ]

        ax.scatter(
            temp["log2FC"],
            temp["-Log10P"],
            s=38,
            c=colors[group_name],
            alpha=0.85,
            edgecolors="none",
            zorder=2
        )


    # --------------------------------------------------------
    # 6. 阈值线
    # --------------------------------------------------------

    line_color = "#555555"

    ax.axhline(
        -np.log10(p_thresh),
        linestyle="--",
        linewidth=1.1,
        color=line_color,
        zorder=1
    )

    ax.axvline(
        fc_thresh,
        linestyle="--",
        linewidth=1.1,
        color=line_color,
        zorder=1
    )

    ax.axvline(
        -fc_thresh,
        linestyle="--",
        linewidth=1.1,
        color=line_color,
        zorder=1
    )


    # --------------------------------------------------------
    # 7. 标记最显著代谢物
    # --------------------------------------------------------

    sig_df = df[
        df["Group"].isin(["Up", "Down"])
    ].copy()

    # 按P值选择最显著的几个
    label_df = sig_df.sort_values(
        "P_value"
    ).head(label_top_n)


    # 7. 标记代谢物
    for _, row in label_df.iterrows():

        x_offset = 6 if row["log2FC"] > 0 else -6
        align = "left" if row["log2FC"] > 0 else "right"

        ax.annotate(
            row["Metabolite"],
            xy=(
                row["log2FC"],
                row["-Log10P"]
            ),
            xytext=(
                x_offset,
                5
            ),
            textcoords="offset points",
            ha=align,
            va="bottom",
            fontsize=10,
            fontfamily="Times New Roman",
            color="black"
        )


    # ============================================================
    # 注意！！！从这里开始已经退出 for 循环
    # ============================================================

    # --------------------------------------------------------
    # 8. 图例
    # --------------------------------------------------------

    legend_items = [
        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            label=f"Up ({up_n})",
            markerfacecolor=colors["Up"],
            markeredgecolor="none",
            markersize=8
        ),

        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            label=f"Down ({down_n})",
            markerfacecolor=colors["Down"],
            markeredgecolor="none",
            markersize=8
        ),

        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            label=f"NS ({insig_n})",
            markerfacecolor=colors["Insignificant"],
            markeredgecolor="none",
            markersize=8
        )
    ]


    # --------------------------------------------------------
    # 图例放到坐标轴右侧外面
    # --------------------------------------------------------

    legend = ax.legend(
        handles=legend_items,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.05, 1.00),
        fontsize=11,
        borderaxespad=0,
        handletextpad=0.6,
        labelspacing=0.7
    )


    # --------------------------------------------------------
    # 9. 坐标轴
    # --------------------------------------------------------

    ax.set_xlabel(
        r"$\log_{2}(\mathrm{Fold\ Change})$",
        fontsize=13
    )

    ax.set_ylabel(
        r"$-\log_{10}(P\mathrm{-value})$",
        fontsize=13
    )

    ax.set_title(
        "Volcano Plot of Differential Metabolites",
        fontsize=14,
        pad=12
    )


    # --------------------------------------------------------
    # 10. SCI风格
    # --------------------------------------------------------

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

    ax.tick_params(
        axis="both",
        direction="out",
        width=1,
        length=4,
        labelsize=11
    )

    ax.grid(False)


    # --------------------------------------------------------
    # 11. 给右侧图例留空间
    # --------------------------------------------------------

    fig.subplots_adjust(
        left=0.14,
        right=0.72,
        bottom=0.13,
        top=0.90
    )


    # --------------------------------------------------------
    # 12. 保存
    # --------------------------------------------------------

    fig.savefig(
        save_path,
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )


    if save_path.lower().endswith(".png"):

        pdf_path = save_path[:-4] + ".pdf"

        fig.savefig(
            pdf_path,
            format="pdf",
            bbox_inches="tight",
            facecolor="white"
        )


    plt.close(fig)

    print(f"Volcano plot saved to: {save_path}")

diff_stats = pd.read_csv(
    "All_Metabolite_Statistics.csv",
    index_col=0
)
plot_volcano_sci(
    stats_df=diff_stats,
    fc_thresh=1.0,
    p_thresh=0.05,
    label_top_n=3,
    save_path="Volcano_Plot.png"
)