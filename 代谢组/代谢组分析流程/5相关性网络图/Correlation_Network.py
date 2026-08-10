# ============================================================
# 04_Correlation_Network.py
# SCI-style Metabolite Correlation Network
#
# 节点颜色：代谢物类别
# 边颜色：
#     红色 = 正相关
#     蓝色 = 负相关
# 节点大小：Degree
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.lines import Line2D
import matplotlib
matplotlib.use("Agg")

# ============================================================
# 1. 全局SCI绘图参数
# ============================================================

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10


# ============================================================
# 2. 相关性网络绘图函数
# ============================================================

def plot_correlation_network(
    sig_df,
    expr_df,
    node_category=None,
    threshold=0.75,
    figsize=(12, 10),
    save_path="Correlation_Network.png"
):
    """
    绘制差异代谢物 Pearson 相关性网络图

    参数
    ----------
    sig_df : pd.DataFrame
        差异代谢物统计结果
        index = 代谢物名称

    expr_df : pd.DataFrame
        表达矩阵
        行 = 代谢物
        列 = 样本

    node_category : dict
        代谢物 -> 类别 映射

    threshold : float
        |Pearson r| 阈值
    """

    # ========================================================
    # 1. 匹配代谢物
    # ========================================================

    common_ids = sig_df.index.intersection(
        expr_df.index
    )

    if len(common_ids) < 2:
        print("No enough metabolites for correlation network.")
        return

    sig_expr = expr_df.loc[common_ids].copy()


    # ========================================================
    # 2. Pearson相关矩阵
    # ========================================================

    corr_matrix = sig_expr.T.corr(
        method="pearson"
    )

    print(
        f"Metabolites used: {corr_matrix.shape[0]}"
    )


    # ========================================================
    # 3. 构建网络
    # ========================================================

    G = nx.Graph()

    metabolites = corr_matrix.columns.tolist()

    # 先加入所有节点
    G.add_nodes_from(metabolites)


    for i in range(len(metabolites)):

        for j in range(i + 1, len(metabolites)):

            node1 = metabolites[i]
            node2 = metabolites[j]

            r = corr_matrix.iloc[i, j]

            if pd.isna(r):
                continue

            if abs(r) >= threshold:

                G.add_edge(
                    node1,
                    node2,
                    correlation=r,
                    weight=abs(r)
                )


    print(
        f"Initial network: "
        f"{G.number_of_nodes()} nodes, "
        f"{G.number_of_edges()} edges"
    )


    # ========================================================
    # 4. 删除孤立节点
    # ========================================================

    isolates = list(
        nx.isolates(G)
    )

    G.remove_nodes_from(
        isolates
    )


    print(
        f"After removing isolates: "
        f"{G.number_of_nodes()} nodes, "
        f"{G.number_of_edges()} edges"
    )


    if G.number_of_edges() == 0:

        print(
            "No edges found. "
            "Try threshold=0.6 or threshold=0.7."
        )

        return


    # ========================================================
    # 5. 类别颜色
    #
    # 同一种类永远使用同一种颜色
    # ========================================================

    category_palette = {

        "Organic acid and Its derivatives":
            "#66C2A5",

        "Heterocyclic compounds":
            "#FFF176",

        "Others":
            "#B39DDB",

        "Amino acid and Its metabolites":
            "#E76F51",

        "Aldehyde,Ketones,Esters":
            "#80B1D3",

        "Benzene and substituted derivatives":
            "#FDB462",

        "GP":
            "#8D6E63",

        "Bile acids":
            "#C44E52",

        "Nucleotide and Its metabolites":
            "#4DB6AC",

        "Phenolic acids":
            "#C5C86D",

        "Alkaloids":
            "#BDBDBD",

        "Alcohol and amines":
            "#F4A6A6",

        "FA":
            "#B2DF8A",

        "Tryptamines,Cholines,Pigments":
            "#FF8C00",

        "Lipid":
            "#6A9FB5"
    }


    # ========================================================
    # 6. 节点颜色
    # ========================================================

    node_colors = []

    for node in G.nodes():

        if node_category is None:

            category = "Others"

        else:

            category = node_category.get(
                node,
                "Others"
            )

        node_colors.append(
            category_palette.get(
                category,
                "#BDBDBD"
            )
        )


    # ========================================================
    # 7. 节点大小
    #
    # 根据 degree 设置
    # ========================================================

    degree_dict = dict(
        G.degree()
    )

    node_sizes = [

        180 + degree_dict[node] * 35

        for node in G.nodes()
    ]


    # ========================================================
    # 8. 网络布局
    #
    # spring_layout 比较适合这种高密度网络
    # ========================================================

    pos = nx.spring_layout(
        G,
        seed=42,
        k=1.5 / np.sqrt(
            G.number_of_nodes()
        ),
        iterations=500,
        weight="weight"
    )


    # ========================================================
    # 9. 获取正负相关边
    # ========================================================

    positive_edges = []

    negative_edges = []


    positive_widths = []

    negative_widths = []


    for u, v, d in G.edges(data=True):

        r = d["correlation"]

        width = 0.5 + 1.8 * abs(r)

        if r > 0:

            positive_edges.append(
                (u, v)
            )

            positive_widths.append(
                width
            )

        else:

            negative_edges.append(
                (u, v)
            )

            negative_widths.append(
                width
            )


    # ========================================================
    # 10. 绘图
    # ========================================================

    fig, ax = plt.subplots(
        figsize=figsize
    )


    # --------------------------------------------------------
    # 正相关边
    # --------------------------------------------------------

    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=positive_edges,
        width=positive_widths,
        edge_color="#FF4F87",
        alpha=0.50,
        ax=ax
    )


    # --------------------------------------------------------
    # 负相关边
    # --------------------------------------------------------

    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=negative_edges,
        width=negative_widths,
        edge_color="#0096D6",
        alpha=0.50,
        ax=ax
    )


    # ========================================================
    # 11. 节点
    # ========================================================

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="#666666",
        linewidths=0.6,
        alpha=0.95,
        ax=ax
    )


    # ========================================================
    # 12. 标签
    # ========================================================

    nx.draw_networkx_labels(
        G,
        pos,
        font_family="Times New Roman",
        font_size=6,
        font_color="black",
        ax=ax
    )


    # ========================================================
    # 13. 代谢物类别图例
    # ========================================================

    if node_category is not None:

        used_categories = []

        for node in G.nodes():

            category = node_category.get(
                node,
                "Others"
            )

            if category not in used_categories:

                used_categories.append(
                    category
                )


        category_legend = []

        for category in used_categories:

            color = category_palette.get(
                category,
                "#BDBDBD"
            )

            category_legend.append(

                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="None",
                    markerfacecolor=color,
                    markeredgecolor="#666666",
                    markeredgewidth=0.6,
                    markersize=8,
                    label=category
                )

            )


        legend1 = ax.legend(
            handles=category_legend,
            title="Metabolite Class",
            loc="upper left",
            bbox_to_anchor=(0.70, 1.02),
            frameon=False,
            fontsize=9,
            title_fontsize=10,
            labelspacing=0.5
        )

        ax.add_artist(
            legend1
        )


    # ========================================================
    # 14. 相关性图例
    # ========================================================

    edge_legend = [

        Line2D(
            [0], [0],
            color="#FF4F87",
            lw=2,
            label="Positive correlation"
        ),

        Line2D(
            [0], [0],
            color="#0096D6",
            lw=2,
            label="Negative correlation"
        )
    ]


    ax.legend(
        handles=edge_legend,
        title="Correlation",
        loc="lower left",
        bbox_to_anchor=(0.70, 0.00),
        frameon=False,
        fontsize=9,
        title_fontsize=10
    )


    # ========================================================
    # 15. 标题
    # ========================================================

    ax.set_title(
        f"Metabolite Correlation Network "
        f"(|r| ≥ {threshold})",
        fontsize=14,
        pad=12
    )


    ax.axis(
        "off"
    )


    # 给右侧图例留位置
    fig.subplots_adjust(
        left=0.03,
        right=0.73,
        bottom=0.03,
        top=0.93
    )


    # ========================================================
    # 16. 保存
    # ========================================================

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

    print(
        f"Correlation network saved to: "
        f"{save_path}"
    )




np.random.seed(42)

sig = pd.read_csv(
    "Differential_Metabolites.csv",
    index_col=0
)

expr_df = pd.read_csv(
    "Expression_Matrix.csv",
    index_col=0
)

classes = [
    "Organic acid and Its derivatives",
    "Heterocyclic compounds",
    "Amino acid and Its metabolites",
    "Nucleotide and Its metabolites",
    "Benzene and substituted derivatives",
    "Phenolic acids",
    "Alcohol and amines",
    "FA"
]

category_dict = {
    metabolite: np.random.choice(classes)
    for metabolite in sig.index
}

plot_correlation_network(
    sig,
    expr_df,
    node_category=category_dict,
    threshold=0.75,
    save_path="04_Correlation_Network.png"
)
