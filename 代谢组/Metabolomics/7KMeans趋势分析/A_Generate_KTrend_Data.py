# ============================================================
# 06A_Generate_KTrend_Data.py
#
# 生成适合 K-Means 趋势分析的模拟代谢组数据
#
# 7个条件:
# T0, T1, T2, T3, T4, T5, T6
#
# 每组4个生物学重复
# ============================================================

import numpy as np
import pandas as pd


np.random.seed(42)


# ============================================================
# 1. 基本设置
# ============================================================

groups = [
    "T0",
    "T1",
    "T2",
    "T3",
    "T4",
    "T5",
    "T6"
]

replicates = 4

n_metabolites = 120


metabolites = [
    f"Met_{i+1}"
    for i in range(
        n_metabolites
    )
]


samples = []

sample_groups = []


for group in groups:

    for rep in range(
        replicates
    ):

        samples.append(
            f"{group}_R{rep+1}"
        )

        sample_groups.append(
            group
        )


# ============================================================
# 2. 四类真实趋势模板
# ============================================================

# Cluster 1：中间出现明显峰值
trend_1 = np.array([
    -1.0,
    -0.5,
     0.3,
     2.0,
     0.4,
    -0.3,
    -0.8
])


# Cluster 2：逐步上升
trend_2 = np.array([
    -1.3,
    -0.9,
    -0.4,
     0.0,
     0.5,
     1.0,
     1.5
])


# Cluster 3：逐步下降
trend_3 = np.array([
     1.5,
     1.0,
     0.5,
     0.0,
    -0.4,
    -0.9,
    -1.4
])


# Cluster 4：双峰
trend_4 = np.array([
    -0.8,
     1.2,
    -0.5,
     1.5,
    -0.7,
     0.8,
    -1.0
])


patterns = [
    trend_1,
    trend_2,
    trend_3,
    trend_4
]


# ============================================================
# 3. 生成表达矩阵
# ============================================================

data = np.zeros(
    (
        n_metabolites,
        len(samples)
    )
)


true_cluster = []


for i in range(
    n_metabolites
):

    cluster_id = (
        i // 30
    )

    pattern = patterns[
        cluster_id
    ]

    true_cluster.append(
        cluster_id + 1
    )

    # 每个代谢物略有不同的动态幅度
    amplitude = np.random.uniform(
        0.8,
        1.4
    )

    baseline = np.random.uniform(
        5,
        9
    )

    column = 0

    for group_idx in range(
        len(groups)
    ):

        for rep in range(
            replicates
        ):

            noise = np.random.normal(
                0,
                0.25
            )

            data[
                i,
                column
            ] = (
                baseline
                +
                amplitude
                * pattern[group_idx]
                +
                noise
            )

            column += 1


# ============================================================
# 4. 保存表达矩阵
# ============================================================

expr_df = pd.DataFrame(
    data,
    index=metabolites,
    columns=samples
)

expr_df.index.name = "Metabolite"


group_series = pd.Series(
    sample_groups,
    index=samples,
    name="Group"
)


true_cluster_df = pd.DataFrame(
    {
        "Metabolite": metabolites,
        "True_cluster": true_cluster
    }
)


expr_df.to_csv(
    "KTrend_Expression_Matrix.csv"
)

group_series.to_csv(
    "KTrend_Sample_Group.csv"
)

true_cluster_df.to_csv(
    "KTrend_True_Cluster.csv",
    index=False
)


print(
    "K-Means trend example data generated."
)

print(
    expr_df.shape
)

print(
    group_series.value_counts()
)
