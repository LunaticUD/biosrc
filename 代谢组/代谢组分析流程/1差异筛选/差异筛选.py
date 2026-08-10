# ============================================================
# Differential_Metabolite_Screening.py
# 差异代谢物筛选独立脚本
#
# 筛选标准：
# VIP > 1
# |log2FC| >= 1
# P < 0.05
#
# 输出：
# 1. All_Metabolite_Statistics.csv
# 2. Differential_Metabolites.csv
# ============================================================

import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from sklearn.cross_decomposition import PLSRegression


# ============================================================
# 1. 模拟代谢组数据
# ============================================================

np.random.seed(42)

n_metabolites = 200

# 两组，每组 8 个样本
n_control = 8
n_treatment = 8

n_samples = n_control + n_treatment

metabolites = [
    f"Met_{i+1}"
    for i in range(n_metabolites)
]

control_samples = [
    f"Control_{i+1}"
    for i in range(n_control)
]

treatment_samples = [
    f"Treatment_{i+1}"
    for i in range(n_treatment)
]

samples = control_samples + treatment_samples


# ------------------------------------------------------------
# 基础表达量
# ------------------------------------------------------------

data = np.random.lognormal(
    mean=1.5,
    sigma=0.35,
    size=(n_metabolites, n_samples)
)


# ------------------------------------------------------------
# 人为制造差异
#
# Met_1 ~ Met_25：Treatment 上调
# Met_26 ~ Met_45：Treatment 下调
# ------------------------------------------------------------

# 上调
data[0:25, n_control:] *= 2.5

# 下调
data[25:45, n_control:] *= 0.40


# ------------------------------------------------------------
# 构建表达矩阵
#
# 行 = 代谢物
# 列 = 样本
# ------------------------------------------------------------

expr_df = pd.DataFrame(
    data,
    index=metabolites,
    columns=samples
)

expr_df.index.name = "Metabolite"


# ------------------------------------------------------------
# 样本分组信息
# ------------------------------------------------------------

group_series = pd.Series(
    ["Control"] * n_control +
    ["Treatment"] * n_treatment,
    index=samples,
    name="Group"
)


print("表达矩阵大小：", expr_df.shape)

print("\n样本分组：")
print(group_series.value_counts())


# ============================================================
# 2. 差异筛选函数
# ============================================================

def differential_metabolite_screening(
    expr_df,
    group_series,
    control_group="Control",
    treatment_group="Treatment",
    vip_thresh=1.0,
    log2fc_thresh=1.0,
    p_thresh=0.05
):
    """
    差异代谢物筛选

    参数
    ----------
    expr_df : pd.DataFrame
        表达矩阵
        行 = 代谢物
        列 = 样本

    group_series : pd.Series
        样本分组信息
        index = 样本名称
        value = 分组名称

    control_group : str
        对照组名称

    treatment_group : str
        处理组名称

    vip_thresh : float
        VIP 阈值

    log2fc_thresh : float
        |log2FC| 阈值

    p_thresh : float
        P 值阈值

    返回
    ----------
    diff_stats : pd.DataFrame
        所有代谢物统计结果

    sig : pd.DataFrame
        差异代谢物结果
    """

    # --------------------------------------------------------
    # 1. 检查样本是否匹配
    # --------------------------------------------------------

    common_samples = expr_df.columns.intersection(
        group_series.index
    )

    if len(common_samples) == 0:
        raise ValueError(
            "expr_df 与 group_series 没有共同样本。"
        )

    expr_df = expr_df.loc[:, common_samples]

    group_series = group_series.loc[
        common_samples
    ]


    # --------------------------------------------------------
    # 2. 提取两组数据
    # --------------------------------------------------------

    ctrl_samples = group_series[
        group_series == control_group
    ].index

    treat_samples = group_series[
        group_series == treatment_group
    ].index


    if len(ctrl_samples) < 2:
        raise ValueError(
            f"{control_group} 组样本数少于2。"
        )

    if len(treat_samples) < 2:
        raise ValueError(
            f"{treatment_group} 组样本数少于2。"
        )


    ctrl = expr_df.loc[:, ctrl_samples]

    treat = expr_df.loc[:, treat_samples]


    print(
        f"\n{control_group} 样本数："
        f"{len(ctrl_samples)}"
    )

    print(
        f"{treatment_group} 样本数："
        f"{len(treat_samples)}"
    )


    # ========================================================
    # 3. 计算 Fold Change
    # ========================================================

    ctrl_mean = ctrl.mean(axis=1)

    treat_mean = treat.mean(axis=1)


    # 防止出现除0
    epsilon = 1e-12

    fold_change = (
        treat_mean + epsilon
    ) / (
        ctrl_mean + epsilon
    )


    log2fc = np.log2(
        fold_change
    )


    # ========================================================
    # 4. 独立样本 t 检验
    # ========================================================

    t_stat, p_values = ttest_ind(
        ctrl.values,
        treat.values,
        axis=1,
        equal_var=False,
        nan_policy="omit"
    )


    p_values = pd.Series(
        p_values,
        index=expr_df.index
    )


    t_stat = pd.Series(
        t_stat,
        index=expr_df.index
    )


    # ========================================================
    # 5. PLS 回归计算近似 VIP
    # ========================================================

    # X:
    # 行 = 样本
    # 列 = 代谢物
    X = expr_df.T.values


    # Control = 0
    # Treatment = 1
    y = (
        group_series == treatment_group
    ).astype(int).values


    pls = PLSRegression(
        n_components=1,
        scale=True
    )

    pls.fit(
        X,
        y
    )


    W = pls.x_weights_

    T = pls.x_scores_


    # --------------------------------------------------------
    # 单成分情况下的近似 VIP
    # --------------------------------------------------------

    S = np.sum(
        T ** 2
    )


    numerator = (
        expr_df.shape[0]
        * (W ** 2).flatten()
        * S
    )


    denominator = np.sum(
        (W ** 2) * S
    )


    vip = np.sqrt(
        numerator / denominator
    )


    vip = pd.Series(
        vip,
        index=expr_df.index
    )


    # ========================================================
    # 6. 汇总统计结果
    # ========================================================

    diff_stats = pd.DataFrame(
        {
            "Control_mean":
                ctrl_mean,

            "Treatment_mean":
                treat_mean,

            "FoldChange":
                fold_change,

            "log2FC":
                log2fc,

            "T_stat":
                t_stat,

            "P_value":
                p_values,

            "-log10P":
                -np.log10(
                    p_values.clip(
                        lower=1e-300
                    )
                ),

            "VIP":
                vip
        }
    )


    # ========================================================
    # 7. 添加 Up / Down / NS 分类
    # ========================================================

    def classify(row):

        if (
            row["VIP"] > vip_thresh
            and
            row["P_value"] < p_thresh
            and
            abs(row["log2FC"]) >= log2fc_thresh
        ):

            if row["log2FC"] > 0:
                return "Up"

            else:
                return "Down"

        return "NS"


    diff_stats["Regulation"] = (
        diff_stats.apply(
            classify,
            axis=1
        )
    )


    # ========================================================
    # 8. 筛选差异代谢物
    # ========================================================

    sig = diff_stats[
        diff_stats["Regulation"].isin(
            ["Up", "Down"]
        )
    ].copy()


    # --------------------------------------------------------
    # 按 |log2FC| 排序
    # --------------------------------------------------------

    sig["abs_log2FC"] = (
        sig["log2FC"].abs()
    )


    sig = sig.sort_values(
        [
            "abs_log2FC",
            "P_value"
        ],
        ascending=[
            False,
            True
        ]
    )


    # ========================================================
    # 9. 输出统计
    # ========================================================

    up_n = (
        diff_stats["Regulation"]
        == "Up"
    ).sum()

    down_n = (
        diff_stats["Regulation"]
        == "Down"
    ).sum()

    ns_n = (
        diff_stats["Regulation"]
        == "NS"
    ).sum()


    print("\n==============================")
    print("差异代谢物筛选结果")
    print("==============================")

    print(
        f"筛选标准："
        f"VIP > {vip_thresh}, "
        f"|log2FC| >= {log2fc_thresh}, "
        f"P < {p_thresh}"
    )

    print(
        f"Up: {up_n}"
    )

    print(
        f"Down: {down_n}"
    )

    print(
        f"NS: {ns_n}"
    )

    print(
        f"Total differential metabolites: "
        f"{sig.shape[0]}"
    )


    return diff_stats, sig


# ============================================================
# 3. 执行差异筛选
# ============================================================

diff_stats, sig = (
    differential_metabolite_screening(
        expr_df=expr_df,
        group_series=group_series,

        control_group="Control",
        treatment_group="Treatment",

        vip_thresh=1.0,
        log2fc_thresh=1.0,
        p_thresh=0.05
    )
)


# ============================================================
# 4. 查看前几行结果
# ============================================================

print("\n所有代谢物统计结果：")

print(
    diff_stats.head()
)


print("\n差异代谢物：")

print(
    sig.head(20)
)


# ============================================================
# 5. 保存结果
# ============================================================

diff_stats.to_csv(
    "All_Metabolite_Statistics.csv",
    encoding="utf-8-sig"
)


sig.to_csv(
    "Differential_Metabolites.csv",
    encoding="utf-8-sig"
)


expr_df.to_csv(
    "Expression_Matrix.csv",
    encoding="utf-8-sig"
)


group_series.to_csv(
    "Sample_Group.csv",
    encoding="utf-8-sig"
)


print(
    "\n✅ 已保存："
    "\nAll_Metabolite_Statistics.csv"
    "\nDifferential_Metabolites.csv"
    "\nExpression_Matrix.csv"
    "\nSample_Group.csv"
)
