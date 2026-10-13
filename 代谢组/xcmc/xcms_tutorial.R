# ============================================================
# XCMS 完整分析流程教程（LC-MS 非靶向代谢组）
# 示例数据：faahKO（12 个样本：6 KO + 6 WT，正离子 QTOF）
#
# 目录
#   1  介绍
#   2  加载R包
#   3  数据下载
#   4  数据导入
#   5  获取mzML的样本信息
#   6  设置全局颜色
#   7  读取和评估原始数据（读取数据 / BPC / TIC / 查看单个特征）
#   8  特征检测（执行特征检测 / 总结统计信息）
#   9  保留时间校正（基于 OBIWarp / 基于 Peak Groups）
#   10 特征分组/分箱（测试区域评估 / 数据集特征分组）
#   11 填充缺失的峰值
#   12 特征总结（绘制特征直方图 / 保存特征数据）
#   13 统计评估特征（评估数据 / PCA / 层次聚类 / PAM / 输出表格）
#   17 差异分析（监督方法：t 检验 + BH-FDR + 火山图）—— 新增，出结论的关键一步
#   14 总结
#   15 系统信息
#   16 参考
# ============================================================

# ------------------------------------------------------------
## 1 介绍
# ------------------------------------------------------------
# 目标：用 XCMS 走完 LC-MS 非靶向代谢组"原始数据 -> 特征矩阵 -> 统计分析"全流程。
# 核心思想：LC-MS 是 m/z × RT × 强度 的三维信号；同一化合物在不同样品中会漂移，
#          故需"峰检测 -> RT校正 -> 分组 -> 补峰"把它归并成可比较的特征表。
# 数据：faahKO（Bioconductor 示例数据，无需联网）。
# 安装依赖：options(BioC_mirror = NULL); BiocManager::install(c("xcms","MSnbase","faahKO"))

# ------------------------------------------------------------
## 2 加载R包
# ------------------------------------------------------------
need <- c("xcms", "MSnbase", "faahKO")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("缺少包: ", paste(miss, collapse = ", "),
                       "\n请运行：options(BioC_mirror = NULL); ",
                       "BiocManager::install(c('xcms','MSnbase','faahKO'))")

suppressPackageStartupMessages({
  library(xcms)
  library(MSnbase)
  library(faahKO)
})

# 并行后端（可选）：
# Windows 上默认多进程(SnowParam)会打印无害警告 "载入时 'package:stats' 可能无用"，
# 那是子进程序列化时重复加载 stats 所致，不影响结果。
# 想要干净输出 -> 串行；想加速 -> 改用 register(SnowParam(4)) 并忽略该警告。
suppressPackageStartupMessages(library(BiocParallel))
register(SerialParam())

outdir <- "xcms_full_output"
dir.create(outdir, showWarnings = FALSE)

# ------------------------------------------------------------
## 3 数据下载
# ------------------------------------------------------------
# 方式A（本教程默认，离线）：faahKO 包自带 CDF 文件，无需下载。
# 方式B：从网络下载 mzML（示例，按需替换 URL）：
#   dir.create("mzML", showWarnings = FALSE)
#   download.file("https://example.org/sample01.mzML",
#                 "mzML/sample01.mzML", mode = "wb")
# 提示：厂商原始格式(.raw/.wiff)必须先用 ProteoWizard 的 msConvert 转成 mzML：
#   msconvert *.raw --mzML --filter "peakPicking vendor msLevel=1-"

# ------------------------------------------------------------
## 4 数据导入
# ------------------------------------------------------------
cdfs <- sort(dir(system.file("cdf", package = "faahKO"),
                 full.names = TRUE, recursive = TRUE))
raw <- readMSData(cdfs, mode = "onDisk")   # onDisk：大数据不全读进内存

# ------------------------------------------------------------
## 5 获取mzML的样本信息
# ------------------------------------------------------------
raw$sample_name  <- sub("\\.CDF$", "", basename(cdfs), ignore.case = TRUE)
raw$sample_group <- ifelse(grepl("^ko", basename(cdfs), ignore.case = TRUE), "KO", "WT")
try(polarity(raw) <- "positive", silent = TRUE)   # netCDF 不存极性，手动设

sample_info <- data.frame(file  = basename(cdfs),
                          name  = raw$sample_name,
                          group = raw$sample_group)
cat("== 样本信息 ==\n"); print(sample_info)

# ------------------------------------------------------------
## 6 设置全局颜色
# ------------------------------------------------------------
cols <- c(KO = "#185FA5", WT = "#993C1D")
gcol <- cols[raw$sample_group]

# ------------------------------------------------------------
## 7 读取和评估原始数据
# ------------------------------------------------------------
## 7.1 读取数据（概览）
print(raw)

## 7.2 绘制基峰色谱图 BPC（每个 RT 取最大强度）
png(file.path(outdir, "01_BPC.png"), 1100, 700)
plot(chromatogram(raw, aggregationFun = "max"), col = gcol, main = "BPC")
dev.off()

## 7.3 绘制总离子流图 TIC（每个 RT 取强度之和）
png(file.path(outdir, "02_TIC.png"), 1100, 700)
plot(chromatogram(raw, aggregationFun = "sum"), col = gcol, main = "TIC")
dev.off()

## 7.4 查看单个特征（提取某 m/z + RT 范围的 EIC）
mzr <- c(334.9, 335.1)      # 质量窗口
rtr <- c(2600, 3000)        # 保留时间窗口(秒)
png(file.path(outdir, "03_EIC.png"), 1100, 700)
plot(chromatogram(raw, mz = mzr, rt = rtr), col = gcol, main = "EIC  m/z 335")
dev.off()

# ------------------------------------------------------------
## 8 特征检测
# ------------------------------------------------------------
## 8.1 执行特征检测（centWave：高分辨数据）
cwp <- CentWaveParam(peakwidth = c(20, 80), snthresh = 10,
                     prefilter = c(3, 10000), noise = 5000, mzdiff = 0.01)
x <- findChromPeaks(raw, param = cwp)

## 8.2 总结统计信息
cat("== 峰总数:", nrow(chromPeaks(x)), "==\n")
cat("== 各样本峰数 ==\n"); print(table(chromPeaks(x)[, "sample"]))
png(file.path(outdir, "04_peaks_sample1.png"), 1100, 700)
plotChromPeaks(x, file = 1); dev.off()

# ------------------------------------------------------------
## 9 保留时间校正
# ------------------------------------------------------------
## 9.1 基于 OBIWarp 算法（不依赖峰组，逐样品非线性对齐）
x_obi <- adjustRtime(x, param = ObiwarpParam(binSize = 0.6))
png(file.path(outdir, "05_rt_obiwarp.png"), 1100, 700)
plotAdjustedRtime(x_obi, col = gcol); dev.off()

## 9.2 基于 Peak Groups 方法（关键：必须先分组！）
# PeakGroups 校正以"峰组(即特征)"为对齐锚点，因此对象里必须先有 featureDefinitions。
# 若直接 adjustRtime 会报：No feature definitions found in 'object'!
x_pg <- groupChromPeaks(x, param = PeakDensityParam(
          sampleGroups = raw$sample_group, bw = 30,
          minFraction = 0.5, minSamples = 1, binSize = 0.025))
x_pg <- adjustRtime(x_pg, param = PeakGroupsParam(minFraction = 0.9, extraPeaks = 1))
png(file.path(outdir, "06_rt_peakgroups.png"), 1100, 700)
plotAdjustedRtime(x_pg, col = gcol); dev.off()

# 两种方法二选一，后续沿用 OBIWarp 的结果
x <- x_obi

# ------------------------------------------------------------
## 10 特征分组/分箱
# ------------------------------------------------------------
pdp <- PeakDensityParam(sampleGroups = raw$sample_group,
                        bw = 30, minFraction = 0.5, minSamples = 1,
                        binSize = 0.025)

## 10.1 测试区域评估（先在某个 m/z 区域看分组是否合理）
# 注意：直接在 XCMSnExp 上调用 plotChromPeakDensity 已被弃用(discouraged)，
# 应先用 chromatogram() 抽出 XChromatograms 对象再画。
mz_test  <- 508.2
chr_test <- chromatogram(x, mz = mz_test + c(-0.01, 0.01))
png(file.path(outdir, "07_density_test_region.png"), 1100, 700)
plotChromPeakDensity(chr_test, param = pdp, col = gcol); dev.off()

## 10.2 数据集特征分组
x <- groupChromPeaks(x, param = pdp)
cat("== 特征数:", nrow(featureDefinitions(x)), "==\n")

# ------------------------------------------------------------
## 11 填充缺失的峰值
# ------------------------------------------------------------
# 【关键】补峰之前先留一份"实测值"底档。
# fillChromPeaks() 会按特征所在的 m/z-RT 区域把缺失格子补上数值，补完再看数据就是
# "全都没缺失"，无法区分哪些是真测到的、哪些是算法算出来的。
# 所以这里在补峰【前】取一次 featureValues，得到实测掩膜 obs_mask。
fv_obs   <- featureValues(x, method = "maxint", value = "into", missing = NA)
obs_mask <- !is.na(fv_obs)          # TRUE = 该特征在该样本有实测峰

x <- fillChromPeaks(x, param = ChromPeakAreaParam())

# ------------------------------------------------------------
## 12 特征总结
# ------------------------------------------------------------
fv <- featureValues(x, method = "maxint", value = "into", missing = NA)
fd <- as.data.frame(featureDefinitions(x))

## 12.1 绘制特征直方图（m/z 与 RT 分布）
png(file.path(outdir, "08_feature_hist.png"), 1100, 700)
par(mfrow = c(1, 2))
hist(fd$mzmed, breaks = 50, main = "特征 m/z 分布", xlab = "m/z")
hist(fd$rtmed, breaks = 50, main = "特征 RT 分布",  xlab = "RT (s)")
dev.off()

## 12.2 保存特征数据（peakidx 是列表列，先展平再写 CSV）
# 注意：feature_intensity.csv 是【补峰后】的矩阵，含填充值。
#       配套的 feature_observed_mask.csv 标出哪些格子是实测的，两者必须一起用。
write.csv(fv,       file.path(outdir, "feature_intensity.csv"))
write.csv(obs_mask, file.path(outdir, "feature_observed_mask.csv"))
fd_out <- fd
if ("peakidx" %in% names(fd_out))
  fd_out$peakidx <- vapply(fd_out$peakidx, paste, character(1), collapse = ";")
write.csv(fd_out, file.path(outdir, "feature_definition.csv"), row.names = FALSE)

## 12.3 填充率诊断：最终矩阵里有多少格子是补出来的？
fill_rate <- 1 - rowMeans(obs_mask)
cat("\n== 填充率诊断 ==\n")
cat(sprintf("总格子 %d：实测 %d (%.1f%%)，填充 %d (%.1f%%)\n",
            length(obs_mask), sum(obs_mask), 100 * mean(obs_mask),
            sum(!obs_mask), 100 * mean(!obs_mask)))
cat("每个特征被填充的格子数分布：\n"); print(table(round(fill_rate * ncol(fv))))
cat("12 个样本全部实测（零填充）的特征数：",
    sum(rowSums(obs_mask) == ncol(obs_mask)), "\n")
# 经验法则：填充率 > 30% 的特征不要拿去做组间比较——
# 那等于在比"算法怎么补的"，而不是比生物学。
cat("填充率 > 30% 的特征数：", sum(fill_rate > 0.3), "\n")

# ------------------------------------------------------------
## 13 统计评估特征
# ------------------------------------------------------------
## 13.1 评估数据：过滤 + log 变换 + 样本间中位数归一化
# 【重要】不要用 rowSums(!is.na(fv)) 过滤——fv 是补峰后的，几乎每行都"无缺失"，
# 那样等于没过滤（本例 294 个特征里会放行 271 个）。必须用补峰前的实测掩膜 obs_mask。
min_obs <- 9                                   # 要求至少 9/12 个样本有实测峰
keep <- rowSums(obs_mask) >= min_obs & rowSums(is.na(fv)) == 0
if (sum(keep) < 2) stop("保留特征过少（", sum(keep), "），请调低 min_obs")
cat(sprintf("\n检出率门槛 >= %d/%d → 保留 %d / %d 个特征\n",
            min_obs, ncol(fv), sum(keep), nrow(fv)))

mat <- log2(fv[keep, , drop = FALSE] + 1)
mat <- sweep(mat, 2, apply(mat, 2, median, na.rm = TRUE), "-")
# 说明：mat 仍含少量填充值（每个特征 <= 12 - min_obs 个），用于 PCA/聚类可视化尚可；
#       要做假设检验请用 §13.6 导出的纯实测矩阵，别用这一张。

## 13.2 主成分分析 PCA
pca <- prcomp(t(mat), scale. = TRUE)
grp <- factor(raw$sample_group, levels = c("KO", "WT"))
png(file.path(outdir, "09_PCA.png"), 900, 700)
plot(pca$x[, 1], pca$x[, 2], col = cols[grp], pch = 19, cex = 1.8,
     xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]),
     main = "PCA")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.8)
legend("topright", legend = levels(grp), col = cols, pch = 19)
dev.off()

## 13.3 层次聚类分析
hc <- hclust(dist(t(mat)), method = "ward.D2")
png(file.path(outdir, "10_hclust.png"), 900, 700)
plot(hc, main = "Hierarchical clustering (ward.D2)", xlab = "", sub = "")
dev.off()

## 13.4 PAM 聚类
suppressPackageStartupMessages(library(cluster))
pamres <- pam(t(mat), k = 2)

# 注意：plot(pamres) 默认用 princomp 画二维图，而 princomp 要求"样本数 > 变量数"；
# 组学数据是 n(样本) << p(特征)，会报 "'princomp'只能在单位比变量多的情况下使用"。
# 解决：改用基于 SVD 的 prcomp 得分来展示聚类（prcomp 无此限制）。
png(file.path(outdir, "11_pam.png"), 900, 700)
pcol <- c("#185FA5", "#993C1D")[pamres$clustering]
plot(pca$x[, 1], pca$x[, 2], col = pcol, pch = 19, cex = 1.8,
     xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]),
     main = "PAM clusters (k=2) on PCA")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.8)
legend("topright", legend = paste("Cluster", seq_len(2)),
       col = c("#185FA5", "#993C1D"), pch = 19)
dev.off()

# 轮廓系数图（不依赖 PCA，安全；用于评估聚类质量）
png(file.path(outdir, "11b_pam_silhouette.png"), 900, 700)
plot(silhouette(pamres), main = "PAM silhouette (k=2)")
dev.off()

## 13.5 输出表格
write.csv(as.data.frame(pca$x), file.path(outdir, "PCA_scores.csv"))
write.csv(mat, file.path(outdir, "feature_matrix_normalized.csv"))
cat("== 聚类结果(PAM) ==\n"); print(table(pamres$clustering, raw$sample_group))

## 13.6 导出"仅实测值"矩阵（下游差异分析请用这一张）
mat_obs <- mat
mat_obs[!obs_mask[keep, , drop = FALSE]] <- NA
write.csv(mat_obs, file.path(outdir, "feature_matrix_observed.csv"))
cat(sprintf("实测矩阵: %d 特征 x %d 样本；实测 %d 格 (%.1f%%)，缺失 %d 格\n",
            nrow(mat_obs), ncol(mat_obs),
            sum(!is.na(mat_obs)), 100 * mean(!is.na(mat_obs)),
            sum(is.na(mat_obs))))

# ------------------------------------------------------------
## 17 差异分析（监督方法：t 检验 + BH-FDR + 火山图）
# ------------------------------------------------------------
# 为什么必须有这一节：
# PCA 和聚类是【无监督】的——它们只最大化方差，并不知道你关心哪两组，
# 所以"分不开"不代表没有差异，只代表差异不是方差最大的那个方向。
# 只有把分组标签喂给算法（有监督），才能直接回答"KO 和 WT 到底有没有差异"。
# 同时必须做多重检验校正：109 个特征各检验一次，α=0.05 时纯靠运气
# 也会"发现" 109 × 0.05 ≈ 5.5 个假差异。

## 17.1 逐特征检验（用 §13.6 的【实测值】矩阵，填充值不参与统计）
grp_v <- as.character(raw$sample_group)     # 依赖"矩阵列顺序 == 样本顺序"
stopifnot(length(grp_v) == ncol(mat_obs))

diff_res <- data.frame(
  feature  = rownames(mat_obs),
  log2FC   = rowMeans(mat_obs[, grp_v == "KO", drop = FALSE], na.rm = TRUE) -
             rowMeans(mat_obs[, grp_v == "WT", drop = FALSE], na.rm = TRUE),
  p_ttest  = NA_real_,
  p_wilcox = NA_real_)

for (i in seq_len(nrow(mat_obs))) {
  a <- mat_obs[i, grp_v == "KO"]; a <- a[!is.na(a)]
  b <- mat_obs[i, grp_v == "WT"]; b <- b[!is.na(b)]
  if (length(a) >= 3 && length(b) >= 3) {           # 每组至少 3 个实测值才检验
    diff_res$p_ttest[i]  <- tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_)
    diff_res$p_wilcox[i] <- tryCatch(suppressWarnings(wilcox.test(a, b)$p.value),
                                     error = function(e) NA_real_)
  }
}

## 17.2 BH-FDR 校正（这一步不能省）
diff_res$p_ttest_BH  <- p.adjust(diff_res$p_ttest,  method = "BH")
diff_res$p_wilcox_BH <- p.adjust(diff_res$p_wilcox, method = "BH")
diff_res$significant <- !is.na(diff_res$p_ttest_BH) &
                        diff_res$p_ttest_BH < 0.05 & abs(diff_res$log2FC) > 1

## 17.3 结果汇总
cat("\n== 差异分析 ==\n")
cat("可检验特征:", sum(!is.na(diff_res$p_ttest)), "/", nrow(diff_res), "\n")
cat(sprintf("t 检验 未校正 p<0.05: %d  (其中约 %.1f 个可能是假阳性)\n",
            sum(diff_res$p_ttest < 0.05, na.rm = TRUE),
            sum(!is.na(diff_res$p_ttest)) * 0.05))
cat("t 检验 BH 校正后 p<0.05 :", sum(diff_res$p_ttest_BH  < 0.05, na.rm = TRUE), "\n")
cat("Wilcoxon BH 校正后 p<0.05:", sum(diff_res$p_wilcox_BH < 0.05, na.rm = TRUE), "\n")
cat("最终显著 (BH<0.05 且 |log2FC|>1):", sum(diff_res$significant), "\n")
print(head(diff_res[order(diff_res$p_ttest_BH),
                    c("feature", "log2FC", "p_ttest", "p_ttest_BH", "significant")], 10),
      row.names = FALSE)

# 为什么 Wilcoxon 在这个样本量下一个都检不出（这不是 bug）：
# n = 6 vs 6 时，Wilcoxon 双尾最小可能 p = 2 / C(12,6) = 2/924 ≈ 0.0022；
# 经 BH 校正后最好情况也只有 0.0022 × 109 ≈ 0.24 > 0.05 —— 数学上不可能显著。
# 结论：小样本代谢组学一般用 t 检验，或 limma 的 moderated t（借全体特征的信息稳定方差）。

## 17.4 火山图
png(file.path(outdir, "12_volcano.png"), 900, 700)
plot(diff_res$log2FC, -log10(diff_res$p_ttest),
     pch = 19, cex = 1.4,
     col = ifelse(diff_res$significant, "#993C1D", "#888780"),
     xlab = expression(log[2] ~ "FC  (KO / WT)"),
     ylab = expression(-log[10] ~ "(p, t 检验)"),
     main = "Volcano plot")
abline(h = -log10(0.05), lty = 2, col = "#B4B2A9")
abline(v = c(-1, 1),     lty = 2, col = "#B4B2A9")
sig <- which(diff_res$significant)
if (length(sig))
  text(diff_res$log2FC[sig], -log10(diff_res$p_ttest[sig]),
       labels = diff_res$feature[sig], pos = 3, cex = 0.85)
legend("topleft", bty = "n", pch = 19, col = c("#993C1D", "#888780"),
       legend = c("显著 (BH<0.05 且 |log2FC|>1)", "不显著"))
dev.off()

## 17.5 导出差异分析结果
write.csv(diff_res, file.path(outdir, "differential_results.csv"), row.names = FALSE)

# ------------------------------------------------------------
## 14 总结
# ------------------------------------------------------------
cat("== 流程完成 ==\n")
cat("样本数:", ncol(fv), " 特征数:", nrow(fv), "\n")
cat("检出率 >=", min_obs, "/", ncol(fv), " 参与统计:", sum(keep), "\n")
cat(sprintf("全部特征填充率 %.1f%%；其中填充率 >30%% 的 %d 个（不建议做组间比较）\n",
            100 * mean(!obs_mask), sum(fill_rate > 0.3)))
cat("下游差异分析请用 feature_matrix_observed.csv（纯实测值）\n")
cat(sprintf("差异分析: 显著特征 %d 个（t 检验 BH<0.05 且 |log2FC|>1）\n",
            sum(diff_res$significant)))
cat("输出目录:", normalizePath(outdir), "\n")
print(list.files(outdir))

# ------------------------------------------------------------
## 15 系统信息
# ------------------------------------------------------------
print(sessionInfo())

# ------------------------------------------------------------
## 16 参考
# ------------------------------------------------------------
# - xcms 官方文档:   https://bioconductor.org/packages/xcms
# - xcms 教程(vignette): 在 R 中运行  vignette("xcms")
# - MSnbase 文档:    https://bioconductor.org/packages/MSnbase
# - MetaboAnalyst（在线统计/通路）: https://www.metaboanalyst.ca
# - ProteoWizard(msConvert): https://proteowizard.sourceforge.io
