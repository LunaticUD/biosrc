#============安装============
# # 安装 BiocManager（用于管理 Bioconductor 上的包）
# if(!require("BiocManager", quietly =TRUE))
#   install.packages("BiocManager")
# # 差异表达分析核心包（来自 Bioconductor）
# BiocManager::install("DESeq2")
# # 数据处理、绘图、质控相关的包（来自 CRAN）
# install.packages(c("data.table","ggplot2","FactoMineR","factoextra","pheatmap","psych"))
# bash
# mamba activate biosrc

# mamba install -c conda-forge -c bioconda \
#     bioconductor-deseq2 \
#     r-data.table \
#     r-ggplot2 \
#     r-factominer \
#     r-factoextra \
#     r-pheatmap \
#     r-psych

# 加载
library(data.table)# 快速读写表格文件
library(DESeq2)# 差异表达分析核心包，同时提供vst标准化
library(ggplot2)# 绘图
library(FactoMineR)# 计算PCA
library(factoextra)# 可视化PCA结果
library(pheatmap)# 画热图（相关性图、聚类图）
library(psych)# 计算样本间相关系数
#============设定工作目录============
setwd("~/data/RNAseq")
# #============计数转矩阵============
# 找到所有样本的 featureCounts 结果文件
count_file_path <-"~/data/RNAseq/GSE259357/aligned"
file_names <- list.files(count_file_path, pattern ="featureCounts_exon\\.txt$",
                         recursive =TRUE, full.names =TRUE)
# 逐个读入每个样本的count文件（featureCounts输出默认第1列是基因名，第7列是count数，跳过前1行注释）
count_list <- lapply(file_names,function(f){
	tmp <- fread(f, skip =1, data.table =FALSE)
	sample_name <- gsub("_featureCounts_exon.*","", basename(f))
	data.frame(Gene = tmp[[1]], Count = tmp[[7]])|> setNames(c("Gene", sample_name))
})
# 以基因名为key，把所有样本的count合并成一张矩阵
counts_raw <- Reduce(function(x, y) merge(x, y, by ="Gene"), count_list)
rownames(counts_raw)<- counts_raw$Gene
counts_raw <- counts_raw[,-1]
# 保存合并好的原始表达矩阵，以后直接读这个文件就行，不用每次重新合并
fwrite(counts_raw,"counts_raw.tsv", sep ="\t", row.names =TRUE)
# #============读取表达矩阵============
count_tbl <- fread("counts_raw.tsv", data.table =FALSE)
# #============把第一列基因名设为行名============
rownames(count_tbl)<- count_tbl[[1]]
count_tbl <- count_tbl[,-1]
dim(count_tbl)# 看看现在是多少个基因 x 多少个样本
colnames(count_tbl)# 确认样本名都对
# #============过滤低表达/不表达的基因============
perc_keep <- 0.8   # 要求至少在80%的样本里检测到表达
gene_keep <- rowSums(count_tbl >0)>=ceiling(perc_keep * ncol(count_tbl))
count_tbl_filtered <- count_tbl[gene_keep,]
dim(count_tbl_filtered)# 过滤后剩下的基因数
# #============创建样本分组信息（Metadata）============
# meta <- data.frame(
#   SampleID = colnames(count_tbl_filtered),
#   Treatment = ifelse(grepl("^control", colnames(count_tbl_filtered)),
# "control","treatment")
# )
# BUG: 定义实验与对照组
meta <- data.frame(
  SampleID = colnames(count_tbl_filtered),
  Treatment = ifelse(colnames(count_tbl_filtered) %in% c("SRR28119110","SRR28119111"),
                     "control", "treatment")
)

rownames(meta)<- meta$SampleID
meta

colnames(count_tbl_filtered)
table(meta$Treatment)
# #============构建 DESeqDataSet 并做 vst 标准化============
dds <- DESeqDataSetFromMatrix(countData =round(count_tbl_filtered),
                              colData = meta,
                              design =~ Treatment)

# vst：方差稳定转换，同时完成文库大小校正
vsd <- vst(dds, blind =TRUE)
# 保存一下，后面差异分析直接读这个对象就行，不用重新跑一遍
saveRDS(dds,"dds.rds")
# 取出标准化后的表达矩阵，后面三张图都用它
exp_vst <- assay(vsd)
# # ============图一：PCA 图============
pca_data <- PCA(t(exp_vst), scale.unit =TRUE, ncp =5, graph =FALSE)
pca_plot <- fviz_pca_ind(pca_data,
                         geom.ind =c("point","text"),
                         pointsize =4,
                         labelsize =3,
                         repel =TRUE,
                         col.ind = meta$Treatment,
                         addEllipses =TRUE,
                         mean.point =FALSE,
                         palette =c("#0074b3","#982b2b"))+
  ggtitle("PCA: Control vs Treatment")+
  theme_bw()
pca_plot
pdf("PCA_Plot.pdf", width =8, height =6)
pca_plot
dev.off()
pca_data$eig   # 查看每个主成分能解释多少比例的变异，PC1+PC2一般看个大概就行
# #============图二：样本相关性热图============
# 挑出变异程度最高的500个基因（用MAD衡量变异程度，比标准差更抗离群值干扰）
top_var_genes <-names(sort(apply(exp_vst,1, mad), decreasing =TRUE)[1:500])
# 计算样本两两之间的Spearman相关系数
exp_cor <- psych::corr.test(exp_vst[top_var_genes,], method ="spearman")
cor_matrix <- exp_cor$r
annotation_col <- data.frame(Treatment = meta$Treatment)
rownames(annotation_col)<- meta$SampleID
correlation_plot <- pheatmap(cor_matrix,
                             scale ="none",
                             fontsize_row =8,
                             fontsize_col =8,
                             angle_col =45,
                             cluster_rows =TRUE,
                             cluster_cols =TRUE,
                             annotation_col = annotation_col,
                             color = colorRampPalette(c("#0f86a9","white","#ed8b10"))(100),
                             main ="Sample Correlation Heatmap")
pdf("correlation_heatmap.pdf", width =10, height =8)
correlation_plot
dev.off()
# #============图三：层次聚类热图============
exp_cluster <- exp_vst[top_var_genes,]
ann_colors <-list(Treatment =c(control ="#0074b3", treatment ="#982b2b"))
cluster_plot <- pheatmap(exp_cluster,
                         scale ="row",
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         fontsize_col =8,
                         angle_col =45,
                         cluster_rows =TRUE,
                         cluster_cols =TRUE,
                         color = colorRampPalette(c("#0f86a9","white","#ed8b10"))(100),
                         show_rownames =FALSE,
                         main ="Top Variable Genes Clustering")
pdf("clustering_heatmap.pdf", width =12, height =8)
cluster_plot
dev.off()
# #============差异表达分析============
# 用真实的design重新跑模型（前面vst用的是blind=TRUE，这里要用真实分组信息估计离散度）
dds <- DESeq(dds)
# 提取"处理组 vs 对照组"的检验结果
res <- results(dds, contrast =c("Treatment","treatment","control"))
# 按adjusted p value从小到大排序
res <- res[order(res$padj),]
deg_tbl <- as.data.frame(res)
fwrite(deg_tbl,"DEG_treatment_vs_control.tsv", sep ="\t", row.names =TRUE)
# #============输出差异基因============
deg_sig <- subset(deg_tbl, padj <0.05&abs(log2FoldChange)>1)
dim(deg_sig)# 看看筛出来多少个显著差异基因
