# mamba install -c conda-forge -c bioconda \
#   r-biocmanager \
#   bioconductor-clusterprofiler \
#   bioconductor-org.mm.eg.db \
#   bioconductor-enrichplot \
#   bioconductor-dose \
#   r-data.table \
#   r-ggplot2
options(timeout = 600)  # 默认60秒下载超时对KEGG不够，调大到600秒
# ======加载======
library(data.table)# 快速读写表格文件
library(clusterProfiler)# 富集分析核心包，GO/KEGG/自定义基因集都靠它
library(org.Mm.eg.db)# 小鼠基因注释数据库
library(ggplot2)# 绘图
library(DOSE)# clusterProfiler底层依赖，提供部分可视化和统计支持
library(enrichplot)# 富集结果专用可视化（dotplot、gseaplot2等）

# ======读取DESeq2 分析生成的结果表======
deg<-fread("~/data/RNAseq/Exon_Mortix_DEG/DEG_treatment_vs_control.tsv")
# 设定 ORA 用的筛选阈值
p_threshold<-0.05
fc_threshold<-2 # 对应 log2FoldChange > log2(2) = 1

# ======GO 富集（GSEA 方法）======
# GSEA 需要一个按表达变化幅度排好序、并且带有基因名的向量作为输入：
# 按log2FoldChange从大到小排序
deg_order_fc<-deg[order(-log2FoldChange)]
logfc<-deg_order_fc$log2FoldChange
names(logfc)<-deg_order_fc$V1 # 基因名作为向量的名字
# 跑GSEA（小鼠用org.Mm.eg.db；人用org.Hs.eg.db；大鼠用org.Rn.eg.db）
enrich_go_gsea<-gseGO(geneList=logfc, OrgDb=org.Mm.eg.db, ont="ALL",
            pvalueCutoff=0.05, keyType="SYMBOL", verbose=FALSE)
# ======保存结果表======
enrich_go_gsea_df<-enrich_go_gsea@result
fwrite(enrich_go_gsea_df,"enrich_go_gsea_df.tsv", sep="\t")
# ont = "ALL"表示同时对 MF、BP、CC 三大类都做富集，如果只关心其中一类，可以改成"BP"、"MF"或"CC"。

# ======气泡图======
dotplot_enrich_go_gsea<-dotplot(enrich_go_gsea, showCategory=10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_gsea.png", dotplot_enrich_go_gsea,
   device="png", units="cm", width=16, height=18)
# ======单通路富集曲线======
# 从上面保存的enrich_go_gsea_df里挑一个你感兴趣的GO term编号
GOID<-"GO:0010634"
gsea_plot<-gseaplot2(enrich_go_gsea, geneSetID=GOID,
           title=subset(enrich_go_gsea_df, ID==GOID)$Description)
ggsave(paste0("gsea_plot_", gsub(":","", GOID),".png"), gsea_plot,
   device="png", units="cm", width=18, height=16)

# ======KEGG 富集（GSEA 方法）======
# KEGG 要求用 Entrez ID，不能直接用基因名
# 把基因名转成Entrez ID（KEGG分析的硬性要求）
logfc_sym<-logfc # 保存符号版logfc，MSigDB基因集为基因符号，GSEA要求两边ID一致
entrez_ids<-mapIds(
 x=org.Mm.eg.db,# 人换org.Hs.eg.db，大鼠换org.Rn.eg.db
 keys=names(logfc),
 keytype="SYMBOL",
 column="ENTREZID",
 multiVals="first"# 一个基因名对应多个Entrez ID时，取第一个
)
names(logfc)<-entrez_ids
logfc<-logfc[!is.na(names(logfc))]# 去掉转换失败（没匹配上）的基因
# KEGG物种代码："mmu"=小鼠, "hsa"=人, "rno"=大鼠
# ---- 离线替代 gseKEGG：rest.kegg.jp 不可达，改用本地 Enrichr 小鼠 KEGG 库 ----
kegg_raw<-fread("KEGG_2019_Mouse.txt", header=FALSE, sep="\t", fill=Inf)
kegg_t2g<-melt(kegg_raw, id.vars="V1", variable.factor=FALSE)[value!="" & !is.na(value), .(term=V1, gene=trimws(value))]
allsym<-keys(org.Mm.eg.db, keytype="SYMBOL"); names(allsym)<-toupper(allsym)
kegg_t2g[, pgs:=unname(allsym[toupper(gene)])]
kegg_t2g[, gene:=unname(mapIds(org.Mm.eg.db, keys=unique(pgs[!is.na(pgs)]), keytype="SYMBOL", column="ENTREZID", multiVals="first")[pgs])]
kegg_t2g[, pgs:=NULL]
kegg_t2g<-unique(kegg_t2g[!is.na(gene) & gene!=""])
enrich_kegg_gsea<-GSEA(logfc, TERM2GENE=kegg_t2g, pvalueCutoff=0.05)
enrich_kegg_gsea_df<-enrich_kegg_gsea@result
fwrite(enrich_kegg_gsea_df,"enrich_kegg_gsea_df.tsv", sep="\t")

# ======GO 富集（ORA 方法）======
# 筛出上调的显著差异基因
deg_up<-deg[log2FoldChange>log2(fc_threshold)&padj<p_threshold]$V1
# deg_dn <- deg[log2FoldChange < -log2(fc_threshold) & padj < p_threshold]$V1  # 下调基因
enrich_go_fet_up<-enrichGO(gene=deg_up, OrgDb=org.Mm.eg.db, keyType="SYMBOL",
              ont="ALL", pvalueCutoff=0.05,
              pAdjustMethod="BH", qvalueCutoff=0.05)
enrich_go_fet_up_df<-enrich_go_fet_up@result
fwrite(enrich_go_fet_up_df,"enrich_go_fet_up_df.tsv", sep="\t") 

# ======可视化======
dotplot_enrich_go_fet_up<-dotplot(enrich_go_fet_up, showCategory=10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_fet_up.png", dotplot_enrich_go_fet_up,
   device="png", units="cm", width=16, height=18)

# ======KEGG 富集（ORA 方法）======
enrich_kegg_fet_up<-enricher(entrez_ids[names(entrez_ids)%in%deg_up], TERM2GENE=kegg_t2g, pvalueCutoff=0.05, qvalueCutoff=0.05)
enrich_kegg_fet_up_df<-enrich_kegg_fet_up@result
fwrite(enrich_kegg_fet_up_df,"enrich_kegg_fet_up_df.tsv", sep="\t")

# ======MSigDB 富集分析（进阶，覆盖面最广）======
# GSEA 方法
msigdb<-read.gmt("./msigdb.v2024.1.Mm.symbols.gmt")# 换成你下载的实际文件名
enrich_msigdb_gsea<-GSEA(logfc_sym, TERM2GENE=msigdb) # logfc此时已是Entrez，改用符号版logfc_sym
enrich_msigdb_gsea_df<-enrich_msigdb_gsea@result
fwrite(enrich_msigdb_gsea_df,"enrich_msigdb_gsea_df.tsv", sep="\t")
# ORA 方法
msigdb<-read.gmt("./msigdb.v2024.1.Mm.symbols.gmt")
enrich_msigdb_fet_up<-enricher(deg_up, TERM2GENE=msigdb)
enrich_msigdb_fet_up_df<-enrich_msigdb_fet_up@result
fwrite(enrich_msigdb_fet_up_df,"enrich_msigdb_fet_up_df.tsv", sep="\t")
