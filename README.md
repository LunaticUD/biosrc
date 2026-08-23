# Biosrc
## 宏基因组
1.  Lefse
2.  差异OTU和差异物种曼哈顿图
3.  物种Alpha多样性箱线图
4.  物种Beta多样性NMDS分析
5.  物种Beta多样性PCoA分析
6.  物种三元图
7.  物种堆积图
8.  物种组成Circos弦图
9.  物种组成分组聚类热图
10. 物种组成气泡图
11. 物种韦恩图
12. [宏基因组去冗余](./宏基因组/宏基因组组装流程.md)
	1) FastQC --- 质控
		fastp、Fastqc+trimmomatic、kneaddata
	2) MAG --- 宏基因组组装
		MEGAHIT、metaSPAdes
	3) Purge_dups --- 宏基因组去冗余
		CD-HIT、MMseqs2、UCLUST、VSEARCH、dRep、评估、可视化、基因目录
## [代谢组](./代谢组/代谢组流程.md)
### Step 1：差异代谢物筛选（VIP + FC + P）
1. 计算 Fold Change & P 值
2. 计算 VIP（OPLS-DA）
3. 近似 VIP（科研中常用简化版）
### Step 2：差异代谢物模式分析
#### 差异代谢物筛选与分析
1. 火山图
2. 差异代谢物条形图
#### 差异代谢物模式分析
1. 相关性热图
2. 相关性网路图
3. 双向聚类热图
4. K-Means 趋势分析
#### 功能注释及富集分析
1. 差异代谢物通路分类图
2. 差异代谢物通路富集图
## Step 3：KEGG富集通路
1. KEGG富集分析柱图
2. KEGG富集分析气泡图
## [转录组](./转录组/转录组分析流程.md)
### Step 1：Fastq到计数表
### Step 2：从计数表到表达矩阵到差异基因
### Step 3：从差异基因到富集通路
