#!/bin/bash
# build_gene_catalog.sh
# 1. 合并所有样本的基因
cat */predicted_genes.faa > all_samples_genes.faa
# 2. 去冗余（构建核心基因集）
cd-hit  -i all_samples_genes.faa \
		-o gene_catalog.faa \
		-c 0.95 \          # 95%相似度
		-n 5 \
		-M 64000 \         # 64GB内存
		-T 32 \
		-g 1 \
		-d 0
# 3. 为每个基因分配唯一ID
awk '/^>/ {printf(">GC_%08d\n", ++i); next} {print}' gene_catalog.faa > gene_catalog_final.faa
# 4. 创建基因信息表
echo "gene_id,original_id,length,cluster" > gene_catalog_info.csv
grep '^>' gene_catalog.faa | sed 's/^>//' | awk -F' ' '{print NR","$1","length($0)}' >> gene_catalog_info.csv
echo "基因目录构建完成！"
echo "总基因数: $(grep -c '^>' gene_catalog_final.faa)"
