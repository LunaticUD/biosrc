#!/bin/bash
# use_gene_catalog.sh
gene_catalog="gene_catalog_final.faa"
sample_genes="sample1_genes.faa"
# 将样本基因比对到基因目录
# 使用diamond进行快速比对
diamond makedb --in $gene_catalog -d gene_catalog_db
diamond blastp  -d gene_catalog_db \
				-q $sample_genes \
				-o sample1_vs_catalog.m8 \
				--max-target-seqs 1 \
				--evalue 1e-5 \
				--threads 8
# 统计基因覆盖度
echo "样本基因在目录中的覆盖度:"
total_genes=$(grep -c '^>' $sample_genes)
matched_genes=$(cut -f1 sample1_vs_catalog.m8 | sort -u | wc -l)
coverage=$(echo "scale=2; $matched_genes/$total_genes*100" | bc)
echo "  总基因数: $total_genes"
echo "  匹配基因: $matched_genes"
echo "  覆盖度: ${coverage}%"