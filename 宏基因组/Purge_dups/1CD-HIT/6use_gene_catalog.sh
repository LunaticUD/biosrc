#!/bin/bash
# use_gene_catalog.sh
# 将样本(全部预测)基因比对到基因目录，统计覆盖度
gene_catalog="gene_catalog_final.faa"
# 优先用命令行参数指定查询文件，默认用合并后的全部预测基因
sample_genes="${1:-all_samples_genes.faa}"
out_m8="all_samples_vs_catalog.m8"

# 输入检查
if [ ! -f "$gene_catalog" ]; then
    echo "错误: 基因目录 $gene_catalog 不存在，请先运行 5build_gene_catalog.sh" >&2
    exit 1
fi
if [ ! -f "$sample_genes" ]; then
    echo "错误: 查询基因文件 $sample_genes 不存在" >&2
    exit 1
fi

# 将样本基因比对到基因目录
diamond makedb --in $gene_catalog -d gene_catalog_db || exit 1
diamond blastp -d gene_catalog_db \
                -q $sample_genes \
                -o $out_m8 \
                --max-target-seqs 1 \
                --evalue 1e-5 \
                --threads 8 || exit 1

# 统计基因覆盖度
echo "样本基因在目录中的覆盖度:"
total_genes=$(grep -c '^>' $sample_genes)
matched_genes=$(cut -f1 $out_m8 | sort -u | wc -l)
if [ -z "$total_genes" ] || [ "$total_genes" -eq 0 ]; then
    echo "错误: 查询文件为空" >&2
    exit 1
fi
coverage=$(echo "scale=2; $matched_genes/$total_genes*100" | bc)
echo "  总基因数: $total_genes"
echo "  匹配基因: $matched_genes"
echo "  覆盖度: ${coverage}%"
