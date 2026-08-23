#!/bin/bash
# build_gene_catalog.sh
# 1. 合并所有样本的基因
# 优先使用各样本子目录的 predicted_genes.faa；否则用顶层的单个 predicted_genes.faa
if [ -s all_samples_genes.faa ]; then
    echo "使用已有的 all_samples_genes.faa"
else
    shopt -s nullglob
    sample_files=( */predicted_genes.faa )
    if [ ${#sample_files[@]} -gt 0 ]; then
        cat "${sample_files[@]}" > all_samples_genes.faa
    elif [ -s predicted_genes.faa ]; then
        cp predicted_genes.faa all_samples_genes.faa
    else
        echo "错误: 未找到预测基因文件 (*/predicted_genes.faa 或 predicted_genes.faa)" >&2
        exit 1
    fi
fi
echo "输入基因数: $(grep -c '^>' all_samples_genes.faa)"

# 2. 去冗余（构建核心基因集）
# 95%相似度
cd-hit -i all_samples_genes.faa \
    -o gene_catalog.faa \
    -c 0.95 \
    -n 5 \
    -M 64000 \
    -T 32 \
    -g 1 \
    -d 0
if [ ! -s gene_catalog.faa ]; then
    echo "错误: cd-hit 执行失败，gene_catalog.faa 未生成" >&2
    exit 1
fi
echo "去冗余后基因数: $(grep -c '^>' gene_catalog.faa)"

# 3. 为每个基因分配唯一ID
awk '/^>/ {printf(">GC_%08d\n", ++i); next} {print}' gene_catalog.faa > gene_catalog_final.faa

# 4. 创建基因信息表
echo "gene_id,original_id,length,cluster" > gene_catalog_info.csv
grep '^>' gene_catalog.faa | sed 's/^>//' | awk -F' ' '{print NR","$1","length($0)}' >> gene_catalog_info.csv

echo "基因目录构建完成！"
echo "总基因数: $(grep -c '^>' gene_catalog_final.faa)"
