#!/bin/bash
# step3_analyze_clusters.sh
# 统计聚类大小分布
clstr_file="${1:-nonredundant_genes.faa.clstr}"
fasta_file="${2:-nonredundant_genes.faa}"

if [ ! -f "$clstr_file" ]; then
    echo "错误: 聚类文件 $clstr_file 不存在" >&2
    echo "用法: bash step3_analyze_clusters.sh [聚类.clstr] [代表序列.faa]" >&2
    exit 1
fi

echo "聚类大小统计（$clstr_file）："
# 每个 >Cluster N 块内的成员行数即聚类大小
awk '/^>/ {if(size>0) print size; size=0} /^[0-9]/ {size++} END {if(size>0) print size}' "$clstr_file" \
    | sort -n | uniq -c | sort -k2n > cluster_size_distribution.txt

echo "聚类大小分布已保存到 cluster_size_distribution.txt"

# 提取代表性序列ID
if [ -f "$fasta_file" ]; then
    grep '^>' "$fasta_file" | sed 's/^>//' > representative_ids.txt
    echo "代表性序列ID已保存到 representative_ids.txt"
else
    echo "警告: $fasta_file 不存在，跳过代表性序列ID提取" >&2
fi
