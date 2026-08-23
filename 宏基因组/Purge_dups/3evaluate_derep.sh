#!/bin/bash
# evaluate_derep.sh
input_faa="predicted_genes.faa"
nr_faa="nonredundant_genes.faa"
echo "=== 去冗余效果评估 ==="
echo "1. 序列数量变化:"
input_count=$(grep -c '^>' $input_faa)
nr_count=$(grep -c '^>' $nr_faa)
reduction_rate=$(echo "scale=2; (1 - $nr_count/$input_count)*100" | bc)
echo "   输入序列: $input_count"
echo "   输出序列: $nr_count"
echo "   减少比例: ${reduction_rate}%"
echo "2. 序列长度统计:"
echo "   输入平均长度: $(seqkit stats $input_faa -T | tail -1 | awk '{print $6}')"
echo "   输出平均长度: $(seqkit stats $nr_faa -T | tail -1 | awk '{print $6}')"
echo "3. 聚类信息:"
if [ -f "cdhit_clusters.clstr" ]; then
	cluster_count=$(grep -c '^>' cdhit_clusters.clstr)
	echo "   聚类数量: $cluster_count"
	echo "   平均聚类大小: $(echo "scale=2; $input_count/$cluster_count" | bc)"
fi