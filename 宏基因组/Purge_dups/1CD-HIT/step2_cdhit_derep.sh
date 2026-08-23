##!/bin/bash
# step2_cdhit_derep.sh
input_genes="predicted_genes.faa"
output_nr="nonredundant_genes.faa"
clustering_info="cdhit_clusters.clstr"
# 设置相似度阈值（通常90-95%用于蛋白质）
similarity=0.95
echo "开始CD-HIT去冗余，相似度阈值: $similarity"
cd-hit -i $input_genes \
   	-o $output_nr \
   	-c $similarity \
   	-n 5 \
   	-M 16000 \
   	-T 8 \
   	-d 0 \
   	-sc 1
echo "去冗余完成！"
echo "输入基因数: $(grep -c '^>' $input_genes)"
echo "去冗余后基因数: $(grep -c '^>' $output_nr)"
echo "聚类信息保存在: $clustering_info"