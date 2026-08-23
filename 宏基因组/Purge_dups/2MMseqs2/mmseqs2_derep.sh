#!/bin/bash
# mmseqs2_derep.sh
input_genes="predicted_genes.faa"
output_prefix="mmseqs2_nr"
tmp_dir="mmseqs2_tmp"
# 创建数据库
mmseqs createdb $input_genes ${output_prefix}DB
# 聚类（90%相似度，80%覆盖度）
mmseqs cluster ${output_prefix}DB ${output_prefix}Cluster $tmp_dir \
	--min-seq-id 0.9 \
	--cov-mode 1 \
	-c 0.8 \
	--threads 8
# 提取代表性序列
mmseqs createtsv ${output_prefix}DB ${output_prefix}DB ${output_prefix}Cluster ${output_prefix}_cluster.tsv
mmseqs result2repseq ${output_prefix}DB ${output_prefix}Cluster ${output_prefix}DB_rep
mmseqs result2flat ${output_prefix}DB ${output_prefix}DB ${output_prefix}DB_rep ${output_prefix}.faa
echo "MMseqs2去冗余完成！"
echo "输出文件: ${output_prefix}.faa"
echo "聚类关系: ${output_prefix}_cluster.tsv"