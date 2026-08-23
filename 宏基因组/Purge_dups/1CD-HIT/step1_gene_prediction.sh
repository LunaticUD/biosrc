##!/bin/bash
#step1_gene_prediction.sh
# 使用Prodigal进行基因预测
input_contigs="contigs.fasta"
output_genes="predicted_genes.faa"
# 预测蛋白质编码基因
prodigal -i $input_contigs \
	 	-a $output_genes \
	 	-p meta \
	 	-q
echo "基因预测完成！共预测 $(grep -c '^>' $output_genes) 个基因"