#!/bin/bash
# run_metaspades_basic.sh
forward="sample_clean_1.fastq.gz"
reverse="sample_clean_2.fastq.gz"
spades.py \
	--meta \
	-1 $forward \
	-2 $reverse \
	-o metaspades_results \
	-t 20 \
	-m 64  # 限制使用64GB内存
echo "metaSPAdes组装完成！"
echo "结果保存在 metaspades_results 目录中"
echo "主要的contigs文件：metaspades_results/contigs.fasta"