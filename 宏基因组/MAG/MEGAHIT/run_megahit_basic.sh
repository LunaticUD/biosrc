#!/bin/bash
# run_megahit_basic.sh
# 输入文件
forward="sample_clean_1.fastq.gz"
reverse="sample_clean_2.fastq.gz"
output_dir="megahit_results"
# 运行MEGAHIT
megahit \
	-1 $forward \
	-2 $reverse \
	-o $output_dir \
	--min-contig-len 500 \
	-t 8  
	# 使用8个线程
echo "MEGAHIT组装完成！"
echo "结果保存在 $output_dir 目录中"
echo "主要的contigs文件：$output_dir/final.contigs.fa"