#!/bin/bash
# run_megahit_advanced.sh
forward="sample_clean_1.fastq.gz"
reverse="sample_clean_2.fastq.gz"
megahit \
	-1 $forward \
	-2 $reverse \
	-o megahit_advanced \
	--k-list 21,29,39,59,79,99,119,141 \
	--min-count 2 \
	--prune-level 2 \
	--min-contig-len 1000 \
	--presets meta-sensitive \
	-t 16 \
	-m 0.9
	# 使用90%的可用内存
echo "高级模式MEGAHIT组装完成！"