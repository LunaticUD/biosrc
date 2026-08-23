#!/bin/bash
# run_metaspades_advanced.sh
forward="sample_clean_1.fastq.gz"
reverse="sample_clean_2.fastq.gz"
spades.py \
	--meta \
	-1 $forward \
	-2 $reverse \
	-o metaspades_advanced \
	-k 21,33,55,77,99,127 \
	--only-assembler \
	--careful \
	-t 32 \
	-m 256
echo "高级metaSPAdes组装完成！"