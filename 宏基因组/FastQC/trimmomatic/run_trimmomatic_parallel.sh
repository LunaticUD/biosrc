#!/bin/bash
# 先批量运行FastQC
echo "批量运行FastQC..."
fastqc *_1.fastq.gz *_2.fastq.gz -t 8
# 然后批量运行Trimmomatic
for forward_file in *_1.fastq.gz
do    
	base_name=${forward_file%_1.fastq.gz}    
	reverse_file="${base_name}_2.fastq.gz"    
	echo "正在用Trimmomatic处理样本：$base_name ..."
	trimmomatic PE -threads 2 \
		"$input_forward" "$input_reverse" \
		"${sample_name}_clean_1.fastq.gz" "${sample_name}_unpaired_1.fastq.gz" \
		"${sample_name}_clean_2.fastq.gz" "${sample_name}_unpaired_2.fastq.gz" \
		ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &
done
wait
echo "所有样本的Trimmomatic质控已完成！"