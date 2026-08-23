#!/bin/bash
db_path="db"
for forward_file in *_1.fastq.gz
do    
	base_name=${forward_file%_1.fastq.gz}    
	reverse_file="${base_name}_2.fastq.gz"    
	output_dir="kneaddata_output_${base_name}" 
	# 为每个样本创建独立的输出目录    
	echo "正在用KneadData处理样本：$base_name ..."    
	kneaddata -i1 “$forward_file” -i2 "$reverse_file" \
		-o $output_dir \
		-db "$db_path/Homo_sapiens_Bowtie2_v1" \
		-t 2 \
		--trimmomatic-options "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50" \
		--remove-intermediate-output &
done
wait
echo "所有样本的Kneaddata分析已完成！"