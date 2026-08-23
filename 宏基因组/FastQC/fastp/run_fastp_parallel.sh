# !/bin/bash
# 文件名命名规则一致：*_1.fastq.gz 和 *_2.fastq.gz
for forward_file in *_1.fastq.gz
do
	base_name=${forward_file%_1.fastq.gz}
	reverse_file="${base_name}_2.fastq.gz"
	output_prefix="${base_name}_clean_"
	echo "处理样本：$base_name"
	fastp -i "$forward_file" -I "$reverse_file" \
          -o "${output_prefix}1.fastq.gz" \
          -O "${output_prefix}2.fastq.gz" \
          -h "${output_prefix}fastp_report.html" \
          -j "${output_prefix}fastp_report.json" \
          --thread 2 &          
          # 后台运行，每个任务用2个线程
done
wait
echo "所有样本Fastp质控已完成！"