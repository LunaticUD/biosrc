#!/bin/bash
# 定义输入文件和输出文件前缀
input_forward="sample_1.fastq.gz"
input_reverse="sample_2.fastq.gz"
output_prefix="clean_"
# 运行fastp
fastp -i "$input_forward" -I "$input_reverse" \
      -o "${output_prefix}${input_forward}" \
      -O "${output_prefix}${input_reverse}" \
      -h "${output_prefix}fastp_report.html" \
      -j "${output_prefix}fastp_report.json" \
      --thread 4
# 使用4个CPU线程加速
echo "Fastp质控完成！报告见 ${output_prefix}fastp_report.html"
