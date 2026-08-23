#!/bin/bash
sample_name="sample"
input_forward="${sample_name}_1.fastq.gz"
input_reverse="${sample_name}_2.fastq.gz"
# Step 1: 使用FastQC进行质控前评估
echo "Step 1: 运行FastQC进行原始数据质量评估..."
fastqc "$input_forward" "$input_reverse" -t 2 # -t 指定线程数
# Step 2: 使用Trimmomatic进行质控
echo "Step 2: 运行Trimmomatic进行数据质控..."
trimmomatic PE -threads 4 \
    "$input_forward" "$input_reverse" \
    "${sample_name}_clean_1.fastq.gz" "${sample_name}_unpaired_1.fastq.gz" \
    "${sample_name}_clean_2.fastq.gz" "${sample_name}_unpaired_2.fastq.gz" \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50    
    # 参数解释：    
    # PE: 双端数据    
    # ILLUMINACLIP: 切除接头。TruSeq3-PE-2.fa是适配器序列文件。    
    # LEADING:3 从 reads 的开头切除质量值低于3的碱基    
    # TRAILING:3 从 reads 的末尾切除质量值低于3的碱基    
    # SLIDINGWINDOW:4:20 使用4碱基宽的滑动窗口，当窗口内平均质量值低于20时，切除其后所有序列    
    # MINLEN:50 丢弃长度小于50bp的reads
echo "Trimmomatic质控完成！"
echo "高质量配对数据：${sample_name}_clean_1.fastq.gz 和 ${sample_name}_clean_2.fastq.gz"
echo "未配对的数据被单独存放。"