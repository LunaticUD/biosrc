#!/bin/bash
# run_batch_assembly.sh
# 批量处理所有样本
for forward_file in *_clean_1.fastq.gz
do    
	# 获取样本名    
	base_name=${forward_file%_clean_1.fastq.gz}    
	reverse_file="${base_name}_clean_2.fastq.gz"
    echo "开始组装样本: $base_name"
    # 使用MEGAHIT组装    
    megahit \
    	-1 $forward_file \
    	-2 $reverse_file \
    	-o "${base_name}_megahit" \
    	--min-contig-len 500 \
    	-t 4 &
    # 如果也想用metaSPAdes，可以取消下面注释    
    # spades.py \
    #     --meta \
    #     -1 $forward_file \
    #     -2 $reverse_file \
    #     -o "${base_name}_metaspades" \
    #     -t 4 &
done
wait
echo "所有样本组装完成！"
