#!/bin/bash
# multi_sample_derep.sh
# 合并所有样本的基因
# cat sample*_contigs.fasta > all_contigs_combined.fasta
# 测试
cat contigs.fasta > all_contigs_combined.fasta
# 基因预测
prodigal -i all_contigs_combined.fasta \
	-a all_genes.faa \
	-p meta \
	-q
# 去冗余（使用更严格的阈值）
# 90%相似度
cd-hit -i all_genes.faa \
	-o all_genes_nr.faa \
	-c 0.9 \
	-n 5 \
	-M 32000 \
	-T 16 \
	-g 1 \
	-d 0
echo "多样本合并去冗余完成！"
echo "总基因数: $(grep -c '^>' all_genes.faa)"
if [ -f all_genes_nr.faa ]; then
	echo "非冗余基因数: $(grep -c '^>' all_genes_nr.faa)"
else
	echo "错误: all_genes_nr.faa 未生成，cd-hit 执行失败！" >&2
	exit 1
fi
