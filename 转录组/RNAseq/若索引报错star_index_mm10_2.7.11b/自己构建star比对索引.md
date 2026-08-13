STAR
    --runThreadN 20 \  # 使用 12 个线程并行计算
    --limitGenomeGenerateRAM 200000000000 \    # 内存上限设为 200GB
    --runMode genomeGenerate \                 # 指定当前模式为"生成基因组索引"
    --genomeDir ~/data/RNAseq/star_index_mm10_2.7.11b \   # 索引文件的输出目录
    --genomeFastaFiles ~/data/RNAseq/reference/GRCm38.primary_assembly.genome.fa \  # 基因组 FASTA 文件
    --sjdbGTFfile ~/data/RNAseq/reference/gencode.vM25.annotation.gtf                  # 基因注释 GTF 文件