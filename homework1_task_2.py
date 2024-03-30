from Bio import SeqIO

# 打开下载的基因组数据文件（FASTA格式）
genome_file = "DataSet/sequenceVD045.fasta"

# 创建一个空列表来存储Gabija基因序列
gabija_sequences = []

# 使用Biopython的SeqIO模块打开基因组数据文件
for record in SeqIO.parse(genome_file, "fasta"):
    # 在注释信息中查找Gabija基因
    if ("GajA" in record.description) or ("GajB" in record.description):
        # 将找到的Gabija基因序列添加到列表中
        gabija_sequences.append(record)

# 打印或保存Gabija基因序列
for seq_record in gabija_sequences:
    print(seq_record.id)
    print(seq_record.seq)

# 如果需要，将Gabija基因序列保存为新的FASTA文件
SeqIO.write(gabija_sequences, "gabija_sequences.fasta", "fasta")