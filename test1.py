"""
date：2024-03-26
"""

from Bio import Entrez, SeqIO

# 设置您的邮箱地址（必须）
Entrez.email = "jawnpythias@gmail.com"

# 根据 GenBank 记录的 accession 号下载 GenBank 文件
accession_number = "GCF_000291215.1"
handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")

# 解析 GenBank 文件并找到 Gabija 基因的序列
for record in SeqIO.parse(handle, "genbank"):
    for feature in record.features:
        if feature.type == "gene" and "Gabija" in feature.qualifiers.get("gene", ""):
            gabija_gene_sequence = feature.location.extract(record.seq)
            print("Gabija Gene Sequence:")
            print(gabija_gene_sequence)
            break  # 找到第一个 Gabija 基因即停止搜索

handle.close()
