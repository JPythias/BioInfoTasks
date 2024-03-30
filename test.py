
"""
date：2024-03-26
"""

from Bio import SeqIO
import gzip

# 从GenBank文件中提取Gabija基因的注释信息和序列
def extract_gabija_gene(genbank_file, gene_name):
    gabija_gene_record = None
    with gzip.open(genbank_file, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene" and feature.qualifiers.get("gene", [""])[0] == gene_name:
                    gabija_gene_record = feature
                    break
    return gabija_gene_record

# 下载的GenBank文件路径
genbank_file = "D:\CurrentClass\BioInfo\homework\homework-1\\GCF_000291215.1_Baci_cere_VD045_V1_genomic.gbff.gz"
# Gabija基因名称
gabija_gene_name = "GajA"

# 提取Gabija基因的注释信息和序列
gabija_gene_record = extract_gabija_gene(genbank_file, gabija_gene_name)
if gabija_gene_record:
    print("Gabija gene found:")
    print(gabija_gene_record)
    print("Gabija gene sequence:")
    print(gabija_gene_record.extract(record.seq))
else:
    print("Gabija gene not found in the GenBank file.")
