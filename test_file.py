from Bio import Entrez
from Bio import SeqIO

def search_download_genome(genus_species):
    """使用NCBI的Entrez工具搜索并下载指定物种的基因组信息"""
    Entrez.email = "jawnpythias@gmail.com"  # 设置你的邮箱地址
    query = f"{genus_species}[Organism] AND complete genome[Title]"
    handle = Entrez.esearch(db="nucleotide", term=query)
    record = Entrez.read(handle)
    handle.close()
    if record['Count'] == '0':
        print(f"No complete genome found for {genus_species}.")
        return
    genome_id = record['IdList'][0]
    handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
    return handle.read()

def extract_gabija_sequence(genome_content):
    """从基因组文件中提取Gabija基因序列"""
    for record in SeqIO.parse(genome_content, "genbank"):
        for feature in record.features:
            if 'gene' in feature.qualifiers and 'Gabija' in feature.qualifiers['gene'][0]:
                print(f"Found Gabija gene sequence: {feature.qualifiers['gene'][0]}")
                return feature.extract(record.seq)
    print("Gabija gene sequence not found.")
    return None

def main():
    # 搜索并下载VD045细菌基因组信息
    genus_species = "VD045"
    genome_content = search_download_genome(genus_species)
    if not genome_content:
        return

    # 提取Gabija基因序列
    gabija_sequence = extract_gabija_sequence(genome_content)
    if gabija_sequence:
        print("Gabija gene sequence:")
        print(gabija_sequence)

if __name__ == "__main__":
    main()
