"""
date：2024-03-25
"""
import os
import csv
from Bio import Entrez

def download_replicon(replicon_id, output_folder):
    """下载指定 replicon_id 的 replicon 文件"""
    Entrez.email = "jawnpythias@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=replicon_id, rettype="gb", retmode="text")
    replicon_content = handle.read().encode()  # 将字符串转换为字节类型
    with open(os.path.join(output_folder, f"{replicon_id}.gb"), "wb") as f:  # 使用二进制模式写入文件
        f.write(replicon_content)
    handle.close()

def download_replicons_from_csv(csv_file, replicon_column, output_folder):
    """从 CSV 文件中读取 replicon 列的值，并下载相应的 replicon 文件"""
    # 创建存储 replicon 文件的文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 读取 CSV 文件，并下载 replicon 文件
    with open(csv_file, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            replicon_id = row[replicon_column]
            download_replicon(replicon_id, output_folder)
            print(f"Replicon {replicon_id} 已下载到 {output_folder} 文件夹中。")

def main():
    csv_file = "D:\CurrentClass\BioInfo\homework\homework-1\\viruses.csv"
    replicon_column = "Replicons"
    output_folder = "D:\CurrentClass\BioInfo\homework\homework-1\data"

    # 下载 replicon 文件
    download_replicons_from_csv(csv_file, replicon_column, output_folder)

if __name__ == "__main__":
    main()
