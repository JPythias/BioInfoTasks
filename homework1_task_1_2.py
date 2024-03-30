import requests

def download_fasta_from_uniprot(uniprot_ids, local_directory):
    try:
        for uniprot_id in uniprot_ids:
            # 构建 UniProt 的 FASTA 文件下载链接
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"

            # 发送 GET 请求获取 FASTA 文件内容
            response = requests.get(url)

            # 检查请求是否成功
            if response.status_code == 200:
                # 构建本地保存路径
                local_file_path = f"{local_directory}/{uniprot_id}.fasta"
                # 将 FASTA 文件内容写入本地文件
                with open(local_file_path, 'w') as f:
                    f.write(response.text)
                print(f"FASTA file for UniProt ID {uniprot_id} downloaded successfully to {local_file_path}")
            else:
                print(f"Failed to fetch the FASTA file for UniProt ID {uniprot_id}!")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    # 输入 UniProt ID 列表和本地保存目录
    uniprot_ids = ['A0A0A7HFG9', 'A9JKQ6', 'A9JKQ8', 'A9JKR0', 'A9JKR2',
                   'A9JKR4', 'A9JKR5', 'A9JKR7', 'A9JKR9', 'A9JKS1']  # 这里是示例 UniProt ID 列表，你可以根据需要添加或修改
    local_directory = 'D:\CurrentClass\BioInfo\homework\homework-1\data2'  # 本地保存目录

    # 调用函数下载 FASTA 文件
    download_fasta_from_uniprot(uniprot_ids, local_directory)

if __name__ == "__main__":
    main()
