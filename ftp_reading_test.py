"""
date：2024-03-24
"""
import csv
from ftplib import FTP
import os

def list_files_ftp(ftp, directory):
    '''获得指定路径下的所有文件名称'''
    try:
        ftp.cwd(directory)  # 切换到指定目录
        files = ftp.nlst()  # 获取目录下的文件列表
        return files
    except Exception as e:
        print("Error: Unable to list files.")
        print(e)
        return []

def download_files_ftp(ftp, directory, localpath):
    '''下载文件到本地文件夹'''
    try:
        files = list_files_ftp(ftp, directory)
        if not files:
            return False
        # 获取FTP链接的最后一个部分作为文件夹名称
        folder_name = directory.rsplit('/', 1)[-1]
        # 创建新文件夹
        folder_path = os.path.join(localpath, folder_name)
        os.makedirs(folder_path, exist_ok=True)
        # 切换到新文件夹
        os.chdir(folder_path)
        # 下载文件
        for file in files:
            with open(os.path.join(folder_path, os.path.basename(file)), 'wb') as f:
                ftp.retrbinary('RETR ' + file, f.write)
        print("Downloaded files from directory '{}' successfully.".format(directory))
        return True
    except Exception as e:
        print("Error: Unable to download files from directory '{}'.".format(directory))
        print(e)
        return False

def ftp_connect(ftpserver, username, password):
    '''连接ftp服务器'''
    try:
        ftp = FTP(ftpserver)
        ftp.login(username, password)
        print("Connected to FTP server.")
        return ftp
    except Exception as e:
        print("Error: Unable to connect to FTP server.")
        print(e)
        return None

def ftp_disconnect(ftp):
    '''关闭ftp连接'''
    try:
        ftp.quit()
        print("Disconnected from FTP server.")
    except Exception as e:
        try:
            ftp.close()
            print("Disconnected from FTP server.")
        except Exception as e:
            print("Error: Unable to disconnect from FTP server.")
            print(e)

def read_ftp_links_from_csv(csv_file):
    '''解析csv下的ftp文件路径'''
    ftp_links = []
    try:
        with open(csv_file, newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                '''注意有两个ftp来源，一个genbank，一个refseq'''
                if row[13]:
                    ftp_link = row[13].replace('ftp://ftp.ncbi.nlm.nih.gov', '')  # 替换掉前缀部分
                    ftp_links.append(ftp_link)  # 加入列表
                elif row[14]:
                    ftp_link = row[14].replace('ftp://ftp.ncbi.nlm.nih.gov', '')  # 替换掉前缀部分
                    ftp_links.append(ftp_link)  # 加入列表
        return ftp_links
    except Exception as e:
        print("Error: Unable to read FTP links from CSV file.")
        print(e)
        return []

def main():
    # 输入参数
    ftpserver = 'ftp.ncbi.nlm.nih.gov'
    username = ''
    password = ''
    localpath = 'D:\CurrentClass\BioInfo\homework\homework-1\data_ftp'  # 本地保存路径
    csv_file = 'D:\CurrentClass\BioInfo\homework\homework-1\\viruses.csv'  # 读取文件路径

    # 读取表格中的目录列表
    file_links = []
    try:
        file_links = read_ftp_links_from_csv(csv_file)
    except Exception as e:
        print("Error: Unable to read directories from csv file.")
        print(e)
        return

    # 连接到FTP服务器
    ftp = ftp_connect(ftpserver, username, password)
    if not ftp:
        return

    # 逐个下载目录中的文件
    for file_link in file_links:
        download_files_ftp(ftp, file_link, localpath)

    # 断开与FTP服务器的连接
    ftp_disconnect(ftp)

if __name__ == "__main__":
    main()


