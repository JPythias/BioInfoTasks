"""
date：2024-03-27
"""
'''
使用非洲猪瘟的A9JKQ6和A9JKQ8
'''
from Bio import SeqIO

MIN = -float("inf")
A = -10
B = -1

def read_fasta(file_path):
    '''直接用了Bio包读fasta文件的方法，返回的值是一个list，并且只有一个元素'''
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def read_blosum():
    '''读取Blosum62矩阵，定义打分方法'''
    with open(r'D:\CurrentClass\BioInfo\PyProject\RefData\BLOSUM62.txt') as f:
        fileContent = f.readlines()
    global alignment_score # 定义全局变量打分矩阵值，下面的处理是让矩阵不带氨基酸索引
    alignment_score = []
    for i in range(1, len(fileContent)):
        # 去掉第一行的氨基酸索引
        alignment_score = alignment_score + [fileContent[i].split()]
    alignment_score = [line[1:] for line in alignment_score] # 去掉第一列的氨基酸索引
    global alignmentIndices # 定义氨基酸序号的索引
    alignmentIndices = fileContent[0].split()
    return fileContent

def get_score(a, b):
    '''获取打分'''
    indice_a = -1
    indice_b = -1
    for i in range(0, len(alignmentIndices)):
        if (a == alignmentIndices[i]):
            indice_a = i
        if (b == alignmentIndices[i]):
            indice_b = i
    return int(alignment_score[indice_a][indice_b])

def initialisation_x(i, j):
    '''初始化E'''
    if i >0 and j == 0:
        # 在最左列且不是第一行，初始化为负无穷
        return MIN
    else:
        if j > 0:
            # 除了最左列的其他列，空位惩罚
            return -10+(-0.5*j)
        else:
            # 第一行，初始化为0
            return 0

def initialisation_y(i, j):
    '''初始化F，与上同理'''
    if j > 0 and i == 0:
        return MIN
    else:
        if i > 0:
            return -10 + (-0.5 * i)
        else:
            return 0

def initialisation_m(i, j):
    '''初始化G，除了边缘其余都为负无穷'''
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return MIN
        else:
            return 0

def matrix_distance(s, t):
    '''计算矩阵距离，NW算法'''
    E = [[initialisation_x(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]
    F = [[initialisation_y(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]
    G = [[initialisation_m(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]

    for j in range(1, len(s)+1):
        for i in range(1, len(t)+1):
            E[i][j] = max((A + B + G[i][j - 1]), (B + E[i][j - 1]), (A + B + F[i][j - 1]))
            F[i][j] = max((A + B + G[i - 1][j]), (A + B + E[i - 1][j]), (B + F[i - 1][j]))
            G[i][j] = max(get_score(t[i - 1], s[j - 1]) + G[i - 1][j - 1], E[i][j], F[i][j])
    return [E, F, G]

def backtracking(s, t, E, F, G):
    '''回溯最佳比对序列'''
    seq1 = ''
    seq2 = ''
    i = len(t)
    j = len(s)
    # 右下角的分数
    score = G[i][j]
    while (i>0 or j>0):
        if (i>0 and j>0 and G[i][j] == G[i-1][j-1] + get_score(t[i - 1], s[j - 1])):
            seq1 += s[j-1]
            seq2 += t[i-1]
            i -= 1
            j -= 1
        elif (i>0 and G[i][j] == F[i][j]):
            seq1 += '_'
            seq2 += t[i-1]
            i -= 1
        elif (j>0 and G[i][j] == E[i][j]):
            seq1 += s[j-1]
            seq2 += '_'
            j -= 1
    seq1r = ' '.join([seq1[j] for j in range(-1, -(len(seq1)+1), -1)])
    seq2r = ' '.join([seq2[j] for j in range(-1, -(len(seq2)+1), -1)])
    return [seq1r, seq2r, score]

def write_strings_to_file(seq1r, seq2r, file_path):
    '''写入比对结果'''
    with open(file_path, "w") as file:
        file.write(seq1r + "\n")
        file.write(seq2r + "\n")

def main():
    # 1读取fasta文件，转变为string类型
    seq1 = read_fasta(r"D:\CurrentClass\BioInfo\homework\homework-1\data2\A9JKQ6.fasta")
    s = seq1[0]
    seq2 = read_fasta(r"D:\CurrentClass\BioInfo\homework\homework-1\data2\A9JKQ8.fasta")
    t = seq2[0]

    # 2. 获取Blosum62矩阵
    read_blosum()

    # 3.进行比对，写入结果
    [E, F, G] = matrix_distance(s, t)
    [seq1r, seq2r, score] = backtracking(s, t, E, F, G)
    print(seq1r)
    print(seq2r)
    print(score)
    write_strings_to_file(seq1r, seq2r, "D:\CurrentClass\BioInfo\homework\homework-2\\result.txt")


if __name__ == "__main__":
    main()
