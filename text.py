import Levenshtein
from pybktree import BKTree
########################
# data/refMotif.txt:
# ATTGG
# ATTGGATTCG
# AATCC
# ATGGAAATCTCTACACCGCTATCTCTATGAAC
# GTTCATAGAGATAGCGGTGTAGAGATTTCCAT
########################

def rotate_strings(s):
    n = len(s)
    return [s[i:] + s[:i] for i in range(n)]

tree = BKTree(Levenshtein.distance)  # 使用 Levenshtein 距离

# 读取数据并添加到 BKTree
with open('data/refMotif.txt', 'r') as motifDB:
    for line in motifDB:
        tree.add(line.strip())  # 向树中添加每个motif

print

# 查询和调整
nondup = ['TGGAA', 'TGGAAA', 'ATGGAAATCTCTACACCGCTATCTCTATGAAC', 
          'GCTATCTCTATGAACATGGAAATCTCTACACC', 'TTTCTACACCGCTATCTCTATGAACATGGAAATC']

nondup_adjusted = []
for motif in nondup:
    min_distance = 1e6
    best_motif = motif
    sequences = rotate_strings(motif)  # 获取所有旋转的变体
    max_distance = int(0.5 * len(motif))  # 设置最大允许的距离
    for seq in sequences:
        ### print(f'############## 当前查询: {seq}')
        matches = tree.find(seq, max_distance)  # 查询匹配的motif
        ### print(matches)
        for dist, ref  in matches:
            ### print(f'匹配: {ref} 距离: {dist}')
            if dist < min_distance:
                best_motif = seq
                min_distance = dist
    nondup_adjusted.append(best_motif)

print("调整后的motif:", nondup_adjusted)