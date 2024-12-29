import argparse
import sourmash
import re
import numpy as np
from Bio import SeqIO
from itertools import compress
from collections import Counter
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

# DEFAULT: python trkscan_test.py testData/testData.fasta testData/testData.annotated.bed
parser = argparse.ArgumentParser(description='get_k_parameter')
parser.add_argument('input', type=str, help='FASTA file you want to annotate')
parser.add_argument('-s', '--sample', type=float, default=1, help='ratio to sample')
parser.add_argument('-w', '--window_size', type=int, default=1000, help='Window size for Fourier transform')
args = parser.parse_args()

fasta_path = args.input
sample_rate = args.sample
window_size = args.window_size

class get_likely_minimun_repeat_length:
    def __init__(self, fasta_path, sample_rate, ksize = [3, 4, 9, 11, 15, 31, 41], min_length=10e3, max_length=50e3):
        self.fasta_path = fasta_path
        self.sample_rate = sample_rate
        self.min_length = min_length
        self.max_length = max_length
        self.win = window_size
        self.sampled_sequence = self.sample()
        self.likely_k = self.get_likely_repeat_length(self.sampled_sequence, ksize)
        print(self.likely_k)

    def sample(self):
        with open(self.fasta_path, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))

        for record in records:
            seq = str(record.seq)
            mask = [base != "N" for base in seq]
            seq = ''.join(compress(seq, mask))
            seq_length = len(seq)
            sample_num = max(int(seq_length * self.sample_rate / window_size) + 1, 1)

            sampled_sequence = []

            if seq_length - window_size + 1 > 0:
                start = np.random.randint(0, seq_length - window_size + 1, size = sample_num)
                sampled_sequence.extend( [seq[s : s + window_size] for s in start] )
            else:
                sampled_sequence.append(seq)
                #print(sampled_sequence)
        return sampled_sequence

    def get_likely_repeat_length(self, sequence, ksize):
        
        for k in ksize:
            for seq in sequence:
                mh = sourmash.MinHash(n = 5000, ksize = k + 1)
                kmers = [kmer for kmer, _ in mh.kmers_and_hashes(seq, force=True) ]
                kmer_count = Counter(kmers)
                # filter
                max_count = kmer_count.most_common(1)[0][1]
                min_count = int(max_count * 0.5)
                kmer_count = [ [k, v] for k, v in kmer_count.items() if v >= min_count ]
                # get key kmer
                kmer_count.sort(key=lambda item: item[1], reverse=True)
                kmer_count = kmer_count[0:5]
                kmers = [kmer[0] for kmer in kmer_count]
                # get index
                for kmer in kmers:
                    matches = re.finditer(kmer, seq)
                    index = np.array([match.start() for match in matches])
                    first_diff = np.diff(index, n=1)
                    print(first_diff)
            

        

            # 使用 DBSCAN 聚类
            dbscan = DBSCAN(eps=1, min_samples=2)  # eps 是距离阈值，min_samples 是每个簇的最小样本数
            dbscan.fit(data)

            # 输出簇标签，-1 表示噪声
            print("聚类标签:", dbscan.labels_)

            # 可视化结果
            plt.scatter(data, [0]*len(data), c=dbscan.labels_)
            plt.show()

        return None
    

if __name__ == "__main__":
    # Example usage
    sequence = args.input  # Replace with your actual file path
    sample_rate = 0.01
    detected_repeats = get_likely_minimun_repeat_length(sequence, sample_rate)
