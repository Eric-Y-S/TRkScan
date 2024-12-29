import argparse
import numpy as np
from Bio import SeqIO
from itertools import compress
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
    def __init__(self, fasta_path, sample_rate, window_size, ksize = [1, 2, 3, 4, 9, 11, 15, 31, 41], min_length=10e3, max_length=50e3):
        self.fasta_path = fasta_path
        self.sample_rate = sample_rate
        self.min_length = min_length
        self.max_length = max_length
        self.win = window_size
        self.sampled_sequence = self.sample()
        self.likely_k = self.sliding_window_fourier(self.sampled_sequence, self.win, ksize)
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

    def sequence_to_numeric(self, sequence, k):
        charaterSet = list(set([char for char in sequence]))
        clen = len(charaterSet)
        mapping = {c : id for id, c in enumerate(charaterSet)}
        signal = []
        for i in range(len(sequence) - k):
            subseq = sequence[i: i + k]
            value = 0
            for j in range(k):
                value *= clen
                value += mapping[subseq[j]]
            signal.append(value)
        print(np.array(signal))
        return np.array(signal)

    def find_repeats(self, sequence, k):
        # Convert sequence to numeric form
        numeric_seq = self.sequence_to_numeric(sequence, k)
        
        # Compute the Fourier transform of the sequence
        magnitudes = np.abs(np.fft.fft(numeric_seq))

        # Consider only the first half of the spectrum (mirrored)
        half_spectrum = magnitudes[:len(magnitudes) // 2]

        # Apply smoothing to the frequency spectrum
        half_spectrum = np.convolve(half_spectrum, np.ones(5) / 5, mode='valid')
        
        #n = len(numeric_seq)
        #interval = 1
        #frequency = np.arange(n / 2) / (n * interval) 
        #plt.plot(frequency, half_spectrum, 'red') 
        #plt.xlabel('Freq (Hz)'), plt.ylabel('Amp. Spectrum') 
        #plt.show() 
        
        # Find the index of the highest peak (ignoring the DC component)
        repeat_index = np.argmax(half_spectrum[1:]) + 1  # Exclude DC component
        
        # Calculate the repeat length based on the frequency
        repeat_length = len(sequence) / repeat_index
        return repeat_length
    
    def sliding_window_fourier(self, sequence, window_size, ksize):
        repeat_lengths_from_k = []
        for k in ksize:
            repeat_lengths_from_window = []

            for seq in sequence:
                tmp = self.find_repeats(seq, k)
                repeat_lengths_from_window.append(tmp)
            print(repeat_lengths_from_window)
            median = np.median(repeat_lengths_from_window)
            repeat_lengths_from_k.append(round(median, 1))
        return repeat_lengths_from_k

if __name__ == "__main__":
    # Example usage
    sequence = args.input  # Replace with your actual file path
    sample_rate = 0.01
    detected_repeats = get_likely_minimun_repeat_length(sequence, sample_rate, window_size)
