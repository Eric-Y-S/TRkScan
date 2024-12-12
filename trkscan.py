import os
import argparse
import pandas as pd
import edlib
from itertools import compress, accumulate
from stringDecompose import Decompose
from TRgenerator import TR_multiMotif
import numpy as np
import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

class BKTree:
    def __init__(self, dist_func):
        self.tree = {}
        self.dist_func = dist_func

    def add(self, node):
        if not self.tree:
            self.tree[node] = {}
        else:
            node_added = False
            for key, children in self.tree.items():
                dist = self.dist_func(key, node)
                if dist not in children:
                    children[dist] = node
                    node_added = True
                    break
            if not node_added:
                self.tree[node] = {}

    def query(self, node, max_distance):
        results = []
        for key, children in self.tree.items():
            dist = self.dist_func(key, node)
            if dist <= max_distance:
                results.append((key, dist))
            for dist, child in children.items():
                if dist <= max_distance:
                    results.append((child, dist))
        return results

def find_N(task):
    pid, total_task, sequence = task
    N_coordinates = []
    start = None
    for i, char in enumerate(sequence):
        if char == 'N' and start is None:
            start = i
        elif char != 'N' and start is not None:
            N_coordinates.append((start, i))
            start = None
    if start is not None:
        N_coordinates.append((start, len(sequence)))
    
    N_coordinate_df = pd.DataFrame(N_coordinates, columns=['start', 'end'])
    print(f'Process N Complete: {pid}/{total_task}')
    return ( pid, total_task, N_coordinate_df )

def decompose_sequence(task):
    pid, total_task, sequence, ksize = task
    seg = Decompose(sequence, ksize)
    seg.count_kmers()
    seg.build_dbg()
    seg.find_motif()
    print(f'Decomposition Complete: {pid}/{total_task}')
    return ( pid, total_task, seg )

def single_motif_annotation(task):
    pid, total_task, seg = task
    seg.annotate_with_motif()
    print(f'Single Motif Annotation Complete: {pid}/{total_task}')
    return ( pid, total_task, seg )

def rolling_same(seq1, seq2):
    if len(seq1) != len(seq2):
        return False
    seq1double = seq1 * 2
    return seq2 in seq1double

def get_distacne_and_cigar(ref, query):
    matches = edlib.align(query, ref, mode = "NW", task = "path")
    return matches['editDistance'], matches['cigar']

def rc(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def calculate_distance(ref, query_list):
    min_dist_list = []
    for query in query_list:
        if len(ref) >= len(query):
            tmp_ref, tmp_query = query, ref
        
        rep = len(tmp_query) // len(tmp_ref)
        extended_ref = tmp_ref * rep
        extended_ref_rc = rc(tmp_ref) * rep
        ### print(extended_ref, extended_ref_rc)
        min_dist1, min_dist2 = len(extended_ref), len(extended_ref_rc)
        for i in range(len(tmp_query)):
            rolling_query = tmp_query[i:] + tmp_query[:i]
            matches_ref = edlib.align(rolling_query, extended_ref, mode="NW", task="distance")
            min_dist1 = min(min_dist1, matches_ref['editDistance'])

            matches_ref_rc = edlib.align(rolling_query, extended_ref_rc, mode="NW", task="distance")
            min_dist2 = min(min_dist2, matches_ref_rc['editDistance'])
        min_dist_list.append(f'{min_dist1},{min_dist2}')
    return min_dist_list

def rotate_strings(s):
    n = len(s)
    return [s[i:] + s[:i] for i in range(n)]

if __name__ == "__main__":
    ##################################
    # configuration
    ##################################
    # DEFAULT: python trkscan_test.py testData/testData.fasta testData/testData.annotated.bed
    script_path = os.path.realpath(__file__)
    script_directory = os.path.dirname(script_path).replace('\\','/')

    parser = argparse.ArgumentParser(description='TRkScan')
    parser.add_argument('input', type = str, help = 'FASTA file you want to annotate')
    parser.add_argument('output', type = str, help = 'output prefix')
    parser.add_argument('-t', '--thread', type = int, default = 1, help = 'number of threads')
    parser.add_argument('-k', '--ksize', type = int, default = 5, help = 'k-mer size for building De Bruijn graph')
    parser.add_argument('-m', '--motif', type = str, default = f'{script_directory}/data/refMotif.txt', help='reference motif set')  ###### need to add reference set
    parser.add_argument('-f', '--force', action='store_true', help="annotate with motif X in given motif set no matter whether motif X is in the sequence")
    args = parser.parse_args()
    
    window_size = 5000
    step_size = 4000

    ##################################
    # generate test data
    ##################################
    if False:
        test_data = TR_multiMotif(['AATGG','CCATT','CTT'], 13000, 0.1, seed)
        test_seq = Seq(test_data.sequence)
        test_record = SeqRecord(test_seq, id = "testSeq", description = '')

        with open("testData/testData.fasta", "w") as handle:
            SeqIO.write(test_record, handle, "fasta")
        df = test_data.annotation
        df.to_csv('testData/testData_annotation.tsv', sep = '\t', index = False)
        df = test_data.annotation_woMut
        df.to_csv('testData/testData_annotation_woMut.tsv', sep = '\t', index = False)


    ##################################
    # read data
    ##################################
    with open(args.input, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    ##################################
    # read reference motif set
    ##################################
    tree = BKTree(Levenshtein.distance)
    with open(args.motif, "r") as motifDB:
        for line in motifDB:
            tree.add(line.strip())



    annotation_list = []
    for record in records:
        seq_name = record.name
        print(f'Start Processing [{seq_name}]')

        ##################################
        # deal with N character          
        ##################################
        seq = str(record.seq)
        seqLenwN = len(seq)
        total_task = max(1, int((seqLenwN - window_size) / step_size + 1) )
        with Pool(processes = args.thread) as pool:
            cur = 0
            tasks = []
            while cur < total_task:
                tasks.append(tuple([cur + 1, total_task, seq[cur*step_size : cur*step_size + step_size]]))
                cur += 1
            results = pool.imap_unordered(find_N, tasks)
            pool.close()    # close the pool and don't receive any new tasks
            pool.join()     # wait for all the tasks to complete
            results = list(results)
        
        merged_list = []
        for pid, _, result in results:
            result[['start', 'end']] += (pid - 1) * window_size
            merged_list.append(result)

        N_coordinate = pd.concat(merged_list, ignore_index=True)
        N_coordinate = N_coordinate.sort_values(by=['start']).reset_index(drop=True)

        if N_coordinate.shape[0] > 0:
            N_coordinate_new = pd.DataFrame(columns=['start', 'end'])

            start = N_coordinate.loc[0, 'start']
            end = N_coordinate.loc[0, 'end']
            for i in range(1, N_coordinate.shape[0]):
                if N_coordinate.loc[i, 'start'] != end:
                    N_coordinate_new = pd.concat([N_coordinate_new, pd.DataFrame({'start': [start], 'end': [end]})], ignore_index=True)
                    start = N_coordinate.loc[i, 'start']
                end = N_coordinate.loc[i, 'end']

            # add the last region
            N_coordinate_new = pd.concat([N_coordinate_new, pd.DataFrame({'start': [start], 'end': [end]})], ignore_index=True)

            # solve with 'old_index' and 'new_index'
            N_coordinate_new['old_index'] = N_coordinate_new['end']
            N_coordinate_new['new_index'] = N_coordinate_new['end']

            # calculate new index
            l = 0
            for i in range(N_coordinate_new.shape[0]):
                l += N_coordinate_new.loc[i,'end'] - N_coordinate_new.loc[i,'start']
                N_coordinate_new.loc[i,'new_index'] -= l

            ### print(N_coordinate_new)  #############################

        ##################################
        # split into windows to decompose
        ##################################
        mask = [base != "N" for base in seq]
        filter_seq = ''.join(compress(seq, mask))
        seqLenwoN = len(filter_seq)
        total_task = max(1, int((seqLenwoN - window_size) / step_size + 1) )

        ##################################
        # get motifs
        ##################################
        with Pool(processes = args.thread) as pool:
            cur = 0
            tasks = []
            while cur < total_task:
                tasks.append(tuple([cur + 1, total_task, filter_seq[cur*step_size : cur*step_size + window_size], args.ksize]))
                cur += 1
            results = pool.imap_unordered(decompose_sequence, tasks)
            
            pool.close()    # close the pool and don't receive any new tasks
            pool.join()     # wait for all the tasks to complete
            results = list(results)

        ##################################
        # merge and decide the representive of motifs
        ##################################
        motifs, nondup = [], []
        for pid, _, result in results:
            motifs.extend(result.motifs['motif'].to_list())
        for idx in range(len(motifs)):
            is_dup = False
            for motif in nondup:
                if rolling_same(motifs[idx], motif):
                    is_dup = True
                    break
            if not is_dup:
                nondup.append(motifs[idx])

        nondup_adjusted = []
        for motif in nondup:
            min_distance = 1e6
            best_motif = motif
            sequences = rotate_strings(motif)
            max_distance = int(0.5 * len(motif))
            for seq in sequences:
                ### print(f'xxxxx:{seq}')
                matches = tree.query(seq, max_distance)
                for ref, dist in matches:
                    ### print(ref, dist)
                    if dist < min_distance:
                        best_motif = seq
                        min_distance = dist
            nondup_adjusted.append(best_motif)
            ### print(best_motif)

        print(nondup_adjusted)

        for pid, _, result in results:
            result.motifs_list = nondup_adjusted

        ##################################
        # annotate single motif
        ##################################
        with Pool(processes = args.thread) as pool:
            cur = 0
            tasks = results
            results = pool.imap_unordered(single_motif_annotation, tasks)
            pool.close()
            pool.join()
            results = list(results)
        
        # merge single motif annotation
        merged_list = []
        for pid, _, result in results:
            df = result.annotation
            df[['start', 'end']] += (pid - 1) * step_size
            merged_list.append(df)

        merged_df = pd.concat(merged_list, ignore_index=True)
        merged_df = merged_df.drop_duplicates().sort_values(by='end').reset_index(drop=True)
        
        del results
        
        ##################################
        # DP to link
        ##################################
        # parameters
        gap_penalty = 1
        distance_penalty = 1
        perfect_bonus = 0.5

        length = len(filter_seq)
        dp = np.zeros(length + 1, dtype=np.float64)  # dp[i] means sep[:i] score sum
        pre = np.full((length + 1, 3), None)  # (pre_i, motif_id, motif)
        idx = 0
        for i in range(1, length + 1):
            # skip one base
            if dp[i-1] - gap_penalty  > 0:
                dp[i] = dp[i-1] - gap_penalty 
                pre[i] = (i-1, None, None)
            
            while idx < merged_df.shape[0] and merged_df.loc[idx, 'end'] <= i:
                if merged_df.loc[idx, 'end'] < i:
                    idx += 1
                    continue
                
                pre_i = merged_df.loc[idx, 'start']
                motif = merged_df.loc[idx,'motif']
                distance = merged_df.loc[idx,'distance']
                bonus = perfect_bonus * len(motif) if distance == 0 else 0
                
                if dp[pre_i] + len(motif) - distance * distance_penalty + bonus >= dp[i]:
                    dp[i] = dp[pre_i] + len(motif) - distance * distance_penalty + bonus
                    pre[i] = (pre_i, idx, motif)

                idx += 1

            if i % 5000 == 0:
                print(f'DP: {i // 1000} kbp is Done!')

        print('DP complete!')

        # retrace
        idx = length
        next = [(None, None) for i in range(length+1)] # (next coordinate, reference motif)
        while idx >= 0:
            if pre[idx][0] != None:
                if pre[idx][1] != None:        # match a motif
                    next[pre[idx][0]] = (idx, merged_df.loc[pre[idx][1],'motif'])
                    idx = pre[idx][0]
                else:                          # skip a base
                    next[idx - 1] = (idx, None)
                    idx -= 1
            else:
                idx -= 1
  
        print('Retrace Complete!')

        # output motif annotation file
        annotation_data = []
        idx = 0
        
        while idx < length:
            if pre[idx][0] == None and next[idx][0] != None:    # is a start site
                idx2 = idx
                start, end, cur_motif, rep_num, score, max_score, cigar_string  = idx, 0, None, 0, 0, 0, ''
                skip_num = 0
                while next[idx2][0] != None:
                    if next[idx2][1]:    # match a motif 
                        ### print(next[idx2][1])   ###########################
                        if cur_motif != next[idx2][1]:       # split the annotation and init
                            ### print(f'{cur_motif} -> {next[idx2][1]}')
                            if cur_motif != None:
                                annotation_data.append(row)
                                start, end, cur_motif, rep_num, score, max_score, cigar_string  = idx2, idx2, next[idx2][1], 0, 0, 0, ''
                                skip_num = 0
                            else:
                                cur_motif = next[idx2][1]

                        if skip_num:
                            cigar_string += f'{skip_num}N'
                            score -= skip_num * gap_penalty
                            skip_num = 0

                        ### print(rep_num)
                        rep_num += 1
                        distance, cigar = get_distacne_and_cigar(next[idx2][1], filter_seq[idx2 : next[idx2][0]])
                        score += len(cur_motif) - distance * distance_penalty
                        cigar_string += cigar + '/'
                        if score >= max_score:
                            row = [seq_name, seqLenwN, start, next[idx2][0], cur_motif, rep_num, score, cigar_string]
                            max_score = score
                    else:               # skip a base
                        if idx2 != 0:
                            ### print(idx2)
                            skip_num += 1
                
                    idx2 = next[idx2][0]
                    end = idx2
                annotation_data.append(row)

            idx += 1

        annotation = pd.DataFrame(annotation_data, columns=['seq','length','start','end','motif','rep_num','score','CIGAR'])
        annotation = annotation.reset_index(drop=True)
        
        ### print(annotation)
        
        # coordinates transformation with N character
        if N_coordinate.shape[0]:
            for i in range(annotation.shape[0]):
                # edit cigar string
                old_cigar = annotation.loc[i, 'CIGAR']
                new_cigar, cot = '', ''
                idx = 0
                cur = annotation.loc[i, 'start']
                while idx < len(old_cigar):
                    symbol = old_cigar[idx]
                    if symbol == '/':
                        new_cigar += '/'
                    elif symbol in ['=','X','I','N']:
                        length = int(cot)
                        tmp = N_coordinate_new[(N_coordinate_new['new_index'] > cur) & (N_coordinate_new['new_index'] <= cur + length)]
                        if tmp.shape[0]:    # if there is N character in this region
                            ### print(f'###{cur} ~ {cur+length}')
                            ### print(tmp)
                            
                            for row in tmp.itertuples():
                                length -= row.new_index - cur
                                new_cigar += f'{row.new_index - cur}{symbol}'
                                new_cigar += f'{row.end - row.start}N'
                                cur += row.new_index - cur
                                ### new_cigar += f'{N_coordinate_new['new_index'] - cur}{symbol}'
                            new_cigar += f'{length}{symbol}' if length != 0 else ''
                            cur += length
                        else:
                            new_cigar += f'{length}{symbol}'
                            cur += length
                        cot = ''
                    elif symbol == 'D':
                        new_cigar += f'{cot}{symbol}'
                        cot = ''
                    else:
                        cot += symbol
                    
                    annotation.loc[i, 'CIGAR'] = new_cigar
                    idx += 1
                
                
                # transform start site
                tmp = N_coordinate_new[N_coordinate_new['new_index'] <= int(annotation.loc[i,'start'])]
                if tmp.shape[0]:
                    tmp = tmp.reset_index(drop=True)
                    l = tmp.shape[0] - 1
                    annotation.loc[i,'start'] += tmp.loc[l,'old_index'] - tmp.loc[l,'new_index']
                # transform end site
                tmp = N_coordinate_new[N_coordinate_new['new_index'] <= int(annotation.loc[i,'end'])]
                if tmp.shape[0]:
                    tmp = tmp.reset_index(drop=True)
                    l = tmp.shape[0] - 1
                    annotation.loc[i,'end'] += tmp.loc[l,'old_index'] - tmp.loc[l,'new_index'] 
        ### print(annotation)
        annotation_list.append(annotation)
        
    merged_annotation = pd.concat(annotation_list, ignore_index=True)
    merged_annotation.to_csv(f'{args.output}.concise.tsv', sep = '\t', index = False)

    ##################################
    # make file - *.annotation.tsv
    ##################################
    detail_annotation = []
    for record in records:
        seq_name = record.name
        sequence = str(record.seq)
        seqLen = len(sequence)
        tmp = merged_annotation[merged_annotation['seq'] == seq_name]
        ### print(tmp)
        
        for index, row in tmp.iterrows():
            idx = row.start
            cigar_string = row.CIGAR
            j = 0
            start, length, actual_motif, sub_cigar, num = row.start, 0, '', '', ''
            while j < len(cigar_string):
                symbol = cigar_string[j]
                if symbol == '/':
                    detail_annotation.append([seq_name, seqLen, start, start + length, row.motif, actual_motif, sub_cigar])
                    start = start + length
                    length, actual_motif, sub_cigar = 0, '', ''
                elif symbol in ['=','X','I']:
                    actual_motif += sequence[start + length : start + length + int(num)]
                    length += int(num)
                    sub_cigar += f'{num}{symbol}'
                    num = ''
                elif symbol == 'D':
                    sub_cigar += f'{num}D'
                    num = ''
                elif symbol == 'N':
                    length += int(num)
                    sub_cigar += f'{num}N'
                    num = ''
                else:   # symbol == number
                    num += symbol
                j += 1
            ### detail_annotation.append([seq_name, seqLen, start, start + length, row.motif, actual_motif, sub_cigar])

    detailed_df = pd.DataFrame(data = detail_annotation, columns = ['seq','length','start','end','motif','actual_motif','CIGAR'])
    detailed_df.to_csv(f'{args.output}.annotation.tsv', sep = '\t', index = False)


    ##################################
    # make file - *.motif.tsv
    ##################################
    rep_num = []
    motifs_list = list(set(merged_annotation['motif'].to_list()))
    for motif in motifs_list:
        tmp = merged_annotation[merged_annotation['motif'] == motif]
        rep_num.append(sum(tmp['rep_num']))

    data = {'motif': motifs_list, 'rep_num': rep_num}
    motif_df = pd.DataFrame(data = data)

    # calculate distance
    for motif in motifs_list:
        motif_df[motif] = calculate_distance(motif, motifs_list)
    motif_df.to_csv(f'{args.output}.motif.tsv', sep = '\t', index = False)





