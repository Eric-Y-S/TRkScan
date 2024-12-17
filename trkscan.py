import os
import argparse
import subprocess
import pandas as pd
import edlib
from pybktree import BKTree
from itertools import compress
from stringDecompose import Decompose
import numpy as np
import Levenshtein
from Bio import SeqIO   # I/O processing
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool    # multi-thread
import time

def find_N(task):
    pid, total_task, sequence = task
    N_coordinates = []
    start = None

    sequence = np.array(list(sequence))
    n_indices = np.where(sequence == 'N')[0]
    if len(n_indices) == 0:
        print(f'Process N Complete: {pid}/{total_task}')
        return pd.DataFrame(columns=['start', 'end'])
    
    for i in range(len(n_indices)):
        if start is None:
            start = n_indices[i]
        if i == len(n_indices) - 1 or n_indices[i+1] != n_indices[i] + 1:
            end = n_indices[i] + 1
            N_coordinates.append((start, end))
            start = None
    
    N_coordinate_df = pd.DataFrame(N_coordinates, columns=['start', 'end'])
    N_coordinate_df[['start', 'end']] += (pid - 1) * window_size
    print(f'Process N Complete: {pid}/{total_task}')

    return N_coordinate_df

def decompose_sequence(task):
    pid, total_task, sequence, ksize, tree = task
    seg = Decompose(sequence, ksize, tree)
    seg.count_kmers()
    seg.build_dbg()
    seg.find_motif()
    print(f'Decomposition Complete: {pid}/{total_task}')
    return ( pid, total_task, seg )

def single_motif_annotation(task):
    pid, total_task, seg = task
    seg.annotate_with_motif()
    df = seg.annotation
    df[['start', 'end']] += (pid - 1) * step_size
    df = df.sort_values(by='end').reset_index(drop=True)
    df.to_csv(f'{args.output}_temp/seg_anno_{pid}.csv', index = False)
    print(f'Single Motif Annotation Complete: {pid}/{total_task}')

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

def calculate_distance(ref, query_list, motif2id):
    df_list = []
    for query in query_list:
        if len(ref) >= len(query):
            tmp_ref, tmp_query = query, ref
        else:
            tmp_ref, tmp_query = ref, query
        
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
        df_list.append([motif2id[tmp_ref], motif2id[tmp_query], min_dist1, False])
        df_list.append([motif2id[tmp_ref], motif2id[tmp_query], min_dist2, True])
    return df_list

def rotate_strings(s):
    n = len(s)
    return [s[i:] + s[:i] for i in range(n)]




##################################
# configuration
##################################
script_path = os.path.realpath(__file__)
script_directory = os.path.dirname(script_path).replace('\\','/')
window_size = 5000
overlap_size = 1000
step_size = window_size - overlap_size
is_test = True

# DEFAULT: python trkscan_test.py testData/testData.fasta testData/testData.annotated.bed
parser = argparse.ArgumentParser(description='TRkScan')
parser.add_argument('input', type = str, help = 'FASTA file you want to annotate')
parser.add_argument('output', type = str, help = 'output prefix')
parser.add_argument('-t', '--thread', type = int, default = 1, help = 'number of threads')
parser.add_argument('-k', '--ksize', type = int, default = 5, help = 'k-mer size for building De Bruijn graph')
parser.add_argument('-m', '--motif', type = str, default = f'{script_directory}/data/refMotif.txt', help='reference motif set')  ###### need to add reference set
parser.add_argument('-n', '--motifnum', type = int, default = 30, help='motif maximum')
parser.add_argument('-f', '--force', action='store_true', help="annotate with motif X in given motif set no matter whether motif X is in the sequence")
args = parser.parse_args()

args.input = args.input.replace('\\','/')
args.output = args.output.replace('\\','/')
os.makedirs(f'{args.output}_temp', exist_ok=True)





if __name__ == "__main__":


    if is_test:
        start_time = time.time()
    ##################################
    # read data
    ##################################
    with open(args.input, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    ##################################
    # read reference motif set
    ##################################
    tree = BKTree(Levenshtein.distance)  # use Levenshtein distance

    with open(args.motif, 'r') as motifDB:
        for line in motifDB:
            tree.add(line.strip())  # add motif
    if is_test:
        end_time = time.time()
        print(f"read data: {end_time - start_time} s")


    annotation_list = []
    for record in records:
        seq_name = record.name
        print(f'Start Processing [{seq_name}]')

        if is_test:
            start_time = time.time()
        ##################################
        # deal with N character          
        ##################################
        seq = str(record.seq)
        seqLenwN = len(seq)
        total_task = max(1, (seqLenwN - window_size) // step_size + 1)
        
        tasks = [(cur + 1, total_task, seq[cur * step_size : (cur + 1) * step_size]) for cur in range(total_task)]
        t = time.time() 
        with Pool(processes = args.thread) as pool:
            results = pool.map(find_N, tasks)
            #pool.close()    # close the pool and don't receive any new tasks
            #pool.join()     # wait for all the tasks to complete
            #results = list(results)
        tt = time.time()
        print(f"slice: {tt - t} s")
        N_coordinate = pd.concat(results, ignore_index=True)
        N_coordinate = N_coordinate.sort_values(by=['start']).reset_index(drop=True)

        ### print(N_coordinate)

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
        
        if is_test:
            end_time = time.time()
            print(f"calculate index transformation: {end_time - start_time} s")

        ##################################
        # split into windows to decompose
        ##################################
        mask = [base != "N" for base in seq]
        filter_seq = ''.join(compress(seq, mask))
        seqLenwoN = len(filter_seq)
        total_task = max(1, (seqLenwoN - overlap_size) // step_size )

        
        if is_test:
            start_time = time.time()
        ##################################
        # get motifs
        ##################################
        with Pool(processes = args.thread) as pool:
            cur = 0
            tasks = []
            while cur < total_task:
                tasks.append(tuple([cur + 1, total_task, filter_seq[cur*step_size : cur*step_size + window_size], args.ksize, tree]))
                cur += 1
            results = pool.imap_unordered(decompose_sequence, tasks)
            
            pool.close()    # close the pool and don't receive any new tasks
            pool.join()     # wait for all the tasks to complete
            results = list(results)

        ##################################
        # merge and decide the representive of motifs
        ##################################
        motifs, motifs_rep, nondup, nondup_rep = [], [], [], []
        for pid, _, result in results:
            motifs.extend(result.motifs['motif'].to_list())
            motifs_rep.extend(result.motifs['value'].to_list())
        for idx in range(len(motifs)):
            is_dup = False
            for idx2, motif in enumerate(nondup):
                if rolling_same(motifs[idx], motif):
                    is_dup = True
                    same_idx= idx2
                    break
            if not is_dup:
                nondup.append(motifs[idx])
                nondup_rep.append(motifs_rep[idx])
            else:
                motifs_rep[same_idx] += motifs_rep[idx]

        print(f'Number of identified motif = {len(nondup)}')

        tmp = pd.DataFrame(columns=['motif','motif_rep'])
        tmp['motif'] = nondup
        tmp['motif_rep'] = nondup_rep
        tmp = tmp.sort_values(by=['motif_rep'], ascending = False).reset_index(drop=True)
        tmp.to_csv(f'{args.output}_motif.txt', sep='\t')
        tmp = tmp.head(args.motifnum)
        nondup = tmp['motif'].to_list()

        for pid, _, result in results:
            result.motifs_list = nondup

        if is_test:
            end_time = time.time()
            print(f"get motif: {end_time - start_time} s")

        if is_test:
            start_time = time.time()
        ##################################
        # annotate single motif
        ##################################
        with Pool(processes = args.thread) as pool:
            cur = 0
            tasks = results
            results = pool.imap_unordered(single_motif_annotation, tasks)
            pool.close()
            pool.join()
                   
        del results
        
        if is_test:
            end_time = time.time()
            print(f"annotate motif: {end_time - start_time} s")

        if is_test:
            start_time = time.time()

        ##################################
        # DP to link
        ##################################
        # parameters
        gap_penalty = 1
        distance_penalty = 1
        perfect_bonus = 0.5

        length = len(filter_seq)
        # Initialize dp and pre arrays
        dp = np.zeros(length + 1, dtype=np.float64)
        pre = np.full((length + 1, 3), None)  # (pre_i, motif_id, motif)

        idx = 0
        df = pd.read_csv(f'{args.output}_temp/seg_anno_1.csv')
        anno_start = df['start'].values
        anno_end = df['end'].values
        anno_motif = df['motif'].values
        anno_distance = df['distance'].values

        # Dynamic programming loop
        for i in range(1, length + 1):
            # read annotation data
            if (i - overlap_size) % step_size == 0 and i > overlap_size:
                pid = (i - overlap_size) // step_size
                idx = 0
                df = pd.read_csv(f'{args.output}_temp/seg_anno_{pid}.csv')
                anno_start = df['start'].values
                anno_end = df['end'].values
                anno_motif = df['motif'].values
                anno_distance = df['distance'].values

            # Skip one base
            if dp[i-1] - gap_penalty > 0:
                dp[i] = dp[i-1] - gap_penalty
                pre[i] = (i-1, None, None)

            # Efficiently iterate over merged_df without repeating checks
            while idx < df.shape[0] and anno_end[idx] <= i:
                if anno_end[idx] < i:
                    idx += 1
                    continue

                pre_i = anno_start[idx]
                motif = anno_motif[idx]
                distance = anno_distance[idx]
                bonus = perfect_bonus * len(motif) if distance == 0 else 0

                score = dp[pre_i] + len(motif) - distance * distance_penalty + bonus
                if score >= dp[i]:
                    dp[i] = score
                    pre[i] = (pre_i, idx, motif)

                idx += 1

            # Print progress every 5000 iterations
            if i % 5000 == 0:
                print(f'DP: {i // 1000} kbp is Done!')

        print('DP complete!')

        # retrace
        idx = length
        next_coords = [None] * (length + 1)  # to store next coordinates
        next_motif = [None] * (length + 1)   # to store motif references

        while idx >= 0:
            if pre[idx][0] is not None:
                if pre[idx][1] is not None:  # matched a motif
                    next_coords[pre[idx][0]] = idx
                    next_motif[pre[idx][0]] = pre[idx][2]
                    idx = pre[idx][0]
                else:  # skipped a base
                    next_coords[idx - 1] = idx
                    next_motif[idx - 1] = None
                    idx -= 1
            else:
                idx -= 1
  
        print('Retrace Complete!')

        if is_test:
            end_time = time.time()
            print(f"DP and retrace: {end_time - start_time} s")

        # output motif annotation file
        annotation_data = []
        idx = 0
        
        while idx < length:
            if pre[idx][0] is None and next_coords[idx] is not None:    # is a start site
                idx2 = idx
                start, end = idx, 0
                cur_motif, rep_num, score, max_score, cigar_string = None, 0, 0, 0, ''
                skip_num = 0
                while next_coords[idx2] != None:
                    motif = next_motif[idx2]
                    if motif is not None:    # match a motif 
                        if cur_motif != motif:       # split the annotation and init
                            if cur_motif != None:
                                annotation_data.append(row)
                            start, cur_motif, rep_num, score, max_score, cigar_string = idx2, motif, 0, 0, 0, ''
                            skip_num = 0

                        if skip_num:
                            cigar_string += f'{skip_num}N'
                            score -= skip_num * gap_penalty
                            skip_num = 0

                        rep_num += 1
                        distance, cigar = get_distacne_and_cigar(motif, filter_seq[idx2: next_coords[idx2]])
                        score += len(cur_motif) - distance * distance_penalty
                        cigar_string += cigar + '/'
                        if score >= max_score:
                            row = [seq_name, seqLenwN, start, next_coords[idx2], motif, rep_num, score, cigar_string]
                            max_score = score
                    else:               # skip a base
                        if idx2 != 0:
                            ### print(idx2)
                            skip_num += 1
                
                    idx2 = next_coords[idx2]
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
    motif_group = merged_annotation.groupby('motif')['rep_num'].sum().reset_index()
    motif_df = motif_group[['motif', 'rep_num']].copy()
    motif_df = motif_df.sort_values(by=['rep_num'], ascending = False).reset_index(drop=True)
    motif_df.index.name = 'id'
    motifs_list = motif_df['motif'].to_list()
    motif_df.to_csv(f'{args.output}.motif.tsv', sep = '\t', index = True)
    ### print([row for row in motif_df.iterrows()])
    motif2id = {row['motif'] : index for index, row in motif_df.iterrows()}
    ### print(motif2id)
    
    ##################################
    # make file - *.dist.tsv
    ##################################
    all_df_row = []
    for motif in motifs_list:
        all_df_row.extend(calculate_distance(motif, motifs_list, motif2id))

    distance_df = pd.DataFrame(all_df_row, columns = ['ref','query','dist','is_rc'])
    distance_df = distance_df[distance_df['ref'] != distance_df['query']]
    distance_df = distance_df.sort_values(by=['dist']).reset_index(drop=True)
    distance_df.to_csv(f'{args.output}.dist.tsv', sep = '\t', index = False)

   





