import argparse
import pandas as pd
import edlib
from itertools import compress, accumulate
from stringDecompose import Decompose
from TRgenerator import TR_multiMotif
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

def find_N(task):
    pid, total_task, sequence = task
    N_coordinate = pd.DataFrame(columns=['start','end'])
    start = None
    for i in range(len(sequence)):
        if sequence[i] == 'N':
            if start is None:
                start = i
        else:
            if start is not None:
                N_coordinate.loc[N_coordinate.shape[0],] = [start, i]
                start = None
    if start is not None:
        N_coordinate.loc[N_coordinate.shape[0],] = [start, len(sequence) - 1]
    print(f'Process N Complete: {pid}/{total_task}')
    return ( pid, total_task, N_coordinate )

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
    for i in range(len(seq1)):
        seq1_new = seq1double[i : i + len(seq1)]
        if seq1_new == seq2:
            return True
    return False

def get_cigar(ref, query):
    matches = edlib.align(query, ref, mode = "NW", task = "path")
    return matches['cigar']

def get_distacne(ref, query):
    matches = edlib.align(query, ref, mode = "NW", task = "path")
    return matches['editDistance']

if __name__ == "__main__":
    ##################################
    # configuration
    ##################################
    # DEFAULT: python trkscan_test.py testData/testData.fasta testData/testData.annotated.bed
    parser = argparse.ArgumentParser(description='TRkScan')
    parser.add_argument('input', type = str, help = 'FASTA file you want to annotate')
    parser.add_argument('output', type = str, help = 'output prefix')
    parser.add_argument('-t', '--thread', type = int, default = 1, help = 'number of threads')
    parser.add_argument('-k', '--ksize', type = int, default = 5, help = 'k-mer size for building De Bruijn graph')
    parser.add_argument('-m', '--motif', type = str, default = '', help='reference motif set')  ###### need to add reference set
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

    for record in records:
        seq_name = record.name
        print(f'Start Processing [{seq_name}]')

        ##################################
        # deal with N character          
        ##################################
        seq = str(record.seq)
        seqLen = len(seq)
        total_task = max(1, int((seqLen - window_size) / step_size + 1) )
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
        
        N_coordinate = pd.DataFrame(columns=['start', 'end'])
        for pid, _, result in results:
            result.loc[:,'start'] += (pid - 1) * window_size
            result.loc[:,'end'] += (pid - 1) * window_size
            N_coordinate = pd.concat([N_coordinate, result], ignore_index = True)

        N_coordinate = N_coordinate.sort_values(by=['start'])
        N_coordinate = N_coordinate.reset_index(drop=True)

        if N_coordinate.shape[0]:
            N_coordinate_new = pd.DataFrame(columns=['start', 'end'])
            start = N_coordinate.loc[0, 'start']
            for i in range(N_coordinate.shape[0]-1):
                if N_coordinate.loc[i, 'end'] != N_coordinate.loc[i+1, 'start']:
                    N_coordinate_new.loc[N_coordinate_new.shape[0], ] = [start, N_coordinate.loc[i, 'end']]
                    start = N_coordinate.loc[i+1, 'start']

            N_coordinate_new.loc[N_coordinate_new.shape[0], ] = [start, N_coordinate.loc[N_coordinate.shape[0] - 1, 'end']]
            N_coordinate_new['old_index'] = N_coordinate_new['end']
            N_coordinate_new['new_index'] = N_coordinate_new['end']
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
        seqLen = len(filter_seq)
        total_task = max(1, int((seqLen - window_size) / step_size + 1) )

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

        for pid, _, result in results:
            result.motifs_list = nondup

        print(nondup)


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
        merged_df = pd.DataFrame(columns=['start', 'end', 'motif', 'seq', 'distance'])
        for pid, _, result in results:
            df = result.annotation
            df.loc[:,'start'] += (pid - 1) * step_size
            df.loc[:,'end'] += (pid - 1) * step_size
            merged_df = pd.concat([merged_df, df], ignore_index = True)

        merged_df = merged_df.drop_duplicates()
        merged_df = merged_df.sort_values(by=['end'])
        merged_df = merged_df.reset_index(drop=True)
        ### print(merged_df.shape[0])
        ### print(merged_df.loc[:20,])
        
        ##################################
        # DP to link
        ##################################
        # parameters
        gap_penalty = 1
        distance_penalty = 1
        perfect_bonus = 0.5    ######## need to change !!!!!!!!!!

        length = len(filter_seq)
        dp = [0] * (length + 1)     # dp[i] means sep[:i] score sum
        pre = [(None, None, None) for i in range(length + 1)]   # (pre_i, motif_id, motif)
        idx = 0
        for i in range(1, length + 1):
            # skip one base
            if dp[i-1] - 1 > 0:
                dp[i] = dp[i-1] - 1
                pre[i] = (i-1, None, None)
            
            while idx < merged_df.shape[0] and merged_df.loc[idx, 'end'] <= i:
                if merged_df.loc[idx, 'end'] < i:
                    idx += 1
                    continue
                
                pre_i = merged_df.loc[idx, 'start']
                motif = merged_df.loc[idx,'motif']
                distance = merged_df.loc[idx,'distance']
                bonus = perfect_bonus * len(motif) if distance == 0 else 0
                
                if dp[merged_df.loc[idx,'start']] + len(motif) - distance * distance_penalty + bonus >= dp[i]:
                    dp[i] = dp[pre_i] + len(motif) - distance * distance_penalty + bonus
                    pre[i] = (pre_i, idx, motif)
                
                if idx % 5000 == 0:
                    print(f'DP {idx/1000}kbp is Done!')

                idx += 1
        
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
        annotation = pd.DataFrame(columns=['seq','start','end','motif','rep_num','score','CIGAR'])
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
                                annotation.loc[annotation.shape[0], ] = row
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
                        score += len(cur_motif) - get_distacne(next[idx2][1], filter_seq[idx2 : next[idx2][0]]) * distance_penalty
                        cigar_string += get_cigar(next[idx2][1], filter_seq[idx2 : next[idx2][0]]) + '/'
                        if score >= max_score:
                            row = [seq_name, start, next[idx2][0], cur_motif, rep_num, score, cigar_string]
                            max_score = score
                    else:               # skip a base
                        if idx2 != 0:
                            ### print(idx2)
                            skip_num += 1
                
                    idx2 = next[idx2][0]
                    end = idx2
                annotation.loc[annotation.shape[0], ] = row

            idx += 1

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
        ### print(annotation.loc[0,'CIGAR'])
        annotation.to_csv(f'{args.output}.concise.tsv', sep = '\t', index = False)
