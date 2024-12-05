
#import xxxxx
import argparse
import pandas as pd
import edlib
import Bio.Seq
from stringDecompose import Decompose
from TRgenerator import TR_multiMotif
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

def decompose_sequence(task):
    pid, total_task, sequence, ksize = task
    seg = Decompose(sequence, ksize)
    seg.count_kmers()
    seg.build_dbg()
    seg.find_motif()
    print(f'work complete: {pid}/{total_task}')
    return ( pid, total_task, seg )

def single_motif_annotation(task):
    pid, total_task, seg = task
    seg.annotate_with_motif()
    print(f'work complete: {pid}/{total_task}')
    return ( pid, total_task, seg )

def rolling_same(seq1, seq2):
    if len(seq1) != len(seq2):
        return False
    print(seq1)
    seq1double = seq1 * 2
    for i in range(len(seq1)):
        ###print(i)
        seq1_new = seq1double[i : i + len(seq1)]
        ###print(seq1_new)
        if seq1_new == seq2:
            return True
    return False

def get_cigar(ref, query):
    matches = edlib.align(query, ref, mode = "NW", task = "path")
    return matches['cigar']

if __name__ == "__main__":
    # configuration
    thread_num = 1
    ksize = 5
    seed = 55889615
    window_size = 5000
    step_size = 4000
    # generate test data
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


    # read data
    with open("testData/testData.fasta", "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    for record in records:
        name = record.name
        ### print(name)
        print(f'start process {name}')
        ### print(record.seq)

        # deal with N character
        # coordicnate transform 

        # split into windows to decompose
        filter_seq = str(record.seq)
        seqLen = len(filter_seq)
        total_task = max(1, int((seqLen - window_size) / step_size + 1) )
        # print(total_task)

        # get motifs
        with Pool(processes = thread_num) as pool:
            cur = 0
            tasks = []
            while cur < total_task:
                ### print(filter_seq[cur*step_size : cur*step_size + window_size])
                tasks.append(tuple([cur + 1, total_task, filter_seq[cur*step_size : cur*step_size + window_size], ksize]))
                cur += 1
            ### print(tasks)
            results = pool.imap_unordered(decompose_sequence, tasks)
            
            pool.close()  # 关闭进程池，不再接受新的任务
            pool.join()
            results = list(results)
            ### print(results)
        
        # merge motifs without reduncy
        motifs, nondup = [], []
        for pid, _, result in results:
            motifs.extend(result.motifs['motif'].to_list())
        for idx in range(len(motifs)):
            is_dup = False
            for motif in nondup:
                ###print('####33')
                if rolling_same(motifs[idx], motif):
                    is_dup = True
                    break
            if not is_dup:
                nondup.append(motifs[idx])

        for pid, _, result in results:
            result.motifs_list = nondup

        print(nondup)
        # decide the representive of motifs
        #######################
        #######################
        #######################

        # annotate single motif
        with Pool(processes = thread_num) as pool:
            cur = 0
            tasks = results
            results = pool.imap_unordered(single_motif_annotation, tasks)
            
            pool.close()  # 关闭进程池，不再接受新的任务
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
        print(merged_df.shape[0])
        print(merged_df.loc[:20,])
        
        # DP to link
        length = len(filter_seq)  ###############
        dp = [0] * (length + 1) # dp[i] means sep[:i] score sum
        pre = [(None, None, None) for i in range(length + 1)] # pre_i, motif_id, motif
        idx = 0
        for i in range(1, length + 1):
            ### print(f'i:{i}')
            if dp[i-1] - 1 >= 0:
                dp[i] = dp[i-1] - 1
                pre[i] = (i-1, None, pre[i-1][2])
            while idx < merged_df.shape[0] and merged_df.loc[idx, 'end'] <= i:
                ### print(f'idx:{idx}')
                pre_i = merged_df.loc[idx, 'start']
                if pre[pre_i][0] != None and pre[pre_i][2] != merged_df.loc[idx, 'motif']: ############## problem!!!!!
                    ### print('dont same motif')
                    idx += 1
                    continue

                if dp[merged_df.loc[idx,'start']] + len(merged_df.loc[idx,'motif']) - merged_df.loc[idx,'distance'] >= dp[i]:
                    ### print(f'{pre_i}->{i}')
                    dp[i] = dp[merged_df.loc[idx,'start']] + len(merged_df.loc[idx,'motif']) - merged_df.loc[idx,'distance'] 
                    pre[i] = (merged_df.loc[idx,'start'], idx, merged_df.loc[idx,'motif'])
                
                idx += 1
        
        ### print(dp[:50])    ###############
        ### print(pre[:50])   ###############

        # retrace
        idx = length
        next = [(None, None) for i in range(length+1)] # next coordinate, reference motif
        while idx >= 0:
            if pre[idx][0] != None:
                if pre[idx][1] != None: # match a motif
                    next[pre[idx][0]] = (idx, merged_df.loc[pre[idx][1],'motif'])
                    idx = pre[idx][0]
                else:  # skip a base
                    next[idx - 1] = (idx, None)
                    idx -= 1
            else:
                idx -= 1

        ### print(next[:40])

        # output motif annotation file
        annotation = pd.DataFrame(columns=['start','end','motif','rep_num','score','CIGAR'])
        '''
        for i in range(length):
            if pre[i][0] is None and next[i][0] != None: # start site
                motif = next[i][1]
                idx, skip_num, rep_num = i, 0, 0
                cigar_string = ''
                while next[idx][0] != None:
                    if next[idx][1]: # motif match
                        if skip_num:
                            cigar_string += f'{skip_num}N'
                            skip_num = 0
                        rep_num += 1
                        cigar_string += get_cigar(next[idx][1], filter_seq[idx : next[idx][0]]) + '/'
                    else: # skip single base
                        skip_num += 1
                    idx = next[idx][0]
                annotation.loc[annotation.shape[0], ] = [i, idx, motif, rep_num, dp[idx], cigar_string]
                print(cigar_string)
        '''
        for i in range(length):
            if pre[i][0] is None and next[i][0] != None: # start site
                motif = next[i][1]
                idx, skip_num, rep_num = i, 0, 0
                cigar_string = ''
                while next[idx][0] != None:
                    if next[idx][1]: # motif match
                        if skip_num:
                            cigar_string += f'{skip_num}N'
                            skip_num = 0
                        rep_num += 1
                        cigar_string += get_cigar(next[idx][1], filter_seq[idx : next[idx][0]]) + '/'
                    else: # skip single base
                        skip_num += 1
                    idx = next[idx][0]
                annotation.loc[annotation.shape[0], ] = [i, idx, motif, rep_num, dp[idx], cigar_string]
                print(cigar_string[:50])


        print(annotation)
            


        # readd N character




