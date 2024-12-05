import sourmash
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from pysuffixarray.core import SuffixArray
import pydivsufsort
import edlib
from TRgenerator import TR_singleMotif, TR_multiMotif

IS_PAINT = True

# calculate layout for graph
def calculate_circular_layout(g):
    offset = 0
    step = 200
    radius = 300
    pos = {}
    num_nodes = len(g.nodes())

    ori_pos = nx.circular_layout(g, scale = 500)

    cycles = list(nx.simple_cycles(g))
    for cycle in cycles:
        angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint = False)  # 计算每个节点的角度
        for idx in range(len(cycle)):
            if (cycle[idx] not in pos.keys()):
                pos[cycle[idx]] = np.array([offset + radius*np.cos(angles[idx]), offset + radius*np.sin(angles[idx])])
        offset += step
    for node in pos.keys():
        ori_pos[node] = pos[node]
    return ori_pos

# paint
def paint_dbg(g, pos, filename):
    plt.figure(figsize=(5, 5))
    #edge_labels = nx.get_edge_attributes(g, 'weight')
    nx.draw_networkx_nodes(g, pos, node_color = '#A7C957', node_size = 300)
    nx.draw_networkx_edges(
        g, pos, width = 1, 
        connectionstyle="arc3,rad=0.2"  # curved edges
    )
    nx.draw_networkx_labels(g, pos, font_size=5)

    # hide frame
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.savefig(filename,dpi=300)

    return 

def find_similar_match(seq:str, motifs:list, max_distances:list):
    seq_bytes = seq.encode('utf-8')
    sa = pydivsufsort.divsufsort(seq_bytes)
    lcp = pydivsufsort.kasai(seq_bytes, sa)

    ### print(sa)
    ### print(lcp)

    motif_match_df = pd.DataFrame(columns=['start', 'end', 'motif', 'seq', 'distance'])
    for motif_id in range(len(motifs)):
        motif = motifs[motif_id]
        ### print(motif) #################
        idx = 0
        while idx < len(sa):
            ### print(idx)  #################
            s = seq[sa[idx]:]
            matches = edlib.align(motif, s, mode = "SHW", task = "path", k = max_distances[motif_id])
            # do not match
            if matches['editDistance'] == -1:
                idx += 1
                continue
            # have matches
            l = 0
            for start, end in matches['locations']:
                l = max(l, end+1)
                motif_match_df.loc[motif_match_df.shape[0]] = [start + sa[idx], end + 1 + sa[idx], motif, s[start: end + 1],matches['editDistance']]
            ### print(matches)  #################
            for idx2 in range(idx, len(lcp)):
                if lcp[idx2] >= l:
                    for start, end in matches['locations']:
                        motif_match_df.loc[motif_match_df.shape[0]] = [start + sa[idx2+1], end + 1 + sa[idx2+1], motif, s[start: end + 1],matches['editDistance']]
                else:
                    break
            idx = max(idx, idx2) + 1
            
    motif_match_df = motif_match_df.sort_values(by=['start'])
    motif_match_df = motif_match_df.reset_index(drop=True)

    return motif_match_df

class Decompose:
    def __init__(self, sequence, ksize):
        self.sequecne = sequence
        self.ksize = ksize
        self.abud_treshold = 0.1
        self.dist_ratio = 0.4
        self.kmer = None
        self.dbg = None
        self.motifs = None
        self.motifs_list = None
        self.annotation = None

    def count_kmers(self):
        if self.kmer is None:
            mh = sourmash.MinHash(n = 1000, ksize = self.ksize + 1)
            kmers = [kmer for kmer, _ in mh.kmers_and_hashes(self.sequecne, force=True) ]
            kmer_count = Counter(kmers)
            max_count = kmer_count.most_common(1)[0][1]
            min_count = max_count * self.abud_treshold
            filtered_kmer_count = {k: v for k, v in kmer_count.items() if v >= min_count}
            self.kmer = filtered_kmer_count
        return self.kmer

    def build_dbg(self):
        if self.dbg is None:
            dbg = nx.DiGraph()
            for k in self.kmer:
                dbg.add_edge(k[:-1], k[1:], weight = self.kmer[k])

            if(IS_PAINT):
                pos = calculate_circular_layout(dbg)
                ### print(pos)
                paint_dbg(dbg, pos, 'test.pdf')
            
            #
            # make compact De Bruijn graph
            #
            nodes_to_merge = [n for n in dbg.nodes() if dbg.in_degree(n) == 1 and dbg.out_degree(n) == 1] 

            cot = 1
            ### print(nodes_to_merge)
            for node in nodes_to_merge:
                if not node in dbg.nodes():
                    continue

                start = node
                added_nodes = [node]

                pre = list(dbg.predecessors(start))[0]
                while dbg.in_degree(pre) == 1 and dbg.out_degree(pre) == 1 and not pre in added_nodes:
                    added_nodes.append(pre)
                    start = pre
                    pre = list(dbg.predecessors(start))[0]
                end = node
                next = list(dbg.successors(end))[0]
                while dbg.in_degree(next) == 1 and dbg.out_degree(next) == 1 and not next in added_nodes:
                    added_nodes.append(next)
                    end = next
                    next = list(dbg.successors(end))[0]

                if start == end:
                    continue
                
                predecessors = list(dbg.predecessors(start))
                successors = list(dbg.successors(end))

                nodes_to_del = [start]
                cur = merged_node = start
                while list(dbg.successors(cur))[0] != end:
                    cur = list(dbg.successors(cur))[0]
                    merged_node += cur[-1]
                    nodes_to_del.append(cur)
                merged_node += end[-1]
                nodes_to_del.append(end)

                ### print(nodes_to_del) 

                if dbg.has_edge(end, start): # just a simple loop
                    dbg.add_edge(merged_node, merged_node, weight = dbg[end][start]['weight'])
                else:
                    for predecessor in predecessors:
                        dbg.add_edge(predecessor, merged_node, weight = dbg[predecessor][start]['weight'])
                    for successor in successors:
                        dbg.add_edge(merged_node, successor, weight = dbg[end][successor]['weight'])
                
                # delete original node
                dbg.remove_nodes_from(nodes_to_del)

                if(IS_PAINT):
                    pos[merged_node] = pos[start]
                    for n in nodes_to_del:
                        pos.pop(n, None)
                    ### print(pos)
                    paint_dbg(dbg, pos, f'test-step{cot}.pdf')

                cot += 1

            self.dbg = dbg

        return self.dbg

    def find_motif(self):
        if self.motifs is None:
            cycles = nx.simple_cycles(self.dbg)
            motifs = pd.DataFrame(columns=['cycle', 'motif', 'value'])
            for cycle in cycles:
                motif = ''
                for node in cycle:
                    motif += node[(self.ksize - 1) : ]
                # calculate counts
                cot = []
                cycle.append(cycle[0])
                for i in range(len(cycle)-1):
                    cot.append(self.dbg[cycle[i]][cycle[i+1]]['weight'])

                motifs.loc[motifs.shape[0]] = [cycle, motif, min(cot)]
            motifs = motifs.sort_values(by=['value'], ascending = False)
            motifs = motifs.reset_index(drop=True)
            self.motifs = motifs
            self.motifs_list = motifs['motif'].to_list()
            ### print(motifs)
        return self.motifs

    def annotate_with_motif(self):
        if self.annotation is None:
            motifs = self.motifs_list
            print(motifs)  ######################################

            max_distances = [ int( len(motif) * self.dist_ratio ) for motif in motifs]
            motif_match_df = find_similar_match(self.sequecne, motifs, max_distances)
    
            self.annotation = motif_match_df
        
        return self.annotation
        '''for motif in motifs:
            matches = sa.match(motif)
            l = len(motif)
            s = matches[0]
            anno_tmp = pd.DataFrame(columns=['start', 'end', 'motif','CIGAR-like'])
            for idx in range(1,len(matches)):
                if matches[idx] != matches[idx - 1] + l:
                    e = matches[idx - 1] + l
                    rep = int((e-s)/l)
                    anno_tmp.loc[anno_tmp.shape[0]] = [s, e, motif, f'{rep}P']
                    s = matches[idx]
            annotation_df = pd.concat([annotation_df, anno_tmp], ignore_index=True)'''

        

        # link adjacent region with the same motif
        '''
        linked_df = pd.DataFrame(columns=['start', 'end', 'motif','CIGAR-like'])
        s = annotation_df.loc[0,'start']
        string = ''
        for i in range(annotation_df.shape[0]):
            if annotation_df.loc[i,'motif'] != annotation_df.loc[i+1,'motif']:
                e = annotation_df.loc[i,'end']
                linked_df.loc[linked_df.shape[0]] = [s, e, motif, f'{rep}P']
                s = annotation_df.loc[i+1, 'start']
            else:
                string += annotation_df.loc[i, 'CIGAR-like']
                string += xxxx()   ################## solve with inter region zone
        
        '''



        # polish border




        # solve with intersect




        # xxxxxxxxxxxxx




        


if __name__ == "__main__":
    ##############################
    # config
    ##############################
    k = 5
    seed = 55889615
    #tr = TR_singleMotif('AATGG', 100, 0.1, seed)
    #sequence = tr.sequence
    #print(sequence)
    
    ##############################
    # count k-mer
    ##############################
    hsat2_hsat3 = TR_multiMotif(['AATGG','CCATT','CTT'], 5000, 0.1, seed)
    print(hsat2_hsat3.sequence[:100])
    example = Decompose(hsat2_hsat3.sequence, k)
    example.count_kmers()
    example.build_dbg()
    example.find_motif()
    example.annotate_with_motif()
    print(example.annotate_with_motif())
    #print(hsat2_hsat3.annotation)
    #example.count_kmers()


