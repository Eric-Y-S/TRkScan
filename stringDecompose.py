import sourmash
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict
import pydivsufsort
import edlib
from TRgenerator import TR_singleMotif, TR_multiMotif

IS_PAINT = False    # if you want to plot dbg, set True; visualization is unless for complex TR region

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

    motif_match_list = []
    motif_positions = list()

    for motif_id, motif in enumerate(motifs):
        max_distance = max_distances[motif_id]
        motif_len = len(motif)

        idx = 0
        while idx < len(sa):
            s = seq[sa[idx]:]
            matches = edlib.align(motif, s, mode="SHW", task="locations", k = max_distance)

            if matches['editDistance'] == -1:   # not match
                idx2 = idx
                while idx2 < len(lcp) and lcp[idx2] >= motif_len + max_distance:
                    idx2 += 1
                idx = idx2 + 1
                continue

            for start, end in matches['locations']:
                motif_positions.append(
                    (start + sa[idx], end + 1 + sa[idx], motif, s[start:end + 1], matches['editDistance'])
                )
            
            max_match_len = max(end + 1 for _, end in matches['locations'])

            idx2 = idx
            while idx2 < len(lcp) and lcp[idx2] >= max_match_len:
                for start, end in matches['locations']:
                    motif_positions.append(
                        (start + sa[idx2 + 1], end + 1 + sa[idx2 + 1], motif, s[start:end + 1], matches['editDistance'])
                    )
                idx2 += 1

            idx = idx2 + 1

    motif_match_df = pd.DataFrame(motif_positions, columns=['start', 'end', 'motif', 'seq', 'distance'])
    motif_match_df = motif_match_df.sort_values(by=['start']).reset_index(drop=True)

    return motif_match_df

class Decompose:
    def __init__(self, sequence, ksize):
        self.sequecne = sequence
        self.ksize = ksize
        self.abud_treshold = 0.01
        self.abud_min = 3
        self.dist_ratio = 0.4
        self.kmer = None
        self.dbg = None
        self.motifs = None
        self.motifs_list = None
        self.annotation = None

    def count_kmers(self):
        if self.kmer is None:
            mh = sourmash.MinHash(n = 5000, ksize = self.ksize + 1)
            kmers = [kmer for kmer, _ in mh.kmers_and_hashes(self.sequecne, force=True) ]
            kmer_count = Counter(kmers)
            max_count = kmer_count.most_common(1)[0][1]
            min_count = max(self.abud_min, max_count * self.abud_treshold)
            filtered_kmer_count = {k: v for k, v in kmer_count.items() if v >= min_count}
            self.kmer = filtered_kmer_count
        return self.kmer

    def build_dbg(self):
        if self.dbg is None:
            dbg = nx.DiGraph()
            for k in self.kmer:
                dbg.add_edge(k[:-1], k[1:], weight = self.kmer[k])

            if IS_PAINT:
                pos = calculate_circular_layout(dbg)
                paint_dbg(dbg, pos, 'initial_graph.pdf')
            
            # make compact De Bruijn graph
            nodes_to_merge = [n for n in dbg.nodes() if dbg.in_degree(n) == 1 and dbg.out_degree(n) == 1] 

            cot = 1
            for node in nodes_to_merge:
                if node not in dbg.nodes():
                    continue

                # find the start
                added_nodes = [node]
                pre = list(dbg.predecessors(node))[0]
                while dbg.in_degree(pre) == 1 and dbg.out_degree(pre) == 1 and not pre in added_nodes:
                    added_nodes.append(pre)
                    pre = list(dbg.predecessors(pre))[0]
                
                added_nodes = added_nodes[::-1]

                # find the end
                next = list(dbg.successors(node))[0]
                while dbg.in_degree(next) == 1 and dbg.out_degree(next) == 1 and not next in added_nodes:
                    added_nodes.append(next)
                    next = list(dbg.successors(next))[0]

                if len(added_nodes) == 1:    # Simple loop; no merging
                    continue
                
                merged_node = added_nodes[0] + ''.join(node[-1] for node in added_nodes[1:])
                start = added_nodes[0]
                end = added_nodes[-1]
                predecessors = list(dbg.predecessors(start))
                successors = list(dbg.successors(end))

                if dbg.has_edge(end, start): # just a simple loop
                    dbg.add_edge(merged_node, merged_node, weight = dbg[end][start]['weight'])
                else:
                    for predecessor in predecessors:
                        dbg.add_edge(predecessor, merged_node, weight = dbg[predecessor][start]['weight'])
                    for successor in successors:
                        dbg.add_edge(merged_node, successor, weight = dbg[end][successor]['weight'])
                
                # delete original node
                dbg.remove_nodes_from(added_nodes)

                if(IS_PAINT):
                    pos[merged_node] = pos[start]
                    for n in added_nodes:
                        pos.pop(n, None)
                    paint_dbg(dbg, pos, f'compacted-step{cot}.pdf')

                cot += 1

            self.dbg = dbg

        return self.dbg

    def find_motif(self):
        if self.motifs is None:
            cycles = nx.simple_cycles(self.dbg)
            motifs = []
            for cycle in cycles:
                motif = ''.join([node[(self.ksize - 1):] for node in cycle])
                cot = [self.dbg[cycle[i]][cycle[i+1]]['weight'] for i in range(len(cycle) - 1)]
                cot.append(self.dbg[cycle[-1]][cycle[0]]['weight'])  # add the weight of the last edge

                motifs.append([cycle, motif, min(cot)])
            motifs_df = pd.DataFrame(motifs, columns=['cycle', 'motif', 'value'])
            motifs_df = motifs_df.sort_values(by=['value'], ascending = False).reset_index(drop=True)
            self.motifs = motifs_df
            self.motifs_list = motifs_df['motif'].to_list()
        return self.motifs

    def annotate_with_motif(self):
        if self.annotation is None:
            motifs = self.motifs_list

            max_distances = [ int( len(motif) * self.dist_ratio ) for motif in motifs ]
            motif_match_df = find_similar_match(self.sequecne, motifs, max_distances)
    
            self.annotation = motif_match_df
        
        return self.annotation       


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


