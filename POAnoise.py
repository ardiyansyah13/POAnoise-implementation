import networkx as nx
import numpy as np
import pandas as pd
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from sklearn.metrics.pairwise import pairwise_distances
from scipy.cluster.hierarchy import linkage, fcluster
import time

# -----------------------------
# Input configuration
# -----------------------------
inputs = [
    {
        "preprocessed_file": "the path to the preprocessed noisy sequence dataset (FASTQ file)",
        "abundance_tsv": "the path to the abundance of teh noisy sequence dataset",
        "origin_map": "the path to the origin map of the preprocessed noisy sequence dataset",
        "tag": "the tag of the dataset (e.g. A, B, or C)"
    }
]


# -----------------------------
# Helper functions
# -----------------------------
def parse_sequences(file_path):
    sequences = []
    first_char = open(file_path).read(1)
    if first_char == '>':
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq))
    elif first_char == '@':
        for record in SeqIO.parse(file_path, "fastq"):
            sequences.append(str(record.seq))
    else:
        raise ValueError("File format not recognized (FASTA/FASTQ).")
    return sequences

def load_abundance_data(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t')
    if 'Sequence' not in df.columns or 'Abundance' not in df.columns:
        raise ValueError("TSV must contain 'Sequence' and 'Abundance'")
    abundance_dict = dict(zip(df['Sequence'], df['Abundance']))
    total_abundance = df['Abundance'].sum()
    return abundance_dict, total_abundance

def load_origin_map(tsv_file):
    """
    Load the origin map TSV file.
    Maps filtered sequences to their true origin sequences.
    """
    df = pd.read_csv(tsv_file, sep='\t')
    origin_map = defaultdict(list)
    for _, row in df.iterrows():
        origin_map[row['Filtered_Read_ID']].append(row['True_Sequence_ID'])
    return origin_map

def score(a, b, match=1, mismatch=-1):
    return match if a == b else mismatch

def dual_affine_gap_score(n, open1, ext1, open2, ext2):
    return max(open1 + n * ext1, open2 + n * ext2)

def initialize_dp(rows, cols):
    shape = (rows + 1, cols + 1)
    dp = {s: np.full(shape, -np.inf) for s in ['M','I1','I2','D1','D2']}
    tb = {s: np.full(shape, None, dtype=object) for s in ['M','I1','I2','D1','D2']}
    return dp, tb

# -----------------------------
# POA alignment functions
# -----------------------------
def align(poa, sequence, match=1, mismatch=-1, open1=-1, ext1=-1, open2=-5, ext2=0):
    sorted_nodes = list(nx.topological_sort(poa))
    node_index = {node: j+1 for j, node in enumerate(sorted_nodes)}  # 1-based
    seq_len = len(sequence)
    dp, tb = initialize_dp(seq_len, len(sorted_nodes))
    dp['M'][0][1] = 0

    for i in range(1, seq_len+1):
        for j, node in enumerate(sorted_nodes,1):
            base_j = poa.nodes[node]['base']
            # MATCH
            best_M = -np.inf
            for pred in poa.predecessors(node):
                pred_j = node_index[pred]
                for s in ['M','I1','I2','D1','D2']:
                    cand = dp[s][i-1][pred_j] + score(sequence[i-1], base_j, match, mismatch)
                    if cand > best_M:
                        best_M = cand
                        tb['M'][i][j] = (i-1, pred_j, s)
            dp['M'][i][j] = best_M

            # INSERTIONS
            if i>=1:
                cand_I1 = dp['M'][i-1][j]+open1
                if cand_I1>dp['I1'][i][j]:
                    dp['I1'][i][j]=cand_I1; tb['I1'][i][j]=(i-1,j,'M')
                cand_I1_ext = dp['I1'][i-1][j]+ext1
                if cand_I1_ext>dp['I1'][i][j]:
                    dp['I1'][i][j]=cand_I1_ext; tb['I1'][i][j]=(i-1,j,'I1')
                cand_I2 = dp['M'][i-1][j]+open2
                if cand_I2>dp['I2'][i][j]:
                    dp['I2'][i][j]=cand_I2; tb['I2'][i][j]=(i-1,j,'M')
                cand_I2_ext = dp['I2'][i-1][j]+ext2
                if cand_I2_ext>dp['I2'][i][j]:
                    dp['I2'][i][j]=cand_I2_ext; tb['I2'][i][j]=(i-1,j,'I2')

            # DELETIONS
            for pred in poa.predecessors(node):
                pred_j = node_index[pred]
                # D1
                cand_D1 = dp['M'][i][pred_j]+open1
                if cand_D1>dp['D1'][i][j]:
                    dp['D1'][i][j]=cand_D1; tb['D1'][i][j]=(i,pred_j,'M')
                cand_D1_ext = dp['D1'][i][pred_j]+ext1
                if cand_D1_ext>dp['D1'][i][j]:
                    dp['D1'][i][j]=cand_D1_ext; tb['D1'][i][j]=(i,pred_j,'D1')
                # D2
                cand_D2 = dp['M'][i][pred_j]+open2
                if cand_D2>dp['D2'][i][j]:
                    dp['D2'][i][j]=cand_D2; tb['D2'][i][j]=(i,pred_j,'M')
                cand_D2_ext = dp['D2'][i][pred_j]+ext2
                if cand_D2_ext>dp['D2'][i][j]:
                    dp['D2'][i][j]=cand_D2_ext; tb['D2'][i][j]=(i,pred_j,'D2')

    # Best score
    best_score=-np.inf; best_pos=None
    for j in range(1,len(sorted_nodes)+1):
        for s in ['M','I1','I2','D1','D2']:
            if dp[s][seq_len][j]>best_score:
                best_score=dp[s][seq_len][j]
                best_pos=(seq_len,j,s)
    return dp, tb, best_pos, node_index, sorted_nodes

def traceback(poa, tb, best_pos, node_index, sequence):
    i,j,state = best_pos
    aligned_seq, aligned_poa, aligned_nodes = [],[],[]
    reverse_index = {v:k for k,v in node_index.items()}
    while tb[state][i][j] is not None:
        prev_i, prev_j, prev_state = tb[state][i][j]
        node = reverse_index[j]
        base_seq = sequence[i-1] if i>0 else None
        if state=='M':
            aligned_seq.append(base_seq)
            aligned_poa.append(poa.nodes[node]['base'])
            aligned_nodes.append(node)
        elif state in ['I1','I2']:
            aligned_seq.append(base_seq)
            aligned_poa.append('-')
            aligned_nodes.append(None)
        elif state in ['D1','D2']:
            aligned_seq.append('-')
            aligned_poa.append(poa.nodes[node]['base'])
            aligned_nodes.append(node)
        i,j,state = prev_i, prev_j, prev_state
    return list(reversed(aligned_seq)), list(reversed(aligned_poa)), list(reversed(aligned_nodes))

def update_poa(poa, sequence, match=1, mismatch=-1, open1=-1, ext1=-1, open2=-5, ext2=0):
    dp,tb,best_pos,node_index,sorted_nodes=align(poa,sequence,match,mismatch,open1,ext1,open2,ext2)
    aln_seq, aln_poa, aln_nodes = traceback(poa,tb,best_pos,node_index,sorted_nodes,sequence)
    prev_node = None
    seq_id = sequence
    next_pos = max([node[0] for node in poa.nodes]+[-1])+1
    for base_seq, base_poa, node in zip(aln_seq,aln_poa,aln_nodes):
        if node is not None:
            if 'reads' not in poa.nodes[node]:
                poa.nodes[node]['reads']=set()
            poa.nodes[node]['reads'].add(seq_id)
            current_node = node
        else:
            node_id=(next_pos, base_seq)
            next_pos+=1
            poa.add_node(node_id, base=base_seq, position=node_id[0], reads=set([seq_id]))
            current_node=node_id
        if prev_node is not None:
            if not poa.has_edge(prev_node,current_node):
                poa.add_edge(prev_node,current_node)
        prev_node=current_node
    return poa

def build_partial_order_graph_dual_affine(sequences, match=1, mismatch=-1, open1=-1, ext1=-1, open2=-5, ext2=0):
    if not sequences: raise ValueError("No sequences provided")
    poa=nx.DiGraph()
    reference_sequence=sequences[0]; reference_id=reference_sequence
    prev_nodes=[]
    for i,base in enumerate(reference_sequence):
        node_id=(i,base)
        poa.add_node(node_id, base=base, position=i, reads=set([reference_id]))
        if prev_nodes:
            poa.add_edge(prev_nodes[-1], node_id)
        prev_nodes.append(node_id)
    for seq in sequences[1:]:
        poa=update_poa(poa, seq, match, mismatch, open1, ext1, open2, ext2)
    return poa

# -----------------------------
# Edge matrix & consensus
# -----------------------------
def calculate_edge_weights_optimized(graph, sequences):
    edge_weights = defaultdict(float)
    for u,v in graph.edges():
        for seq in sequences:
            if u[0]<len(seq) and v[0]<len(seq):
                if seq[u[0]]==u[1] and seq[v[0]]==v[1]:
                    edge_weights[(u,v)]+=1
    return edge_weights

def heaviest_bundle(graph, edge_weights):
    node_scores=defaultdict(float); traceback_map={}
    for node in nx.topological_sort(graph):
        max_score,best_pred=0,None
        for pred in graph.predecessors(node):
            edge_score=edge_weights.get((pred,node),0)
            total_score=node_scores[pred]+edge_score
            if total_score>max_score:
                max_score=total_score; best_pred=pred
        node_scores[node]=max_score
        if best_pred is not None: traceback_map[node]=best_pred
    max_score=max(node_scores.values(),default=0)
    end_nodes=[node for node,score in node_scores.items() if score==max_score]
    all_paths=[]
    for end_node in end_nodes:
        path=[]; current=end_node
        while current in traceback_map:
            path.append(current); current=traceback_map[current]
        path.append(current)
        path.reverse(); all_paths.append(path)
    return all_paths

def generate_consensus_using_heaviest_bundle(pog, edge_weights):
    heaviest_paths = heaviest_bundle(pog, edge_weights)
    consensus_sequences = []
    for path in heaviest_paths:
        cons = ''.join([pog.nodes[n]['base'] for n in path])
        weight = sum([edge_weights.get((path[i], path[i+1]), 0) for i in range(len(path)-1)])
        contributing_reads = set()
        for node in path:
            contributing_reads.update(pog.nodes[node].get('reads', set()))
        consensus_sequences.append((cons, weight, contributing_reads))
    return consensus_sequences

def get_edge_sequence_matrix(poa, sequences, abundance_dict):
    edge_list=list(poa.edges()); edge_to_index={e:i for i,e in enumerate(edge_list)}
    seq_to_index={seq:i for i,seq in enumerate(sequences)}
    M=np.zeros((len(sequences), len(edge_list)), dtype=float)
    for e in edge_list:
        src,dst=e
        reads_src=poa.nodes[src].get("reads",set())
        reads_dst=poa.nodes[dst].get("reads",set())
        sequences_on_edge=reads_src & reads_dst
        for seq in sequences_on_edge:
            i=seq_to_index[seq]; j=edge_to_index[e]
            M[i,j]=abundance_dict.get(seq,1)
    return M, edge_list

def cluster_sequences_dynamic(distance_matrix, cluster_threshold, method):
    linkage_matrix = linkage(distance_matrix, method)
    cluster_ids = fcluster(linkage_matrix, cluster_threshold, criterion='distance')
    return cluster_ids

def build_consensus_from_lists(sequences, cluster_ids, abundance_dict):
    clusters = defaultdict(list)
    for seq, cid in zip(sequences, cluster_ids):
        clusters[cid].append(seq)

    consensus_records = []
    for cid, seq_list in clusters.items():
        if not seq_list:
            continue
        # build POA for this cluster
        poa = build_partial_order_graph_dual_affine(seq_list)
        # calculate edge weights
        weights = calculate_edge_weights_optimized(poa, seq_list, abundance_dict)
        # generate consensus sequences with contributing reads
        cons_list = generate_consensus_using_heaviest_bundle(poa, weights)
        # add to final list
        consensus_records.extend(cons_list)

    # deduplicate consensus sequences (merge contributing reads for same sequence)
    return consensus_records

def save_consensus_to_file(consensus_list, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file,'w') as f:
        for i,(seq,weight) in enumerate(consensus_list,1):
            f.write(f">Consensus_{i}_Weight_{weight}\n{seq}\n")
    print(f"Consensus sequences saved to {output_file}")

def save_consensus_table(consensus_list, abundance_dict, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file,'w', newline='') as f:
        writer=csv.writer(f, delimiter='\t')
        writer.writerow(["Consensus","Abundance"])
        for seq,weight in consensus_list:
            ab=sum(v for k,v in abundance_dict.items() if k==seq)
            writer.writerow([seq, ab if ab>0 else weight])
    print(f"Consensus table saved to {output_file}")

def save_consensus_with_origin(consensus_list, seq_to_origin, fasta_file):
    os.makedirs(os.path.dirname(fasta_file), exist_ok=True)

    with open(fasta_file, 'w') as f_fasta:
        for i, (seq, weight, contributing_reads) in enumerate(consensus_list, 1):
            origin_set = set()
            for read in contributing_reads:
                origin_set.update(seq_to_origin.get(read, []))
            origin = ",".join(origin_set) if origin_set else "NEW"
            
            # FASTA header
            f_fasta.write(f">Consensus_{i}_Weight_{weight}_Origin_{origin}\n{seq}\n")
            

    print(f"Consensus sequences with origin saved to {fasta_file}")


# -----------------------------
# Main workflow
# -----------------------------

cluster_thresholds = [0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56,  0.61,0.66,  0.71, 0.76, 0.81, 0.86, 0.91, 0.96]
linkage_methods = [ "single"]

for entry in inputs:
    tag = entry["tag"]
    cluster_file = entry["cluster_file"]
    abundance_tsv = entry["abundance_tsv"]
    origin_map_file = entry["origin_map"]

    print(f"\nProcessing input {tag}...")

    sequences = parse_sequences(cluster_file)
    abundance_dict, total_abundance = load_abundance_data(abundance_tsv)
    origin_map = load_origin_map(origin_map_file)

    # Build POAs
    poa = build_partial_order_graph_dual_affine(sequences)
    M, edge_list = get_edge_sequence_matrix(poa, sequences, abundance_dict)

    distance_matrix = pairwise_distances(M, metric='cosine')
    normalized_matrix = (distance_matrix - distance_matrix.min()) / (
        distance_matrix.max() - distance_matrix.min()
    )

    dists = normalized_matrix[np.triu_indices_from(normalized_matrix, k=1)]

    seq_to_origin = {}
    for filtered_read, true_ids in origin_map.items():
        seq_to_origin[filtered_read] = true_ids

    # -----------------------------
    # PARAMETER SWEEP
    # -----------------------------
    for threshold in cluster_thresholds:
        for method in linkage_methods:

            print(f"  Running threshold={threshold}, method={method}")

            cluster_ids = cluster_sequences_dynamic(dists, threshold, method)
            consensus_records = build_consensus_from_lists(sequences, cluster_ids, abundance_dict)

            # Output filenames include parameters
            #Assuming the base path for the processed noisy sequence dataset is xxx
            output_consensus_file = (
                f"xxx/consensus_sequences_ITS_{tag}_thr{threshold}_meth{method}_abundance.fasta"
            )

            save_consensus_with_origin(
                consensus_records,
                seq_to_origin,
                output_consensus_file
            )

            print(f"Finished")