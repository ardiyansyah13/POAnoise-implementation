

"""
Compute sequence similarity between consensus and reference sequences.
Features -------- 
- Pairwise global alignment 
- Best-match extraction 
- Similarity matrix export 
- Alignment report generation 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio import SeqIO
from matplotlib.colors import LinearSegmentedColormap
from typing import Dict
import os
import re
from glob import glob


# ---------------------------
# Helper functions
# ---------------------------
def write_fasta(sequences: Dict[str, str], output_file: str):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


def highlight_mismatches(aligned_consensus: str, aligned_reference: str):
    """Return a list of mismatch positions (ignoring gaps)."""
    mismatches = []
    for i, (a, b) in enumerate(zip(aligned_consensus, aligned_reference)):
        if a != b and a != '-' and b != '-':
            mismatches.append(i)
    return mismatches


def align_sequences(consensus_seq: str, reference_seq: str):
    """Perform global alignment and compute similarity percentage."""
    primer_sequence = "enter the primer you are interested in"  

    match = 1
    mismatch = -1
    gap_open = -2
    gap_extend = -1

    alignment = pairwise2.align.globalms(
        consensus_seq, reference_seq,
        match, mismatch, gap_open, gap_extend,
        one_alignment_only=True
    )[0]

    aligned_consensus, aligned_reference, score, start, end = alignment

    # Fill first x bases if missing
    for i in range(min(x, len(aligned_consensus))):
        if aligned_consensus[i] == '-':
            aligned_consensus = aligned_consensus[:i] + primer_sequence[i] + aligned_consensus[i+1:]

    # Trim overhanging tails
    i = len(aligned_consensus)
    while i > 0:
        if aligned_consensus[i-1] != '-' and aligned_reference[i-1] != '-':
            break
        i -= 1
    aligned_consensus = aligned_consensus[:i]
    aligned_reference = aligned_reference[:i]

    # Compute similarity
    if len(aligned_consensus) == 0:
        similarity_percentage = 0.0
    else:
        num_matches = sum(1 for a, b in zip(aligned_consensus, aligned_reference) if a == b)
        similarity_percentage = (num_matches / len(aligned_consensus)) * 100

    return similarity_percentage, aligned_consensus, aligned_reference


def write_highlighted_alignment(consensus_id, ref_id, aligned_c, aligned_r, mismatches, similarity, f):
    """Write alignment with mismatches highlighted and similarity."""
    highlighted_c = ''.join([f"**{aligned_c[i]}**" if i in mismatches else aligned_c[i] for i in range(len(aligned_c))])
    highlighted_r = ''.join([f"**{aligned_r[i]}**" if i in mismatches else aligned_r[i] for i in range(len(aligned_r))])

    f.write(f"Consensus ID: {consensus_id}\n")
    f.write(f"Best Match Reference ID: {ref_id}\n")
    f.write(f"Similarity: {similarity:.2f}%\n")
    f.write(f"Aligned Consensus:   {highlighted_c}\n")
    f.write(f"Aligned Reference:   {highlighted_r}\n\n")


def read_fasta(file_path: str):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}


def read_fastq(file_path: str):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fastq")}


# ---------------------------
# Main computation
# ---------------------------

def compute_similarity_matrix(
    consensus_file: str,
    reference_file: str,
    output_matrix_file: str,
    output_best_match_file: str
):
    consensus_sequences = read_fasta(consensus_file)
    print(f"Loaded {len(consensus_sequences)} consensus sequences.")

    reference_sequences = read_fasta(reference_file)

    similarity_matrix = pd.DataFrame(index=consensus_sequences.keys(), columns=reference_sequences.keys())
    alignment_cache = {}
    updated_consensus_sequences = {}

    # Compute all pairwise similarities
    for cid, cseq in consensus_sequences.items():
        for rid, rseq in reference_sequences.items():
            sim, aligned_c, aligned_r = align_sequences(cseq, rseq)
            similarity_matrix.at[cid, rid] = sim
            alignment_cache[(cid, rid)] = (sim, aligned_c, aligned_r)

        # Save best match alignment
        best_rid = similarity_matrix.loc[cid].astype(float).idxmax()
        updated_consensus_sequences[cid] = alignment_cache[(cid, best_rid)][1]

    similarity_matrix = similarity_matrix.astype(float)
    similarity_matrix.to_csv(output_matrix_file)
    print(f"Similarity matrix saved to {output_matrix_file}")

    # Identify best matches
    match_list = []
    best_match_ids = set()
    for cid, row in similarity_matrix.iterrows():
        best_rid = row.idxmax()
        best_score = row.max()
        match_list.append((best_score, cid, best_rid))
        best_match_ids.add(best_rid)

    match_list.sort(reverse=True, key=lambda x: x[0])

    # Write alignments to file
    count_perfect = 0
    count_one_mismatch = 0
    with open(output_best_match_file, "w") as f:
        for _, cid, rid in match_list:
            sim, aligned_c, aligned_r = alignment_cache[(cid, rid)]
            mismatches = highlight_mismatches(aligned_c, aligned_r)
            write_highlighted_alignment(cid, rid, aligned_c, aligned_r, mismatches, sim, f)
            if sim == 100.0:
                count_perfect += 1
            elif len(mismatches) == 1 and sum(1 for a, b in zip(aligned_c, aligned_r) if a == '-' or b == '-') == 0:
                count_one_mismatch += 1
        f.write("\n--- Summary ---\n")
        f.write(f"Total perfect matches (100%): {count_perfect}\n")
        f.write(f"Total matches with exactly 1 mismatch (no gaps): {count_one_mismatch}\n")
    print(f"Best matches written to {output_best_match_file}")

    # ---------------------------
    # Extract high-quality reference sequences (>=99% sim)
    # ---------------------------
    best_reference_sequences = {}
    for cid, rid in [(cid, rid) for _, cid, rid in match_list]:
        sim, aligned_c, aligned_r = alignment_cache[(cid, rid)]
        if sim >= 100.0:
            best_reference_sequences[rid] = reference_sequences[rid]

    best_ref_fasta_file = output_best_match_file.replace(".txt", "_best_refs_high_quality.fasta")
    write_fasta(best_reference_sequences, best_ref_fasta_file)
    print(f"High-quality best-match reference sequences saved to {best_ref_fasta_file}")



if __name__ == "__main__":

    # ---------------------------
    # Folders
    # ---------------------------
    input_folder = "the path to the output sequences from POAnoise"
    reference_file = "the path to the clean sequence dataset"

    output_matrix_folder = "the path to the similarity matrix folder"
    output_bestmatch_folder = "the path to the best match file folder"

    os.makedirs(output_matrix_folder, exist_ok=True)
    os.makedirs(output_bestmatch_folder, exist_ok=True)

    # ---------------------------
    # Find all output files
    # ---------------------------
    #The following serve an example of how to find all output sequences from POAnoise 
    #for the ITS dataset (ITS_A, ITS_B, ITS_C) with clustering threshold thr and clustering method meth
    pattern = os.path.join(input_folder, "consensus_sequences_ITS_*_thr*_methsingle.fasta")
    input_files = glob(pattern)

    print(f"Found {len(input_files)} files.")

    # ---------------------------
    # Process each file
    # ---------------------------
    for denoised_file in input_files:

        filename = os.path.basename(denoised_file)

        # Extract parameters using regex
        match = re.match(r"consensus_sequences_ITS_([ABC])_thr0.(\d+)_methsingle_abundance.fasta", filename)

        if not match:
            print(f"Skipping {filename} (pattern mismatch)")
            continue

        X, Y = match.groups()

        tag = f"ITS_{X}_thr0.{Y}_methsingle"

        print(f"\nProcessing {tag}")

        output_matrix_file = os.path.join(
            output_matrix_folder,
            f"similarity_matrix_{tag}.tsv"
        )

        output_best_match_file = os.path.join(
            output_bestmatch_folder,
            f"best_matches_{tag}.txt"
        )

        compute_similarity_matrix(
            denoised_file,
            reference_file,
            output_matrix_file,
            output_best_match_file
        )

    print("\nAll analyses completed.")

   
