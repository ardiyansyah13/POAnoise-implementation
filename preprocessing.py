import subprocess
import re
import csv
import os
from Bio import SeqIO
from collections import defaultdict

usearch_exe = "the path to the executable usearch"

def usearch_fastq_filter_with_abundance(
        input_file,
        filtered_fastq,
        abundance_tsv,
        origin_map_input,
        filtered_origin_output,
        derep_origin_output,
        maxee):

    os.makedirs(os.path.dirname(filtered_fastq), exist_ok=True)

    # -----------------------------
    # COUNT BEFORE
    # -----------------------------
    reads_before = sum(1 for _ in open(input_file)) // 4
    print(f"Number of sequences before filtering: {reads_before}")

    # -----------------------------
    # FILTER
    # -----------------------------
    subprocess.run([
        usearch_exe, "-fastq_filter", input_file,
        "-fastq_maxee", str(maxee),
        "-fastqout", filtered_fastq
    ], check=True)

    reads_after = sum(1 for _ in open(filtered_fastq)) // 4
    print(f"Number of sequences after filtering: {reads_after}")

    # -----------------------------
    # LOAD FULL ORIGIN MAP
    # -----------------------------
    full_origin_map = {}
    with open(origin_map_input) as f:
        next(f)
        for line in f:
            read_id, true_id = line.strip().split("\t")
            full_origin_map[read_id] = true_id

    # -----------------------------
    # BUILD FILTERED-ONLY ORIGIN MAP
    # -----------------------------
    filtered_origin_map = {}
    seq_to_true = defaultdict(list)

    for record in SeqIO.parse(filtered_fastq, "fastq"):
        read_id = record.id

        if read_id not in full_origin_map:
            continue

        true_id = full_origin_map[read_id]
        filtered_origin_map[read_id] = true_id

        seq_to_true[str(record.seq)].append(true_id)

    # -----------------------------
    # SAVE FILTERED ORIGIN MAP
    # -----------------------------
    with open(filtered_origin_output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Filtered_Read_ID", "True_Sequence_ID"])
        for rid, tid in filtered_origin_map.items():
            writer.writerow([rid, tid])

    print(f"Filtered origin map saved to {filtered_origin_output}")

    # -----------------------------
    # DEREPLICATE (FILTERED ONLY)
    # -----------------------------
    derep_fasta = filtered_fastq.replace(".fastq", "_uniques.fasta")

    subprocess.run([
        usearch_exe, "-fastx_uniques", filtered_fastq,
        "-fastaout", derep_fasta,
        "-sizeout"
    ], check=True)

    # -----------------------------
    # PARSE DEREP FASTA
    # -----------------------------
    derep_data = []

    with open(derep_fasta) as f:
        seq = ""
        size = 1
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if seq:
                    derep_data.append((seq, size))
                    seq = ""
                match = re.search(r"size=(\d+)", line)
                size = int(match.group(1)) if match else 1
            else:
                seq += line

        if seq:
            derep_data.append((seq, size))

    # -----------------------------
    # SAVE ABUNDANCE TSV (FILTERED ONLY)
    # -----------------------------
    with open(abundance_tsv, "w") as f_out:
        f_out.write("Sequence\tAbundance\n")
        for seq, count in derep_data:
            f_out.write(f"{seq}\t{count}\n")

    print(f"Abundance TSV saved to {abundance_tsv}")

    # -----------------------------
    # SAVE DEREP ORIGIN MAP (FILTERED ONLY)
    # -----------------------------
    with open(derep_origin_output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Derep_Sequence", "Abundance", "True_Origin_Counts"])

        for seq, size in derep_data:
            true_counts = defaultdict(int)
            for t in seq_to_true.get(seq, []):
                true_counts[t] += 1

            writer.writerow([seq, size, dict(true_counts)])

    print(f"Dereplicated origin map saved to {derep_origin_output}")
    print(f"Dereplicated FASTA saved to {derep_fasta}")

# ---------------------------------------
# AUTOMATIC PROCESSING FOR:
# ITS_A, ITS_B, ITS_C,
# 16S_A, 16S_B, 16S_C
# ---------------------------------------

DATASETS = ["ITS", "16S"]

SUFFIXES = ["A", "B", "C"]

for dataset in DATASETS:

    for suffix in SUFFIXES:

        dataset_name = f"{dataset}_{suffix}"

        print(f"\nProcessing dataset: {dataset_name}")

        # ---------------------------------------
        # FILE PATHS
        #Assuming the base path for the noisys sequence dataset is xxx
        # ---------------------------------------
        input_file = (
            f"xxx/noisy_sequences_{dataset_name}.fastq"
        )

        filtered_fastq = (
            f"xxx/noisy_preprocessed_{dataset_name}.fastq"
        )

        abundance_tsv = (
            f"xxx/noisy_preprocessed_abundance_{dataset_name}.tsv"
        )

        origin_map_input = (
            f"xxx/noisy_origin_map_{dataset_name}.tsv"
        )

        filtered_origin_output = (
            f"xxx/filtered_origin_map_{dataset_name}.tsv"
        )

        derep_origin_output = (
            f"xxx/derep_origin_map_{dataset_name}.tsv"
        )

        # ---------------------------------------
        # RUN PREPROCESSING
        # ---------------------------------------
        usearch_fastq_filter_with_abundance(
            input_file=input_file,
            filtered_fastq=filtered_fastq,
            abundance_tsv=abundance_tsv,
            origin_map_input=origin_map_input,
            filtered_origin_output=filtered_origin_output,
            derep_origin_output=derep_origin_output,
            maxee=2
        )

        print(f"Finished processing {dataset_name}")