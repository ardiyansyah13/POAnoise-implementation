import random
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# -----------------------------
# CONFIGURATION
# -----------------------------
input_fasta = "the path to your clean sequence dataset (FASTA file)"

output_fastq = "the path to the output noisy sequence daataset (FASTQ file)"

real_fastq = "the path to the empirical FASTQ file for sampling quality score profile"

abundance_file = "the path to the file enumerating the abundance of the output noisy sequence dataset"

origin_file = "the path to the file describing the map connecting each noisy sequence in the output to the sequence in the clean sequence dataset"

# abundance classes for the 3 subgroups (enter either [1,4,16], [4,16,64], or [16,64,256])
ABUNDANCE_CLASSES = [16,64, 256]

# -----------------------------
# FUNCTIONS
# -----------------------------

def introduce_noise(sequence, phred_scores):
    """
    Introduce substitution errors according
    to the sampled quality scores.
    """

    bases = ['A', 'C', 'G', 'T']

    noisy_sequence = []

    for base, q in zip(sequence, phred_scores):

        error_prob = 10 ** (-q / 10)

        if random.random() < error_prob:

            noisy_base = random.choice(
                [b for b in bases if b != base]
            )

            noisy_sequence.append(noisy_base)

        else:

            noisy_sequence.append(base)

    return ''.join(noisy_sequence)


def parse_fastq_quality_pool(file_path):
    """
    Extract all quality profiles
    from a real FASTQ file.
    """

    quality_pool = []

    for record in SeqIO.parse(file_path, "fastq"):

        qs = record.letter_annotations["phred_quality"]

        quality_pool.append(qs)

    return quality_pool


def sample_quality_scores(length, quality_pool):
    """
    Randomly sample ONE real quality profile.
    """

    sampled_qs = random.choice(quality_pool)

    # Match target sequence length
    if len(sampled_qs) >= length:

        return sampled_qs[:length]

    else:

        extra = random.choices(
            sampled_qs,
            k=length - len(sampled_qs)
        )

        return sampled_qs + extra


# -----------------------------
# LOAD REAL QUALITY SCORES
# -----------------------------
quality_scores_pool = parse_fastq_quality_pool(real_fastq)

print(f"Loaded {len(quality_scores_pool)} real quality profiles")

# -----------------------------
# LOAD CLEAN SEQUENCES
# -----------------------------
clean_records = list(
    SeqIO.parse(input_fasta, "fasta")
)

n = len(clean_records)

print(f"Loaded {n} clean sequences")

# -----------------------------
# SPLIT INTO 3 SUBGROUPS
# -----------------------------
group_sizes = [
    n // 3,
    n // 3,
    n - 2 * (n // 3)
]

groups = []

start = 0

for size in group_sizes:

    groups.append(
        clean_records[start:start + size]
    )

    start += size

# -----------------------------
# GENERATE NOISY READS
# -----------------------------
noisy_records = []

abundance_dict = defaultdict(int)

for group, abundance_class in zip(groups, ABUNDANCE_CLASSES):

    print(f"\nProcessing abundance class {abundance_class}")

    for record in group:

        true_id = record.id

        true_seq = str(record.seq)

        # Generate copies according to abundance
        for _ in range(abundance_class):

            # -----------------------------------
            # STEP 1:
            # sample quality profile
            # -----------------------------------
            sampled_qs = sample_quality_scores(
                len(true_seq),
                quality_scores_pool
            )

            # -----------------------------------
            # STEP 2:
            # introduce sequencing noise
            # -----------------------------------
            noisy_seq = introduce_noise(
                true_seq,
                sampled_qs
            )

            # Create FASTQ record
            noisy_record = SeqRecord(
                Seq(noisy_seq),
                id="temp",
                description=(
                    f"origin={true_id};"
                    f"abundance_class={abundance_class}"
                ),
                letter_annotations={
                    "phred_quality": sampled_qs
                }
            )

            noisy_records.append(noisy_record)

            abundance_dict[noisy_seq] += 1

# -----------------------------
# SHUFFLE READS
# -----------------------------
random.shuffle(noisy_records)

# -----------------------------
# RENUMBER READS
# -----------------------------
origin_tracking = []

for idx, record in enumerate(noisy_records, start=1):

    new_id = f"read_{idx}"

    record.id = new_id
    record.name = new_id

    # Recover original true sequence ID
    true_id = (
        record.description
        .split("origin=")[1]
        .split(";")[0]
    )

    origin_tracking.append(
        (new_id, true_id)
    )

# -----------------------------
# WRITE NOISY FASTQ
# -----------------------------
with open(output_fastq, "w") as out_f:

    SeqIO.write(noisy_records, out_f, "fastq")

print(f"\nNoisy FASTQ saved to:")
print(output_fastq)

# -----------------------------
# SAVE ORIGIN MAP
# -----------------------------
with open(origin_file, "w", newline="") as f:

    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        "Noisy_Read_ID",
        "True_Sequence_ID"
    ])

    writer.writerows(origin_tracking)

print(f"\nOrigin map saved to:")
print(origin_file)

# -----------------------------
# SAVE ABUNDANCE TABLE
# -----------------------------
with open(abundance_file, "w", newline="") as tsvfile:

    writer = csv.writer(tsvfile, delimiter='\t')

    writer.writerow([
        "Sequence",
        "Abundance"
    ])

    for seq, abundance in sorted(
        abundance_dict.items(),
        key=lambda x: x[1],
        reverse=True
    ):

        writer.writerow([seq, abundance])

print(f"\nAbundance table saved to:")
print(abundance_file)

# -----------------------------
# REPORT
# -----------------------------
print("\nSUMMARY")
print(f"Total noisy reads: {len(noisy_records)}")
print(f"Unique noisy sequences: {len(abundance_dict)}")

