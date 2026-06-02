# POAnoise: A Graph-based Denoising Pipeline for Amplicon Sequencing Data

POAnoise is a graph-based denoising framework for amplicon sequencing data that combines Partial Order Alignment (POA), abundance-aware clustering, and heaviest-bundle consensus extraction to reconstruct biological sequences from noisy sequencing reads.

This repository includes:

- Synthetic dataset generation
- Read preprocessing and filtering
- POAnoise implementation
- Consensus reconstruction
- Sequence similarity evaluation
- Benchmarking across datasets and parameters

---

# Datasets

Two clean reference datasets are used:

- ITS: `clean_sequences_ITS_75.fasta`
- 16S: `clean_sequences_16S_75.fasta`

Each dataset contains 75 reference sequences.

---

# Quality Score Pool

Real sequencing quality profiles are sampled from:

`quality_scores.fastq.gz`

These empirical Phred quality scores are used to simulate realistic sequencing errors.

---

# Synthetic Data Generation

Script:
`noise_generation.py`

This script generates noisy sequencing reads from clean reference sequences.

## Abundance Classes

Three abundance regimes are used:

- Dataset A: [1, 4, 16]
- Dataset B: [4, 16, 64]
- Dataset C: [16, 64, 256]

Each reference sequence is assigned to one abundance class.

## Error Model

For each base with quality score Q:

P(error) = 10^(-Q/10)

If an error occurs, the base is replaced with a random alternative nucleotide.

Only substitution errors are simulated.

---

# Generated Files

## Noisy Reads
`noisy_sequences_<dataset>.fastq`

## Origin Mapping
`noisy_origin_map_<dataset>.tsv`

Maps noisy reads to original reference sequences.

## Abundance Table
`noisy_abundance_<dataset>.tsv`

Contains counts of unique noisy sequences.

---

# Preprocessing

Script:
`preprocessing.py`

Preprocessing is performed using USEARCH.

## Quality Filtering

Reads are filtered using:

- fastq_maxee = 2

Command:
usearch -fastq_filter -fastq_maxee 2

## Dereplication

Performed using:

usearch -fastx_uniques -sizeout

---

# Preprocessing Outputs

- `noisy_preprocessed_<dataset>.fastq`
- `noisy_preprocessed_abundance_<dataset>.tsv`
- `filtered_origin_map_<dataset>.tsv`
- `derep_origin_map_<dataset>.tsv`

---

# POAnoise Algorithm

Script:
`POAnoise.py`

## Overview

POAnoise constructs a Partial Order Graph (POG) and performs iterative sequence integration using dual-affine Needleman–Wunsch alignment.

---

## Alignment Parameters

- Match = 1  
- Mismatch = -1  
- Gap open 1 = -1  
- Gap extend 1 = -1  
- Gap open 2 = -5  
- Gap extend 2 = 0  

---

## Workflow

### 1. Build POA Graph
The first sequence initializes the graph, and all other sequences are added via alignment.

### 2. Edge-Sequence Matrix
Construct abundance-weighted edge support matrix.

### 3. Distance Matrix
Compute cosine distances between sequences.

### 4. Clustering
Hierarchical clustering using single linkage.

Thresholds:
0.01 to 0.96 (grid search)

### 5. Consensus Extraction
Consensus sequences are extracted using the heaviest-bundle algorithm.

---

## Output Consensus Files

`consensus_sequences_<dataset>_thr<threshold>_methsingle_abundance.fasta`

Each FASTA header contains:

- Consensus ID
- Weight
- Origin sequence IDs

---

# Sequence Similarity Evaluation

Script:
`sequence_similarity.py`

This script compares POAnoise consensus sequences to clean reference sequences.

---

## Alignment

Global alignment using:

- Match = 1
- Mismatch = -1
- Gap open = -2
- Gap extend = -1

---

## Outputs

### Similarity Matrix
`similarity_matrix_<dataset>.tsv`

Pairwise similarity between consensus and reference sequences.

### Best Match Report
`best_matches_<dataset>.txt`

Includes:
- Best matching reference
- Similarity score
- Full alignment visualization

### High-Quality References
`best_matches_<dataset>_best_refs_high_quality.fasta`

Contains references reconstructed with 100% similarity.

---

# Requirements

Install dependencies:

pip install networkx numpy pandas biopython scipy scikit-learn matplotlib

External dependency:
- USEARCH v11+

---

# Reproduction Pipeline

1. Generate noisy reads:
python noise_generation.py

2. Preprocess reads:
python preprocessing.py

3. Run POAnoise:
python POAnoise.py

4. Compute similarity:
python sequence_similarity.py

---

# Citation

If you use this repository, please cite:

POAnoise: A Graph-based Denoising Pipeline for Amplicon Sequencing Data (manuscript in preparation)
