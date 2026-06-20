# Genomic Beacons Under Siege
### Quantifying Re-Identification Risk from Summary Statistics and the Performance Cost of Homomorphic Encryption

> **IEEE BIBM 2026 Submission**
> **Dataset:** 1000 Genomes Project Phase 3, Chromosome 22 (N = 2,504)
> **FHE Library:** TenSEAL (CKKS scheme)

---

## Overview

This repository contains the full experimental pipeline characterizing the "privacy cliff" in genomic beacons and a domain-optimized FHE mitigation via parameter compaction.

| # | Experiment | Script | Key Result |
|---|---|---|---|
| 1 | Membership Inference Attack (MIA) | `src/exp1_membership_inference/attack.py` | AUC phase transition across k-SNP sweep |
| 2 | FHE Overhead Baseline | `src/exp2_fhe_overhead/benchmark.py` | End-to-end latency and dual overhead metrics |
| 3 | Parameter Compaction + Factorial | `src/exp3_parameter_compaction/` | Latency reduction attributed via 2×2 factorial |

---

## Repository Structure

```
GenomicBeaconsUnderSiege/
├── requirements.txt
├── data/
│   ├── raw/          # 1000G chr22 VCF or panel file goes here
│   └── processed/    # derived/filtered genotype matrices
├── src/
│   ├── exp1_membership_inference/
│   │   ├── attack.py              # main MIA attack (logistic regression + baselines)
│   │   ├── baselines.py           # Shringarpure-Bustamante LRT, Raisaro score
│   │   └── cohort_construction.py # external (metadata-derived) label construction
│   ├── exp2_fhe_overhead/
│   │   └── benchmark.py           # CKKS baseline overhead benchmark
│   ├── exp3_parameter_compaction/
│   │   ├── compaction.py          # compacted config (L=1, scale=2^21)
│   │   └── factorial.py           # 2×2 factorial: {L=1,L=2} × {scale=2^21,2^40}
│   └── analysis/
│       └── cue.py                 # CUE table with cross-experiment sanity checks
├── results/
│   ├── exp1/
│   ├── exp2/
│   ├── exp3/
│   └── analysis/
├── figures/
└── tests/
    └── test_no_circularity.py     # regression guard against circular cohort construction
```

---

## Experiment Descriptions

### Experiment 1 — Membership Inference Attack

Measures the beacon re-identification risk (AUC) of a logistic-regression MIA as a
function of the number of SNPs queried (k). Also evaluates the two canonical prior
methods as baselines:

- **Shringarpure & Bustamante (2015) AJHG** likelihood ratio test
- **Raisaro et al. (2017) JAMIA** score function

**Cohort construction:** Membership labels are derived exclusively from 1000 Genomes
population metadata (the panel file mapping sample IDs to superpopulation, e.g. EUR).
A single superpopulation is selected and randomly split 50/50 into study-group
(label=1) and control (label=0). The label has **zero positional or LD relationship**
to any SNP indices used as features — the label cannot be recovered from the genotype
matrix, only from external metadata.

A deprecated circular label method (deriving y from SNP signal at positions 40–60,
then drawing features from LD-adjacent positions 60+) is preserved in
`cohort_construction.py` marked `DEPRECATED`, for ablation comparison only.

**SNP sweep:** k = [5, 10, 20, 50, 60, 70, 80, 90, 100, 110, 200] (finer grid
around the inflection region per reviewer request).

### Experiment 2 — FHE Overhead Baseline

Benchmarks CKKS homomorphic encryption overhead for genomic PRS computation
at N=2,504 individuals × 200 SNPs, using the baseline parameter set
(depth L=2, scale 2^40, coeff_mod=[60,40,40,60]).

**Overhead reporting:** Two overhead figures are computed and reported explicitly:
- **Computation-only overhead:** `total_comp_time / pt_total_time`
- **End-to-end overhead (enc + comp):** `(total_enc_time + total_comp_time) / pt_total_time`

The abstract's cited overhead figure should reference end-to-end overhead consistently.

**Data note:** This benchmark uses synthetic genotype-shaped data
(`np.random.randint(0, 3, ...)`) rather than loading real 1000G genotypes.
CKKS timing is determined by ciphertext structure and cryptographic parameters,
not by the plaintext values being encrypted — this is standard and accepted
practice in FHE benchmarking. The data has the same shape and value range
({0,1,2} diploid dosage) as real 1000G chr22 genotypes. Set `USE_REAL_DATA=True`
in `benchmark.py` if the real genotype matrix is available.

### Experiment 3 — Parameter Compaction and Factorial Decomposition

Two scripts:

**`compaction.py`** — benchmarks the compacted configuration (L=1, scale=2^21,
coeff_mod=[40,21,40]) in isolation.

**`factorial.py`** — runs all four configurations of the 2×2 factorial
{L=1, L=2} × {scale=2^21, 2^40} and computes marginal attribution:

| Config | Depth (L) | Modulus chain | Scale |
|---|---|---|---|
| A (baseline) | 2 | [60,40,40,60] | 2^40 |
| B (depth only) | 1 | [40,40,40] | 2^40 |
| C (scale only) | 2 | [60,21,21,60] | 2^21 |
| D (compacted) | 1 | [40,21,40] | 2^21 |

This decomposes how much of the total latency reduction and MAE increase
comes from depth reduction vs. scale reduction (and their interaction).

---

## Setup

### Requirements

- Python 3.10+
- RAM ≥ 8 GB (16 GB recommended for full-scale FHE benchmarks)

```bash
python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

Pinned versions:
```
tenseal==0.3.14, pandas==2.2.2, numpy==1.26.4,
scikit-learn==1.5.0, matplotlib==3.9.0, tqdm==4.66.4
```

### Data Setup (Experiment 1 only)

Download the 1000 Genomes chr22 genotype matrix and panel file:

```bash
# Panel file (sample -> superpopulation mapping)
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
    -O data/raw/integrated_call_samples_v3.20130502.ALL.panel

# Genotype matrix (VCF → TSV conversion, chr22)
# See data/raw/README for conversion steps
```

---

## Reproducing Results

Run experiments in order (Exp 2 and 3 are independent; Exp 1 requires the panel file):

```bash
# Experiment 2 — FHE overhead baseline (no data needed)
python -m src.exp2_fhe_overhead.benchmark

# Experiment 3 — Compaction (no data needed)
python -m src.exp3_parameter_compaction.compaction

# Experiment 3 — Full factorial decomposition (no data needed)
python -m src.exp3_parameter_compaction.factorial

# Experiment 1 — MIA attack (requires panel file and genotype TSV)
python -m src.exp1_membership_inference.attack \
    --genotype-tsv data/genotype_matrix.tsv \
    --panel data/raw/integrated_call_samples_v3.20130502.ALL.panel \
    --superpopulation EUR

# CUE table (requires results from all three experiments)
python -m src.analysis.cue
```

### Running Tests

```bash
python -m pytest tests/ -v
```

The circularity regression test (`tests/test_no_circularity.py`) will fail loudly
if the deprecated circular cohort construction method is used in place of the
external label construction.
