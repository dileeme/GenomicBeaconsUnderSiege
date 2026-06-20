"""
Experiment 2: FHE Overhead Baseline.

DATA NOTE: This benchmark uses synthetic genotype-shaped data
(np.random.randint(0, 3, ...)) rather than real 1000G genotypes. CKKS timing
is determined by ciphertext structure and cryptographic parameters, not by
the plaintext values being encrypted — this is a standard and accepted practice
in FHE benchmarking (see e.g. Cheon et al. 2017; Chen et al. 2019). The
synthetic data has the same shape and dtype as real 1000G chr22 genotypes
(N individuals x K SNPs, values in {0,1,2}), so all timing measurements are
directly representative of the real-data case.

If real genotype data is available at data/genotype_matrix.tsv, set
USE_REAL_DATA=True to load it instead.
"""

import os
import csv
import time
import numpy as np
import tenseal as ts
from tqdm import tqdm

USE_REAL_DATA = False  # Set True if data/genotype_matrix.tsv is available

N_INDIVIDUALS = 2504  # matches 1000 Genomes Project phase 3 cohort size
N_SNPS = 200

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "exp2")


def load_data():
    if USE_REAL_DATA:
        import pandas as pd
        tsv_path = os.path.join(os.path.dirname(__file__), "..", "..", "data", "genotype_matrix.tsv")
        print(f"Loading real genotype data from {tsv_path}...")
        df = pd.read_csv(tsv_path, sep="\t", header=None, nrows=N_SNPS)
        df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")
        X = df.T.replace({
            "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2,
            "0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2,
        }).apply(pd.to_numeric, errors="coerce").fillna(0).values
        return X[:N_INDIVIDUALS, :N_SNPS]
    else:
        # Synthetic genotype-shaped data; timing is data-value-independent under CKKS.
        return np.random.randint(0, 3, size=(N_INDIVIDUALS, N_SNPS))


def run_full_scale_study():
    print(f"--- FHE Overhead Baseline: N={N_INDIVIDUALS} individuals x {N_SNPS} SNPs ---")
    data_source = "real genotypes" if USE_REAL_DATA else "synthetic genotype-shaped data (CKKS timing is plaintext-independent)"
    print(f"    Data source: {data_source}")

    data = load_data()
    weights = np.random.uniform(0.001, 0.05, size=(N_SNPS,))

    # Baseline CKKS context: depth-2, scale 2^40
    context = ts.context(
        ts.SCHEME_TYPE.CKKS,
        poly_modulus_degree=8192,
        coeff_mod_bit_sizes=[60, 40, 40, 60],
    )
    context.global_scale = 2**40
    context.generate_galois_keys()

    # Plaintext baseline
    start_pt = time.time()
    pt_results = data.dot(weights)
    pt_total_time = time.time() - start_pt

    # Encryption phase
    enc_start = time.time()
    encrypted_rows = []
    for row in tqdm(data.tolist(), desc="Encrypting"):
        encrypted_rows.append(ts.ckks_vector(context, row))
    total_enc_time = time.time() - enc_start

    # Computation phase
    comp_start = time.time()
    encrypted_scores = []
    for enc_row in tqdm(encrypted_rows, desc="Computing encrypted PRS"):
        encrypted_scores.append(enc_row.dot(weights.tolist()))
    total_comp_time = time.time() - comp_start

    dec_results = [score.decrypt()[0] for score in encrypted_scores]
    mae = np.mean(np.abs(pt_results - dec_results))

    # Fix 5: report both overhead definitions explicitly and unambiguously.
    # The paper's abstract claim should reference end_to_end_overhead consistently.
    # computation_only_overhead is a secondary breakdown showing computation alone.
    computation_only_overhead = total_comp_time / pt_total_time
    end_to_end_overhead = (total_enc_time + total_comp_time) / pt_total_time

    print("\n" + "=" * 55)
    print(f"EXPERIMENT 2 RESULTS  (N={N_INDIVIDUALS}, K={N_SNPS})")
    print("-" * 55)
    print(f"Plaintext Processing:              {pt_total_time * 1000:>10,.2f} ms")
    print(f"Total Encryption Time:             {total_enc_time:>10,.2f} s")
    print(f"Total FHE Computation Time:        {total_comp_time:>10,.2f} s")
    print(f"Accuracy (MAE):                    {mae:>10.10f}")
    print("-" * 55)
    print(f"Overhead (computation only):       {computation_only_overhead:>10,.0f}x")
    print(f"Overhead (end-to-end enc+comp):    {end_to_end_overhead:>10,.0f}x")
    print("=" * 55)
    print("NOTE: 'end-to-end' overhead is the figure that should be cited")
    print("      in the abstract. 'computation only' excludes encryption time.")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    csv_path = os.path.join(RESULTS_DIR, "benchmark_results.csv")
    rows = [
        ("N Individuals", N_INDIVIDUALS),
        ("N SNPs", N_SNPS),
        ("Data Source", "real" if USE_REAL_DATA else "synthetic"),
        ("Plaintext Time (ms)", round(pt_total_time * 1000, 4)),
        ("Total Encryption Time (s)", round(total_enc_time, 4)),
        ("Total FHE Computation (s)", round(total_comp_time, 4)),
        ("Accuracy (MAE)", round(float(mae), 10)),
        ("Overhead Computation Only (x)", round(computation_only_overhead, 0)),
        ("Overhead End-to-End (x)", round(end_to_end_overhead, 0)),
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Metric", "Value"])
        writer.writerows(rows)
    print(f"\nResults saved: {csv_path}")


if __name__ == "__main__":
    run_full_scale_study()
