"""
Cohort and membership label construction for Experiment 1.

Two strategies are provided:
  - construct_label_external (CORRECT): derives membership from 1000G population metadata.
    Labels have zero positional or LD relationship to any SNP indices used as features.
  - construct_label_circular_DEPRECATED (WRONG): the original circular method kept solely
    for ablation comparison. DO NOT use in main experiments.
"""

import numpy as np
import pandas as pd


def construct_label_external(panel_path: str, superpopulation: str = "EUR", seed: int = 42):
    """
    Build a binary membership label from 1000 Genomes population metadata.

    Loads the panel file (tab-separated, columns: sample, pop, super_pop, gender),
    filters to the requested superpopulation, then randomly splits those samples
    50/50 into study-group (label=1) and control (label=0).

    Returns
    -------
    sample_ids : list[str]
        Ordered sample IDs matching the returned label vector.
    y : np.ndarray of int, shape (n_samples_in_superpop,)
        Binary membership label. Independent of all SNP values.
    """
    panel = pd.read_csv(panel_path, sep="\t")
    # Column names vary across 1000G panel files; normalise to lowercase.
    panel.columns = [c.lower() for c in panel.columns]
    # Accept 'super_pop' or 'superpopulation'
    sp_col = next(c for c in panel.columns if "super" in c)
    sample_col = next(c for c in panel.columns if "sample" in c)

    subset = panel[panel[sp_col] == superpopulation].copy()
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(subset))
    half = len(subset) // 2
    labels = np.zeros(len(subset), dtype=int)
    labels[idx[:half]] = 1  # study group

    return subset[sample_col].tolist(), labels


# ---------------------------------------------------------------------------
# DEPRECATED — kept for ablation/regression comparison only.
# This method derives y from SNP signal in columns 40–60, then draws features
# from immediately adjacent columns 60+, which are in LD with 40–60 by
# construction. The classifier recovers the label rule via LD, not beacon risk.
# ---------------------------------------------------------------------------

def construct_label_circular_DEPRECATED(X_all: np.ndarray):
    """
    DEPRECATED circular cohort construction. For ablation comparison only.

    Do NOT use in main experiments. Label y is derived from SNPs 40–60;
    features are drawn from SNPs 60+ which are LD-adjacent — this inflates
    AUC and does not measure beacon disclosure risk.
    """
    import warnings
    warnings.warn(
        "construct_label_circular_DEPRECATED is for ablation only and must NOT "
        "be used in main experiments. It produces artificially inflated AUC via LD.",
        DeprecationWarning,
        stacklevel=2,
    )
    signal_snps = X_all[:, 40:60].sum(axis=1)
    threshold = np.percentile(signal_snps, 50)
    y = (signal_snps > threshold).astype(int)
    return y


# ---------------------------------------------------------------------------
# Constants used by attack.py to define non-overlapping feature windows
# ---------------------------------------------------------------------------

# Label-construction SNP range (DEPRECATED method only — documented here so
# test_no_circularity.py can assert features never touch this range).
CIRCULAR_LABEL_SNP_RANGE = (40, 60)
