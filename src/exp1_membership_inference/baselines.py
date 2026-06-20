"""
Baseline membership-inference methods for Experiment 1.

Implements the two canonical prior methods that the paper cites:
  - Shringarpure & Bustamante (2015) AJHG likelihood ratio test (LRT)
  - Raisaro et al. (2017) JAMIA score function

Both accept the same (maf_vector, target_genotype) interface and return a
scalar score, higher = more likely member, compatible with sklearn's AUC
evaluation across the same k-SNP sweep as the logistic regression attack.
"""

import numpy as np
from typing import Sequence


def shringarpure_bustamante_lrt(
    maf_vector: np.ndarray,
    target_genotype: np.ndarray,
) -> float:
    """
    Likelihood ratio test statistic per Shringarpure & Bustamante, AJHG 2015.

    Under H1 (target is member of the reference panel), the genotype likelihood
    at each SNP j is a function of the allele frequency *with* the target included.
    Under H0 (target is not a member), it is the population MAF.

    For a diploid genotype g_j in {0,1,2} and reference MAF p_j, the LRT
    statistic is the log likelihood ratio summed across all SNPs:

        LRT = sum_j [ log P(g_j | p_j_with) - log P(g_j | p_j_without) ]

    We approximate p_j_with = (2*N*p_j + g_j) / (2*N + 2) for large N.
    For beacon queries the caller passes the beacon MAF as p_j (computed
    from the reference panel excluding the query target where possible).

    Parameters
    ----------
    maf_vector : array of shape (k,)
        Minor allele frequencies estimated from the reference cohort.
    target_genotype : array of shape (k,)
        Dosage genotype of the query individual, values in {0, 1, 2}.

    Returns
    -------
    float
        LRT score. Higher values indicate stronger evidence of membership.
    """
    p = np.clip(maf_vector, 1e-6, 1 - 1e-6)
    g = np.asarray(target_genotype, dtype=float)

    def _diploid_log_prob(dosage, freq):
        q = 1.0 - freq
        probs = np.where(
            dosage == 0, q**2,
            np.where(dosage == 1, 2 * freq * q, freq**2)
        )
        return np.log(np.clip(probs, 1e-300, None))

    # H1: frequency estimated as if target is included (N assumed large → minimal shift)
    # We use a conservative N=2504 to match 1000G cohort size; callers may override
    # by passing beacon-adjusted MAF directly.
    N = 2504
    p_with = (2 * N * p + g) / (2 * N + 2)
    p_with = np.clip(p_with, 1e-6, 1 - 1e-6)

    ll_h1 = _diploid_log_prob(g, p_with)
    ll_h0 = _diploid_log_prob(g, p)

    return float(np.sum(ll_h1 - ll_h0))


def raisaro_score(
    maf_vector: np.ndarray,
    target_genotype: np.ndarray,
) -> float:
    """
    Score function per Raisaro et al., JAMIA 2017.

    The Raisaro score sums, over each SNP, the log-ratio of the probability
    of observing the target allele under the hypothesis that the beacon
    returned TRUE (target is a member) versus FALSE:

        S = sum_j [ g_j * log(p_j / (1 - p_j)) ]

    This is equivalent to a log-odds weighting of the target's allele dosage
    by the population log-odds at each SNP, and can be viewed as a
    log-likelihood score under a logistic model with fixed per-SNP weights.

    Parameters
    ----------
    maf_vector : array of shape (k,)
        Minor allele frequencies from the reference beacon cohort.
    target_genotype : array of shape (k,)
        Dosage genotype of the query individual, values in {0, 1, 2}.

    Returns
    -------
    float
        Raisaro score. Higher values indicate stronger evidence of membership.
    """
    p = np.clip(maf_vector, 1e-6, 1 - 1e-6)
    g = np.asarray(target_genotype, dtype=float)
    log_odds = np.log(p / (1.0 - p))
    return float(np.dot(g, log_odds))


def compute_baseline_aucs(
    X: np.ndarray,
    y: np.ndarray,
    maf_vector: np.ndarray,
) -> dict:
    """
    Compute AUC for both baseline methods on a given feature matrix.

    Parameters
    ----------
    X : array of shape (n_samples, k)
        Genotype feature matrix (dosage values 0/1/2).
    y : array of shape (n_samples,)
        Binary membership labels.
    maf_vector : array of shape (k,)
        Per-SNP minor allele frequencies.

    Returns
    -------
    dict with keys 'sb_lrt_auc' and 'raisaro_auc'.
    """
    from sklearn.metrics import roc_auc_score

    sb_scores = np.array([shringarpure_bustamante_lrt(maf_vector, row) for row in X])
    ra_scores = np.array([raisaro_score(maf_vector, row) for row in X])

    return {
        "sb_lrt_auc": roc_auc_score(y, sb_scores),
        "raisaro_auc": roc_auc_score(y, ra_scores),
    }
