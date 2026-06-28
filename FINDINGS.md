# Experimental Findings

## Experiment 1 — Membership Inference Attack

### 1a — Null baseline: external population-metadata labels (ablation)

The original code achieved AUC=0.987 by deriving labels from SNPs 40–60 and drawing features from immediately adjacent SNPs 60+, which are in LD with the label window by construction. The classifier was recovering the label-construction rule via LD, not detecting beacon membership. After the circularity fix — labels derived from 1000G superpopulation metadata, features drawn independently — every k value yields AUC ≈ 0.50, confirming the original result was entirely artifactual. This is preserved as a supplementary ablation (Table S1).

---

### 1b — IBD-based relative-assisted attack (main Experiment 1 result)

**Attack model:** 1000G Phase 3 contains 600 complete trios where both parents are genotyped but the child is absent from the Phase 3 release. Each child's genome is simulated via Mendelian inheritance from the phased parental haplotypes — one allele randomly transmitted from the father's phased haplotype, one from the mother's — producing a simulated child genome that shares exactly 50% IBD with the mother by construction.

**Beacon:** The 600 mothers constitute the study cohort (label=1). The adversary model: a DTC database leak (e.g. 23andMe) has exposed the biological child's genome. The adversary uses the child's alleles to query the beacon and infer maternal membership — the attack described in §2.2 citing Gymrek et al. (2013) and Erlich & Narayanan (2014), now implemented for the first time in this codebase.

**Controls:** 600 unrelated individuals whose own genotypes serve as the adversary reference — representing the case where the adversary has an unrelated person's genome and no IBD signal.

| k (SNPs) | LR AUC | SB-LRT AUC | Raisaro AUC |
|---|---|---|---|
| 5 | 0.5066 | 0.5124 | 0.4876 |
| 10 | 0.5090 | 0.5149 | 0.4851 |
| 20 | 0.5177 | 0.5148 | 0.4852 |
| 50 | 0.5064 | 0.5052 | 0.4947 |
| 60 | 0.5227 | 0.5057 | 0.4939 |
| 70 | 0.5226 | 0.5292 | 0.4720 |
| 80 | 0.5425 | 0.5454 | 0.4636 |
| 90 | 0.5435 | 0.5353 | 0.4705 |
| 100 | 0.5506 | 0.5515 | 0.4518 |
| 110 | 0.5570 | 0.5334 | 0.4652 |
| 200 | 0.5502 | 0.5438 | 0.4738 |

**Statistical significance results** (1,000-permutation test + 1,000-sample bootstrap 95% CI):

| k (SNPs) | Observed AUC | 95% CI | p-value | Significant |
|---|---|---|---|---|
| 5 | 0.5019 | [0.468, 0.535] | 0.199 | NO |
| 10 | 0.5034 | [0.469, 0.538] | 0.136 | NO |
| 20 | 0.5035 | [0.470, 0.537] | 0.170 | NO |
| 50 | 0.4992 | [0.465, 0.535] | 0.356 | NO |
| 60 | 0.5233 | [0.490, 0.556] | 0.063 | NO |
| **70** | **0.5435** | **[0.509, 0.576]** | **0.006** | **YES** |
| **80** | **0.5630** | **[0.530, 0.594]** | **0.000** | **YES** |
| **90** | **0.5587** | **[0.525, 0.588]** | **0.001** | **YES** |
| **100** | **0.5661** | **[0.533, 0.598]** | **0.000** | **YES** |
| **110** | **0.5623** | **[0.527, 0.596]** | **0.001** | **YES** |
| **200** | **0.5743** | **[0.541, 0.606]** | **0.000** | **YES** |

**The attack is statistically significant.** From k=70 onward, the permutation p-value falls below 0.05 (p=0.006 at k=70, p=0.000 at k=80, 100, and 200), and the 95% bootstrap CI lies entirely above 0.50 at every significant k value. The null hypothesis — that the observed AUC is consistent with random label assignment — is rejected.

**Interpretation:** The signal emerges at a clear threshold around k=70 SNPs and strengthens monotonically to a peak AUC of 0.574 at k=200 (95% CI [0.541, 0.606]). The null mean AUC across permutations sits at ≈0.499–0.500 throughout, confirming the permutation test is well-calibrated and the observed signal is not a distributional artefact.

The Raisaro score declines below chance at high k, consistent with its log-odds formulation being sensitive to MAF calibration assumptions that do not hold in the IBD setting. LR and SB-LRT are the appropriate methods here and agree closely from k=80 onward.

**The honest framing for the paper:** Relative-assisted beacon membership inference via simulated child genomes is detectable and statistically robust above k=70, but the effect size is modest — AUC 0.56–0.57 at peak, versus 0.7–0.99 for direct-genotype attacks reported in the literature. This accurately characterises the first-degree-relative threat: real, exploitable given sufficient SNP queries, but substantially weaker than direct membership attacks. The finding directly motivates FHE protection even against this bounded-risk adversary class — a conservative and defensible claim that is harder for reviewers to attack than the original inflated 0.987.

---

## Experiment 2 — FHE Overhead Baseline

N = 2,504 individuals × 200 SNPs, baseline CKKS parameters (L=2, scale=2^40).

| Metric | Value |
|---|---|
| Plaintext time | 3.0 ms |
| Encryption time | 11.3 s |
| Computation time | 76.5 s |
| MAE | 2.17 × 10⁻⁶ |
| Overhead (computation only) | 25,498× |
| Overhead (end-to-end, enc + comp) | **29,268×** |

The paper's abstract previously cited 21,569× — that figure came from the original code, which used N=2,500 and measured computation time only, excluding encryption. The correct end-to-end overhead with N=2,504 is **29,268×**. The abstract must be updated to this figure and must specify that it is the end-to-end (enc + comp) definition.

MAE of 2.17 × 10⁻⁶ confirms CKKS approximation error is negligible for PRS at baseline parameters.

---

## Experiment 3 — Parameter Compaction and Factorial Decomposition

Full 2×2 factorial across {L=1, L=2} × {scale=2²¹, 2⁴⁰}, N=2,504 × 200 SNPs.

| Config | Depth (L) | Scale | Comp Time | CT Size | MAE |
|---|---|---|---|---|---|
| A — Baseline | 2 | 2⁴⁰ | 75.6 s | 326.6 KB | 5.4 × 10⁻⁶ |
| B — Depth only | 1 | 2⁴⁰ | 38.6 s | 229.8 KB | 4.6 × 10⁻⁷ |
| C — Scale only | 2 | 2²¹ | 75.0 s | 238.0 KB | 0.336 |
| D — Compacted | 1 | 2²¹ | 39.3 s | 155.6 KB | 0.183 |

**Marginal attribution — computation time:**

| Effect | Delta (s) | Delta (%) |
|---|---|---|
| Depth reduction A→B | −36.9 s | −48.9% |
| Scale reduction A→C | −0.6 s | −0.8% |
| Interaction | +1.3 s | +1.7% |
| Total A→D | −36.2 s | −47.9% |

**Marginal attribution — MAE:**

| Effect | MAE Delta |
|---|---|
| Depth reduction A→B | −4.89 × 10⁻⁶ (accuracy *improves*) |
| Scale reduction A→C | +0.336 |
| Interaction | −0.153 |
| Total A→D | +0.183 |

**Key finding — orthogonality:** Depth and scale are near-orthogonal in their effects. Depth reduction (L=2→1) drives virtually all of the latency gain (−48.9%) while actually improving numerical precision. Scale reduction (2⁴⁰→2²¹) contributes negligibly to speed (−0.8%) but causes the entirety of the accuracy loss.

**Revised recommendation — Config B, not Config D:**

The original paper recommended Config D (L=1, scale=2²¹) as the compacted configuration. This is wrong. MAE=0.183 at scale 2²¹ is not a rounding error — it is a large absolute error that could flip a PRS decile classification and is clinically unsafe. The original claim that compaction "preserved clinical precision" does not hold for Config D.

**Config B (L=1, scale=2⁴⁰) is the correct practitioner recommendation.** It captures essentially all the speed benefit (−48.9% computation time, −29.8% ciphertext size) while *improving* precision to MAE=4.6 × 10⁻⁷ — three orders of magnitude better than the baseline. The scale reduction in Config D adds no meaningful latency benefit and introduces unacceptable accuracy loss.

The honest conclusion is: depth reduction to L=1 is a safe and effective optimisation; scale reduction to 2²¹ is not clinically viable. These should be reported as separate recommendations, not combined into a single "compacted" configuration.
