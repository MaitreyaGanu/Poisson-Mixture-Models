<!-- PROJECT HEADER -->

<div align="center">

# 🧬 Hierarchical Bayesian Poisson Mixture Models  
### Empirical & Theoretical Analysis for Single-Cell Count Data

**Maitreya Sameer Ganu**  
*Indian Institute of Science Education and Research (IISER), Thiruvananthapuram*  

<i>Supervisor: Dr. Leelavati Narlikar , Indian Institute of Science Education and Research (IISER), Pune</i>

---

<img src="https://img.shields.io/badge/Language-R%20%7C%20Python-blue?style=for-the-badge&logo=r" />
<img src="https://img.shields.io/badge/Statistics-Bayesian%20Inference-orange?style=for-the-badge" />
<img src="https://img.shields.io/badge/Model-Poisson%20Mixtures-green?style=for-the-badge" />
<img src="https://img.shields.io/badge/Status-Completed-success?style=for-the-badge" />

---

[Abstract](#project-abstract) •
[Core Hypothesis](#gamma-variance-hypothesis) •
[Codes](#code-description) •
[Repository](#repository-structure) •
[Usage](#usage) •
[Notes](#notes)

</div>

---

## Project Abstract

This repository contains the **complete codebase** associated with a research project on **Hierarchical Bayesian Poisson Mixture Models (PMMs)**, developed as part of a Summer Internship 2025.

The project investigates why conventional clustering approaches fail on **discrete, sparse, and over-dispersed count data**, particularly in **single-cell genomics (scRNA-seq, scATAC-seq)**. We focus on Poisson mixture models with Gamma priors and show both empirically and theoretically that **prior variance plays a decisive role in cluster separability**.

All mathematical derivations, proofs, and detailed experimental discussion are provided in the **accompanying report (PDF)**. This repository is intentionally code-focused.

---

## Gamma Variance Hypothesis

A central finding of this work is the sensitivity of Poisson Mixture Models to the variance of the Gamma prior on rate parameters.

**Key Insight**

- **Low Gamma variance** forces Poisson rates across clusters to collapse, leading to posterior overlap and clustering failure.
- **High Gamma variance** enables distinct rate parameters and well-separated mixture components.

| Prior Regime | Empirical Outcome |
|-------------|------------------|
| Low Variance | Degenerate clusters, ARI ≈ 0 |
| High Variance | Clear separation, ARI → 1 |

These observations are consistently reproduced across synthetic data, PBMC scRNA-seq data, and scATAC-seq–derived gene activity matrices.

---

## Code Description

The repository consists of **stand-alone R scripts and notebooks**, each corresponding to a specific experiment, figure, or analysis reported in the thesis.

Major categories include:

- **Synthetic data generation** for Poisson mixtures
- **Gamma prior variance experiments**
- **EM vs Gibbs sampling comparisons**
- **PBMC scRNA-seq clustering**
- **scATAC-seq robustness analysis**
- **Visualization utilities (heatmaps, ARI plots, 3D λ surfaces)**

Each script is self-contained and reproducible.

---

## Repository Structure

```
│
├── 100_poisson_datasets.R
├── Synthetic_data_generation_for_poisson_mixture_models.R
├── Cluster_Separability_and_Lambda_Variance.R
├── alphaBeta3Dplot.R
├── alphaBeta3Dplot_copy.R
├── EM_and_Gibbs_Comparison.R
├── Flexmix_vs_PMM.R
├── ARIandLLrandomVSTrue.png
├── PBMC_different_clustering.R
├── PBMC.ipynb
├── pbmc_counts.csv
├── pbmc_labels.csv
├── ATACfinal.R
├── heatmaps.R
├── Mixture_of_Poissons.ipynb
```

**Note:**  
All theoretical derivations, proofs, and result interpretations are documented in the **report**, not duplicated here.

---

## Usage

### Requirements

- R ≥ 4.0  
- Recommended packages:
  - `flexmix`
  - `mclust`
  - `ggplot2`
  - `matrixStats`
  - `cluster`
  - `aricode`

Install common dependencies:
```r
install.packages(c(
  "flexmix", "mclust", "ggplot2",
  "matrixStats", "cluster", "aricode"
))
```
