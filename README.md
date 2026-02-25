# Reproducing Code and Data for: "Error control for microbiome mediator discovery: benchmark and remedy"

## 📖 Overview

Identifying microbial taxa that mediate the effect of exposures on health outcomes is central to understanding causal pathways. However, microbiome sequencing typically quantifies relative abundances (RA). Testing taxon-level mediation effects on the RA scale can yield distorted null behavior under compositionality, leading to excess false discoveries. 

To address this, our paper presents:
* **A Systematic Benchmark:** We benchmarked nine methods using extensive simulations grounded in experimentally recovered absolute abundance (AA) templates. We demonstrate that existing methods exhibit inflated false discovery rates (FDR) on RA count inputs, even in complete null settings.
* **The Proposed Method (CAMRA):** We introduce **CAMRA (Causal Absolute-abundance Mediation from Relative-Abundance data)**, a principled remedy that infers and tests AA-level mediation effects using standard RA inputs. CAMRA rigorously controls FDR, substantially improves discovery power and computational efficiency, and yields a smaller, more interpretable set of candidate mediators in real-data applications.

## 📦 Prerequisites

To run the scripts in this repository, you first need to install the `miMediation` R package, which contains the core `CAMRA` function. You can install the package directly from GitHub using `devtools`:

```R
# Install devtools if you haven't already
# install.packages("devtools")

# Install the miMediation package
devtools::install_github("tangzheng1/miMediation")
```
The details of this package can be seen in our package website https://github.com/tangzheng1/miMediation.

## 📊 Figure 1: Illustration

Figure 1 in our paper illustrates that the commonly used pseudo-AA construction (multiplying relative abundance by an independent, permuted microbial load) fails to reproduce key features of real absolute abundance (AA) data.

<img width="800" height="354" alt="Screenshot 2026-02-24 at 11 54 51 PM" src="https://github.com/user-attachments/assets/f371e67b-d8d3-4c17-99fd-1a8e22060267" />

The datasets (`GALAXY_mOTUs_v25.tsv`, `GALAXY_load.tsv`) used to generate Figure 1 and construct the AA templates were downloaded from Zenodo Record 14280080 (https://zenodo.org/records/14280080). You can run the script located in `figure1_illustration/figure1.R` to reproduce the plot.

## 🧐 Real Data Analysis

## 📚 Reference
Wang Q, Li Y, Peng Y, Zhang H, Tang, ZZ (2026). *Error control in microbiome mediator discovery: benchmark and remedy*. Submitted

