# Reproducibility Archive for Error control for microbiome mediator discovery: benchmark and remedy

This repository contains the data, code, and instructions necessary to reproduce all figures, tables, and numerical results presented in our paper (Wang Q, Li Y, Peng Y, Zhang H, Tang, ZZ (2026). Error control in microbiome mediator discovery: benchmark and remedy. Submitted). 

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
The details of this package is contained in https://github.com/tangzheng1/miMediation.

## 📊 Simulation Study

## 🧐 Real Data Analysis
