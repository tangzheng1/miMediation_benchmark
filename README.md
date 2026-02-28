# Error control for microbiome mediator discovery: benchmark and remedy.

## Overview

Identifying microbial taxa that mediate the effect of exposures on health outcomes is central to understanding causal pathways. However, microbiome sequencing typically quantifies relative abundances (RA). Testing taxon-level mediation effects on the RA scale can yield distorted null behavior under compositionality, leading to excess false discoveries. 

To address this, our paper presents:
* **A Systematic Benchmark:** We benchmarked nine methods using extensive simulations grounded in experimentally recovered absolute abundance (AA) templates. We demonstrate that existing methods exhibit inflated false discovery rates (FDR) on RA count inputs, even in complete null settings.
* **The Proposed Method (CAMRA):** We introduce **CAMRA (Causal Absolute-abundance Mediation from Relative-Abundance data)**, a principled remedy that infers and tests AA-level mediation effects using standard RA inputs. CAMRA rigorously controls FDR, substantially improves discovery power and computational efficiency, and yields a smaller, more interpretable set of candidate mediators in real-data applications.

## Structure

This repository contains the data, code, and instructions necessary to reproduce all figures, tables, and numerical results presented in our paper “Error control in microbiome mediator discovery: benchmark and remedy” by Qiyu Wang, Yiluan Li, Yunfei Peng, and Zheng-Zheng Tang.

The subfolder `data/` contains:
- The template dataset `GALAXYMicrobLiver_study.RData`, which is used to generate simulated datasets under various experimental settings.
- The processed dataset `real_data.RData`, which is used for the real data analysis presented in our paper.

The subfolder `simulations/` contains the simulation R codes and submission scripts for high throught computing, and the R scripts for aggregating and summarizing the simulation outputs.

The subfolder `real_data_analysis\` contains 
- `read_data.R`: R script for cleaning and preprocessing the raw microbiome dataset.
- `data_analysis.R`: R script for performing the mediation analysis on the cleaned dataset.

Executing these scripts can reproduce the results described in our paper under "Results – Gut microbiome mediation of between-country differences in BMI".

The subfolder `figures and tables` contains the R scripts for reproducing all the figures and tables in our paper and supplementary information, including Figure 1~5, Figure S1~S5, and Table S1~S4. 

## Reference
Wang Q, Li Y, Peng Y, Tang, ZZ (2026). *Error control in microbiome mediator discovery: benchmark and remedy*. Submitted

