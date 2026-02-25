# Reproducibility Archive for Error control for microbiome mediator discovery: benchmark and remedy

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

Figure 1 in our paper are used to illustrate thatthe commonly used pseudo-AA construction (multiplying relative abundance by an independent, permuted microbial load) fails to reproduce key features of real absolute abundance (AA) data.

The datasets (`GALAXY_mOTUs_v25.tsv`, `GALAXY_load.tsv`) used to generate Figure 1 and construct the AA templates were downloaded from Zenodo Record 14280080 (https://zenodo.org/records/14280080). You can run the script located in `figure1_illustration/figure1.R` to reproduce the plot.

## 📉 Simulation Study

### 🎯 Tasks & Methods
We benchmarked nine existing microbiome mediation methods across two distinct inferential tasks, alongside our proposed method, CAMRA:

1. **Taxon-Level Mediator Discovery**

- **Goal**: Pinpoint individual mediating taxa and evaluate empirical False Discovery Rate (FDR) and discovery power.
- **Methods Evaluated**: `CAMRA`, `CRAmed`, `LDM-med`, `MarZIC`, `microHIMA`, and `multimedia`.

2. **Global Community Mediation Testing**

- **Goal**: Assess whether the microbiome community as a whole mediates the pathway, focusing on Type I error rate under null settings.
- **Methods Evaluated**: `CMM`, `LDM-med`, `MedTest`, `MODIMA`, and `PERMANOVA-med`.

### ⚙️ Simulation Settings

We employed a template-based resampling approach grounded in the real Absolute Abundance (AA) data from the GALAXY/MicrobLiver cohort (the template data `GALAXYMicrobLiver_study.RData` can be found in the `data/` folder). 

To systematically evaluate the methods, we generated 500 simulated datasets per scenario and varied the following key configurations:

- **Data Dimensions**: Sample sizes (n=200/400/800) and microbiome feature sizes (p=200/400).
- **Effect Direction Balance**: Balancing factor (d=0.5/0.9). When both exposure and outcome associations were present, we evaluated scenarios with Balanced +/- signs (50% positive/negative) and Dominant + signs (90% positive/10% negative) to assess the impact of compositional distortion.
  
We designed carefully controlled causal scenarios by assigning exposure- and outcome-associations to specific subsets of taxa (fixed at 10 taxa per path when nonempty):

1. **Type I Error Evaluation (Null Scenarios)**

We designed four distinct null settings to test if methods can suppress false positives when no true mediators exist:

- **Complete Null**: No taxa are associated with either the exposure or the outcome.
- **Exposure-only**: Taxa are associated only with the exposure.
- **Outcome-only**: Taxa are associated only with the outcome.
- **Disjoint Null**: Both paths are present (10 taxa each), but they occur in completely different, non-overlapping taxa sets. 

2. **FDR and Power Evaluation (True Mediation)**

- **True Mediation**: A specific subset of taxa acts as valid mediators by having nonzero AA association effects on both paths (i.e., overlap exists between the two 10-taxa sets). We varied the number of true mediators (**overlap size** num2 = 3/5/7/9).

### 🏃 Reproduction Steps

The R scripts (`taxon_level_simulation.R` and `global_test_simulation.R`) used to run the simulations on the UW-Madison CHTC cluster are located in `simulations/code/` folder. For each simulation, you can get a RDS file named as "template_%s_n_%d_p%d_d%s_num1A_%d_num1B_%d_num2_%d_seed_%d.rds". 

The RDS files from taxon level simulations contain the index of true mediators, q-value matrix and running time for different methods. The RDS files from global test simulations contain the global p-values and running time for different methods. 

After having all RDS files, you can use the R scripts (`taxon_level_summary.R` and `global_test_summary.R`) located in `simulations/code/` folder to summarize the data, making plots and summary tables. 

## 🧐 Real Data Analysis
