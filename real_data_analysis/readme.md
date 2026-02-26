# Real Data Analysis
This module contains the data processing pipeline and analysis scripts to investigate whether microbial taxa mediate the association between country of origin (United States vs. China) and body mass index (BMI).

## 🔧 Data Preprocessing
The raw metagenomic data is pulled from the `curatedMetagenomicData` R package. To ensure a rigorous causal evaluation and reduce confounding in cross-country comparisons, the study population is constructed in three stages:

1. **Initial Screening:** We extract stool metagenomes, excluding participants who were pregnant or reported current alcohol consumption/smoking. We focus on the two countries with the largest sample sizes: China and the United States.
2. **Propensity Score Matching:** We perform 1:1 propensity score matching (imposing exact matching on age, gender, and antibiotic use) to balance the covariates between the two countries. 
3. **Final Filtering:** We retain samples with complete BMI information and filter out rare taxa by requiring a non-zero proportion of at least 10% (i.e., zero proportion < 90%).

The final analytical dataset consists of 863 samples and 186 taxa, which can be found in the `data/` folder. You can also rerun `read_raw.R` to produce it.

## 🧐 Data Analysis
We apply both global mediation tests and taxon-level tests to the cleaned dataset. To reproduce the real data application tables and figures (including Figure 5 and Supplementary Table S4), you can run `data_analysis.R`.
