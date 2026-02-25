# Simulation Study

## 🎯 Tasks & Methods
We benchmarked nine existing microbiome mediation methods across two distinct inferential tasks, alongside our proposed method, CAMRA:

1. **Taxon-Level Mediator Discovery**

- **Goal**: Pinpoint individual mediating taxa and evaluate empirical False Discovery Rate (FDR) and discovery power.
- **Methods Evaluated**: `CAMRA`, `CRAmed`, `LDM-med`, `MarZIC`, `microHIMA`, and `multimedia`.

2. **Global Community Mediation Testing**

- **Goal**: Assess whether the microbiome community as a whole mediates the pathway, focusing on Type I error rate under null settings.
- **Methods Evaluated**: `CMM`, `LDM-med`, `MedTest`, `MODIMA`, and `PERMANOVA-med`.

## ⚙️ Simulation Settings

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

## 🏃 Reproduction Steps

The R scripts (`taxon_level_simulation.R` and `global_test_simulation.R`) used to run the simulations are located in `simulations/code/` folder. 

```{r}
## Example for taxon level test
kk <- runone_simulation(
  n = 200, 
  p = 200, 
  num1_A = 10, 
  num1_B = 10, 
  num2 = 7,
  beta_treat = log(5), 
  beta_outcome = 1, 
  d = 0.5,
  template = "GALAXYMicrobLiver_study",
  template_dir = ".",
  save_dir = ".",
  seed = 1
)
kk
```

The simulations in our paper were executed through High Throughput Computing (HTC) environments, so we also provide the corresponding CHTC submission scripts in the `simulations/CHTC/` folder.

For each simulation, you can get a RDS file named as "template_%s_n_%d_p%d_d%s_num1A_%d_num1B_%d_num2_%d_seed_%d.rds". The RDS files from taxon level simulations contain the index of true mediators, q-value matrix and running time for different methods. The RDS files from global test simulations contain the global p-values and running time for different methods. 

After having all RDS files, you can use the R scripts (`taxon_level_summary.R` and `global_test_summary.R`) located in `simulations/code/` folder to summarize the data, making plots and summary tables. 
