# Monte Carlo Simulation of DIF Effects on Test Information Function

This repository contains the simulation code and design files for the manuscript:  

> *When DIF Goes Unmodeled: Assessing the Viability of Random Forest for Diagnostic Classification*  

The study investigates the impact of **Differential Item Functioning (DIF)** on test information and diagnostic classification, comparing **Item Response Theory (IRT)** and **Random Forest (RF)** approaches.

---

## Repository Structure

```
├── simulation.R        # Main simulation script (annotated)
├── allConditions.csv   # Design file: defines simulation conditions
├── results/            # Folder where simulation results are saved
└── README.md           # This file
```

---

## Requirements

The script was written in **R (≥ 4.2)**. Install required packages:

```r
install.packages(c(
  "mirt", "mclust", "pROC", "randomForest", 
  "dplyr", "MASS", "caret", "foreach", 
  "doParallel", "parallel"
))
```

---

## Running the Simulation

Run the main script:

```r
source("parallelSimCode.R")
```

- The script will iterate through all conditions specified in `allConditions.csv`.  
- Parallel processing is used automatically (outer parallel loop over replications, optional inner parallel loops for probability generation).  
- Results are written to `results/` as `.csv` files.  

---

## Design File (`allConditions.csv`)

Each row defines one simulation condition. Columns:

| Column          | Definition |
|-----------------|------------|
| `condition`     | Integer ID for the condition row |
| `refGroupSize`  | Number of simulees in the reference group |
| `focalGroupSize`| Number of simulees in the focal group |
| `nItems`        | Number of items |
| `prevRate`      | Prevalence rate of diagnosis (proportion in class = 1) |
| `difRate`       | Proportion of items containing DIF |
| `difType`       | Parameters manipulated: `"none"`, `"a"`, `"b"`, or `"both"` |
| `balanced`      | `"TRUE"` → DIF balanced (half easier, half harder); `"FALSE"` → unidirectional |
| `aStrength`     | Strength of discrimination DIF. For item *i* in DIF set *D*: <br> • increase: \(a_i^{foc} = a_i^{ref}(1 + aStrength)\) <br> • decrease: \(a_i^{foc} = a_i^{ref}(1 - aStrength)\) |
| `bStrength`     | Strength of threshold DIF. For each threshold \(b_{ik}\): <br> • increase: \(b_{ik}^{foc} = b_{ik}^{ref} + bStrength\) <br> • decrease: \(b_{ik}^{foc} = b_{ik}^{ref} - bStrength\) |

---

## Data Generation

- A single latent trait \(\theta_{\text{Diagnosis}} \sim N(0,1)\) is generated for each simulee.  
- Diagnosis classes are assigned using a prevalence-based threshold on \(\theta_{\text{Diagnosis}}\).  
- Item responses are generated with a **graded response model (GRM)** using group-specific parameters.  
- DIF is applied by adjusting **discrimination (a)** and/or **thresholds (b)** for focal group items.  

---

## Analyses Performed

For each replication:

1. **IRT model**: 1-factor GRM fit to training data with `mirt`.
   - Latent scores estimated with EAP.
   - Classification cut-point selected via **Youden’s J index**.

2. **Random Forest**: trained on the same data with prevalence-adjusted weights.
   - Cut-point on class-1 probability also selected via **Youden’s J index**.

3. **Evaluation metrics**:
   - Accuracy, Sensitivity, Specificity, Precision, NPV, F1-score.
   - Results reported for **full sample, reference group, focal group**.

---

## Output

For each condition:

- **`results/condition<id>FinalResults.csv`**  
  Contains classification metrics for all methods and groups.  

- **`results/condition<id>FinalVarImp.csv`**  
  Contains Random Forest variable importance values.  

---

## Reproducibility

- A **global seed** is set in the script (`TOP_LEVEL_SEED`).  
- Each replication derives its own seed to ensure reproducibility with variation across reps.  

---

## Analyses Performed

- This folder contains R code used to perform analysis on the simulation results

---

## Supplementary Materials (Optional)

- A copy of supplementary materials also submitted directly to the journal. 
