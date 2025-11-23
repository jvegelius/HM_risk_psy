# README: Hierarchical WOF + IDAS Example Repository

This repository provides a fully reproducible demonstration of the hierarchical variational-inference method for jointly estimating decision-making parameters from a Wheel-of-Fortune (WOF) task and latent psychological symptom factors from IDAS-like data.

It contains:

- **An example synthetic dataset** that mimics the structure of the empirical dataset used in the manuscript.
- **A fully annotated step-by-step estimation script**, showing how to reproduce the analytical pipeline.
- **Functions for the hierarchical model**, including data-generation and parameter-estimation routines.
- **Simulation code** used to evaluate accuracy, bias, and partial correlations in controlled settings.

The empirical dataset from the manuscript **cannot** be shared due to ethical and legal restrictions, but this repository provides:
- Synthetic data with identical structure.
- A script to generate new datasets.
- A complete example showing how to run the hierarchical model.

---

## Repository Contents

### **1. `Functions_for_hierarchical_WOF_psych.R`**
Contains all core functions required for:
- Data generation (`gen_data_fun()`)
- Computing choice probabilities
- Running the hierarchical estimation (`est_fun()`)
- Computing log-evidence, gradients, Hessians
- Extracting standard errors

These functions are sourced by all other scripts.

---

### **2. `Generate_data_N100.R`**
Generates a **fully synthetic example dataset** of:
- 100 subjects
- 80 WOF trials per subject
- 17 IDAS-like symptom scores

It outputs two CSV files:

- `example_wof_data.csv` — trial-level gambling data
- `example_idas_data.csv` — subject-level symptom data

These serve as inputs to the estimation script.

---

### **3. `step_by_step_guide.R`**
A fully annotated script showing exactly how to:

1. Load example data
2. Build subject-level trial lists
3. Specify initial values and parameter restrictions
4. Run the hierarchical estimation
5. Inspect and interpret the results

This script is intended as the main entry point for researchers who want to:
- Apply the hierarchical model to their own data
- Understand the required data formats
- Replicate the workflow used in the manuscript

---

### **4. `Simulation_N100_1000Reps.R`**
Contains all simulation procedures used to:
- Reproduce simulation results from the manuscript (for N = 100)
- Estimate gamma parameters repeatedly
- Compare hierarchical model (HM) vs. Independent-Estimation approach (IE)
- Compute bias, RMSE, partial correlations, and coverage

This script is more technical and intended for readers who want to verify statistical properties.

---

## How to Run the Example

### **1. Generate example data**
If the CSV files are not already included, run:
```r
source("Generate_data_N100.R")
```
This will produce:
- `example_wof_data.csv`
- `example_idas_data.csv`

### **2. Run the hierarchical estimation**
```r
source("step_by_step_guide.R")
```
This loads:
- Example data
- Model functions
- Parameter restrictions and initial values
- Runs `est_fun()`

Output includes:
- Estimated gamma matrix
- Estimated covariance matrices
- Subject-level xi estimates
- Standard errors
- Convergence info

---

## Example File Naming Conventions

- `Functions_for_hierarchical_WOF_psych.R` — all core model functions.
- `Generate_data_N100.R` — generate synthetic dataset of 100 subjects.
- `example_wof_data.csv` — synthetic WOF trial-level data.
- `example_idas_data.csv` — synthetic IDAS symptoms.
- `step_by_step_guide.R` — walk-through estimation example.
- `Simulation_N100.R` — simulation experiments.

---

## Citation
If you use this code in your own research, please cite the manuscript:

**“Hierarchical modeling of the relationships between internalizing psychiatric symptom factors and attitudes to risk and ambiguity using a variational inference approach”**

Manuscript Number: **YJMPS-D-25-00006**.

---

## Contact
For questions, please contact the corresponding author.

---

