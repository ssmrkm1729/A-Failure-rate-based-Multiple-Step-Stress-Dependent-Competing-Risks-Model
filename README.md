# Supplementary Code for "A Failure-rate based Multiple Step-Stress Dependent Competing Risks Model with Applications to Aerospace Electrical Connectors Data"

This repository contains the **R scripts, data, and sample results** used as supplementary material for the paper:

**"A Failure-rate based Multiple Step-Stress Dependent Competing Risks Model with Applications to Aerospace Electrical Connectors Data"**


## Authors

- **Sarat Sindhu Mukhopadhyay**  
  SQC and OR Unit, Indian Statistical Institute, Bangalore, Karnataka, India  
  Email: ssmrkm1729@gmail.com  

- **Ayan Pal**  
  Department of Statistics, The University of Burdwan, West Bengal, India  
  Email: ayan.pal33@gmail.com  

- **Gijo E.V.**   
  SQC and OR Unit, Indian Statistical Institute, Bangalore, Karnataka, India  
  Email: gijoev@gmail.com  



## Repository Structure

- `code/` : Contains 4 independent R scripts:
  - `MOBBXII.R`
  - `MOBCh.R`
  - `MOBE.R`
  - `MOBW.R`
  - `dataset_ssm.xlsx` : Required Excel dataset for all scripts
- `results/` : Contains outputs corresponding to the scripts.

## How to Use

1. Open any `.R` script in **RStudio** or any R environment.  
2. Ensure `dataset_ssm.xlsx` is located in the `code/` folder.  
3. Run the script to reproduce analyses (results will be displayed or saved if the script is designed to do so).

## Required R Packages

The scripts require the following R packages:

- `numDeriv`  
- `openxlsx`  
- `english`  
- `SMPracticals`  
- `readxl`  
- `plotly`  
- `RColorBrewer`  
- `survival`  
- `processx`  
- `orca`  
- `BoutrosLab.plotting.general`  

