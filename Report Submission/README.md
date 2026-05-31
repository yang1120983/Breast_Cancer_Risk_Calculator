# DATA3888 Biomed28: Breast Cancer Risk Prediction

Clinical decision-support project using gene expression (GEO microarray) data to (1) classify
tumour vs. normal breast tissue (binary) and (2) predict PAM50 molecular subtype (multiclass),
with a deployed Shiny risk calculator.

This folder contains everything needed to reproduce the report, the model outputs, and the app.

---

## 1. Folder Structure

```
Report Submission/
├── report.qmd                      # Main Quarto report (renders to HTML)
├── references.bib                  # Bibliography (BibTeX)
├── apa.csl                         # APA citation style
├── edit_report.css                 # Report styling
├── README.md                       # This file
│
├── model_scripts/                  # Analysis pipelines (run these to regenerate outputs)
│   ├── preprocessing.qmd           # GEO download, QC, normalisation, label inference
│   ├── binary_pipeline.R           # Binary classifier → writes RDS_Files/binary/
│   └── multiclass_pipeline.R       # Multiclass PAM50 classifier → writes RDS_Files/multiclass/
│
├── RDS_Files/                      # Generated model outputs consumed by report.qmd
│   ├── binary/                     #   (created by binary_pipeline.R)
│   ├── data/                       # Contains preprocessed data
│   └── multiclass/                 #   (created by multiclass_pipeline.R)
│   
│
├── images/                         # Contains images used for the report
│
└── ShinyApp/
    └── BreastCancerRiskApp/
        ├── app.R                   # Shiny risk calculator
        └── app_data/               # Trained models, results, and GEO platform cache

```

> **Note:** `RDS_Files/` is produced by the pipeline scripts. If it is missing, run the
> pipelines (Section 4) before rendering the report — `report.qmd` reads its tables and
> figures from `RDS_Files/binary/` and `RDS_Files/multiclass/`.

## 2. Prerequisites

- **R** ≥ 4.2
- **Quarto** ≥ 1.3 (for rendering `report.qmd`)
- **RStudio** recommended — the pipeline scripts auto-set their working directory via
  `rstudioapi`. If running outside RStudio, set the working directory manually (Section 4).

### R packages

```r
# CRAN
install.packages(c(
  "tidyverse", "ggplot2", "patchwork", "forcats", "knitr",
  "dplyr", "tidyr", "tibble", "data.table", "matrixStats",
  "glmnet", "e1071", "randomForest", "nnet", "MASS", "class",
  "xgboost", "caret", "janitor", "survival", "pbapply", "pheatmap",
  "shiny", "rsconnect", "DT", "rmarkdown"
))

# Bioconductor
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase", "limma", "sva"))
```

---

## 3. How to Run

### Step 1 — Generate model outputs (RDS files)

Run each pipeline once to populate `RDS_Files/`. In RStudio, open the script and source it
(the script sets its own working directory to `model_scripts/` and writes to `../RDS_Files/`):

```r
# In RStudio: open then Source, or
source("model_scripts/binary_pipeline.R")
source("model_scripts/multiclass_pipeline.R")
```

If running outside RStudio, set the working directory to `model_scripts/` first:

```r
setwd("path/to/Report Submission/model_scripts")
source("binary_pipeline.R")
source("multiclass_pipeline.R")
```

These scripts download datasets from GEO on first run (cached locally afterwards), perform LODO
cross-validation, gene-panel sweeps, and bootstrap uncertainty estimation, then save all tables
and figures used by the report into `RDS_Files/binary/` and `RDS_Files/multiclass/`.

> `preprocessing.qmd` documents the data download, QC, and label-inference steps. It is reference
> documentation for the biomedical pipeline; the two `.R` scripts are self-contained for
> reproducing the modelling outputs.

### Step 2 — Render the report

From the `Report Submission/` directory:

```bash
quarto render report.qmd
```

This produces `report.html`. The render reads precomputed results from `RDS_Files/`, so Step 1
must have completed first.

### Step 3 — Run the Shiny application (optional)

```r
shiny::runApp("ShinyApp/BreastCancerRiskApp")
```

Upload a gene-expression dataset (CSV, TXT, or GEO `.txt.gz`), select a patient, and the app
returns a cancer-risk score with a bootstrap confidence interval and, for medium/high-risk
patients, the predicted PAM50 subtype. The app loads its trained models from
`ShinyApp/BreastCancerRiskApp/app_data/`.

A live deployment is available at:
<https://data3888-group28.shinyapps.io/RiskCalculator-UI-DEMO/>

---