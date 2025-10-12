# 🧬 Defining the Limits of Pre-Transplant Risk Prediction in AML: Evidence from Machine Learning and Regression Models.

> Reproducible pipeline to compare and externally validate AML post-transplant survival models  
> *(Cox with time-varying covariates, Elastic-Net Cox, HCT-CI/ELN baselines, and Random Survival Forests)*  
> Includes discrimination, calibration, accuracy, clinical utility, reclassification, and individual-level heterogeneity.

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue?logo=r)
![Reproducible](https://img.shields.io/badge/Reproducible-Yes-brightgreen?logo=github)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Mac%20%7C%20Linux-lightgrey)

---

## ✨ Highlights

- 🔄 **End-to-end R workflow** — from raw scores to figures, tables, and Excel report  
- 🧠 **Models:** Cox-TVC, Cox-EN, HCT-CI, ELN, and **RSF** (if RSF workspace present)  
- 📊 **Metrics:** Apparent & optimism-corrected **C-index**, **time-dependent AUC/ROC (12/24/36m)**, **IPCW Brier & IPA**, **IPCW DCA**, **IPCW NRI (24/36m)**  
- ⚙️ **Calibration:** Pre- and post-recalibration on the logit scale  
- 🌍 **External validation:** Bootstrapped CIs and publication-ready plots  
- 🧩 **Individual heterogeneity:** Within-strata dispersion & prediction spectra  
- 🧾 **Reproducibility:** Fixed seed (`147201`) and fully scripted outputs

---

## Table of Contents

1. [📁 Repository Layout](#repository-layout)  
2. [⚙️ Prerequisites](#prerequisites)  
3. [📄 Data Contracts](#data-contracts)  
4. [🚀 Quick Start](#quick-start)  
5. [📊 Outputs](#outputs)  
6. [🧮 Methods](#methods)  
7. [🧱 Reproducibility & Logging](#reproducibility--logging)  
8. [🧰 Troubleshooting](#troubleshooting)

---

## Repository Layout

```
R code/
  ├─ New_Complete_Analysis_Code_4.R          # Main internal model comparison pipeline
  ├─ RSF_code.R                               # Variant focused on RSF-enabled runs
  ├─ external_validation.R                    # External cohort validation (C-index + plot)
  ├─ Individual_patient_heterogeneity.R       # Figure 2: within-strata heterogeneity
  └─ Stata_code.R                             # Legacy name; R code duplicate of main pipeline
```

**Expected data/workspaces (not included):**
```
01_Raw_Data/
  ├─ internal_scores.xlsx      # Internal cohort scores/labels for the main pipeline
  ├─ rsf_workspace.RData       # Optional: RSF fitted objects (list 'final_models')
  └─ External cohort.xlsx      # External cohort for out-of-sample validation
```

> 💡 *Tip:* Replace absolute Windows paths with relative project roots for cross-OS portability.

---

## Prerequisites

- 🧩 **R ≥ 4.2** recommended
- 📦 Core packages:
  `readxl`, `survival`, `survcomp`, `timeROC`, `pec`, `riskRegression`, `ggplot2`,
  `gridExtra`, `dplyr`, `tidyr`, `boot`, `reshape2`, `openxlsx`, `cowplot`  
  For external validation when using RSF predictions: `randomForestSRC`.

- 🌲 Optional but encouraged for reproducibility:
```r
install.packages("renv")
renv::init()
renv::install(c(
  "readxl","survival","survcomp","timeROC","pec","riskRegression","ggplot2",
  "gridExtra","dplyr","tidyr","boot","reshape2","openxlsx","cowplot","randomForestSRC"
))
renv::snapshot()
```

---

## Data Contracts

### 📘 `internal_scores.xlsx` (for the main pipeline)
- **Outcomes**: `Time` (months), `Status` (1=event, 0=censored)
- **Model scores**: `Cox_TVC`, `Cox_EN`; RSF is read from `rsf_workspace.RData`
- **Baselines**: `HCTCI` (rounded to integer categories), `ELN` (factor)

> The scripts cast `HCTCI` and `ELN` to factors:
> ```r
> scores_data$HCTCI <- as.factor(round(HCTCI))
> scores_data$ELN   <- as.factor(ELN)
> ```

### 🧩 `rsf_workspace.RData` (optional)
- Must contain list **`final_models`** (from `randomForestSRC`) exposing:
  `predicted.oob`, `err.rate`, `ntree`, `survival.oob`, `time.interest`

### 🌍 `External cohort.xlsx` (for `external_validation.R`)
- Predictors used to reconstruct linear predictors / categories:
  `age, sex, race, kps, hctci, amlgp, WC, elngroup, time_to_CR, mrd, donor, cond, graftype, lymphodepletion, Time, Status`
- Optional coefficient files: `cox_tvc_coefficients.xlsx`, `cox_en_coefficients.xlsx`
- Optional RSF workspace: `rsf_workspace.RData` (to generate RSF scores on the external set)

---

## Quick Start

### 🧭 1) Configure paths

Edit the *very top* of each script; replace hard‑coded Windows paths with your project root. Example:

```r
# At the top of each script:
proj <- normalizePath(".")
setwd(proj)
dir.create(file.path(proj, "model_comparison", "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(proj, "model_comparison", "Tables"),  recursive = TRUE, showWarnings = FALSE)
```

Place data in `01_Raw_Data/` as shown in **Data Contracts**.

### ⚡ 2) Run analyses

```r
# Full internal comparison pipeline (figures + tables + Excel report)
source("R code/New_Complete_Analysis_Code_4.R")

# RSF-focused variant (optional)
source("R code/RSF_code.R")

# External validation on a held-out cohort (C-index + plot + log + Excel)
source("R code/external_validation.R")

# Individual-level outcome heterogeneity (two figures + summaries)
source("R code/Individual_patient_heterogeneity.R")
```

> Seed is fixed at **`147201`** inside scripts for reproducibility.

---

## Outputs

### 📈 From `New_Complete_Analysis_Code_4.R` (main)

**Tables →** `model_comparison/Tables/`
- `cindex_results.csv` — apparent C-index for HCTCI/ELN, optimism‑corrected for Cox models; RSF OOB if available  
- `auc_results.csv`, `auc_results_lower.csv`, `auc_results_upper.csv` — time‑dependent AUCs (12/24/36m) with bootstrap CIs  
- `brier_results*.csv` — IPCW Brier & **IPA** at 12/24/36m with CIs  
- `dca_results_24m.csv`, `dca_results_36m.csv` — net benefit vs threshold  
- `nri_results.csv` — **IPCW NRI** comparing Cox_TVC/Cox_EN/RSF vs HCTCI/ELN (24/36m)  
- **Excel report**: `Complete_Model_Comparison_Results.xlsx` (Summary, C-index, AUC, Brier, NRI, Calibration)

**Figures →** `model_comparison/Figures/`
- `ROC_12months.png`, `ROC_24months.png`, `ROC_36months.png`  
- `PreCalibration_24months.png`, `PostCalibration_24months.png`  
- `DCA_24months.png`, `DCA_36months.png`

**Workspace & logs**
- `model_comparison/complete_analysis.RData`  
- `model_comparison/Complete_Analysis_Log.txt`

### 🌍 From `external_validation.R`

**Artifacts →** `03_Results/External_Validation_of_Models/`
- `cindex_results.csv`, `results.xlsx`, `results_summary.txt`, `cindex_plot.png`  
- `validation_workspace.RData`  
- Run log: `validation_log_YYYYMMDD.txt`

### 🧠 From `Individual_patient_heterogeneity.R`

**Artifacts →** `Individual_Patient_Analysis/`
- `Figure_2A_Individual_Outcome_Heterogeneity.emf` (density by risk strata & model)  
- `Figure_2B_Individual_Outcome_Heterogeneity.emf` (probability vs time scatter, log‑scaled)  
- `Figure_2_Combined_Individual_Outcome_Heterogeneity.emf`  
- `Individual_Outcome_Heterogeneity_Summary.csv`  
- `Individual_Outcome_Heterogeneity_KM_Results.csv`  
- `Individual_Outcome_Heterogeneity_Log.txt`

---

## Methods

### 🎯 Discrimination
- **C‑index** via `survcomp::concordance.index` (Noether method).  
  Cox models report *pre‑computed optimism‑corrected* values embedded in the script.  
- **Time‑dependent AUC/ROC** at 12/24/36 months via `timeROC::timeROC`, with **bootstrap CIs**.
- Score direction is checked and flipped if needed so higher scores imply higher risk.

### ⚖️ Calibration
- Recalibration on the **logit** scale: fit `glm(outcome ~ logit(pred))` then inverse‑logit to obtain **recalibrated probabilities**.  
- Plots aggregate by deciles (or unique levels for categorical baselines).

### 🎯 Accuracy
- **IPCW Brier** with **IPA** improvement over null risk `p0(1-p0)`; censoring weights from KM of censoring.

### 💉 Clinical Utility
- **IPCW Decision Curve Analysis** (net benefit) across thresholds 0.10–0.45 at 24/36m; includes **Treat All** and **Treat None** strategies.

### 🔄 Reclassification
- **Time‑to‑event NRI with IPCW** comparing (Cox_TVC, Cox_EN, RSF) vs (HCTCI, ELN) at 24/36 months; bootstrap CI and p‑values.

### 🌍 External Validation
- Reconstruct linear predictors from coefficient tables; RSF predictions via loaded forest objects; compute **C‑index** with 1000‑bootstrap CI; export CI bar plot.

---

## Reproducibility & Logging

- 🔢 **Seed**: `set.seed(147201)` used throughout (bootstrap CIs, resampling).  
- 📋 **Determinism**: RSF OOB metrics reflect saved forest objects; to re‑fit forests, fix RNG seeds during training.  
- 🧾 **Comprehensive logs** accompany every run; inspect them if metrics look unusual.

---

## Troubleshooting

- ❌ **Paths**: Replace hard‑coded `C:/Users/...` with project‑relative paths (`normalizePath(".")` + `file.path(...)`).  
- 🧩 **RSF missing**: Ensure `01_Raw_Data/rsf_workspace.RData` exists and contains list `final_models` with required members.  
- ⚠️ **Factor levels**: Inputs must be clean; scripts cast `HCTCI` and `ELN` to factors.  
- 🧩 **NA AUC**: `timeROC` can return `NA` with too few events by the evaluation time; see the run log for details.  
- 🧮 **Censoring weights**: IPCW assumes non‑informative censoring; extreme censoring can destabilize weights—check logs.
