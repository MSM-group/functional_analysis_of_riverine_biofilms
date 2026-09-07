# Functional analysis of riverine biofilms
This repository contains the Python scripts and datasets required to reproduce the statistical analyses described in the following study:

> **Assessing Field-Scale Effects of Treated Wastewater on Riverine Biofilms: Community Composition, Functional Genes, and Micropollutant Biotransformation** 
> > *Authors: 
Martina Kalt; Victoria Poltorak; Eleonora Mastrorilli; Elia Ceppi; Yaochun Yu; Shinichi Sunagawa; Michael Zimmermann; Serina L. Robinson;
> Corresponding Author: Dr. Kathrin Fenner]* > **DOI:** [Insert DOI Link Here]


## Repository Structure

| Directory                  | Description                                              |
|:---------------------------|:---------------------------------------------------------|
| `scripts/python/`          | Python scripts for data processing and modeling.         |
| `scripts/R/`               | Python scripts for amplicon data analysis with FAPROTAX  |
| `scripts/bash/`            | bash script for FAPROTAX execution                       |
| `data/`                    | Input files (feature tables, EC abundances).             |
| `data/amplican_analysis`   | Input files for FAPROTAX analysis                        |
| `output/`                  | including correlation values, scripts output destination |
| `output/amplicon_analysis` | including correlation values, scripts output destination |

---

## Scripts

### 1. Feature Table Correlation
**File:** `scripts/python/20260424_mann_whitney_u_test_above_vs_below_lod.py`  
Mann-Whitney U-test for all the unique EC classes, comparing biotransformin samples vs non-biotransforming samples across samples, for individual compound.

### 2. Spearman Correlation (EC vs. Rates)
**File:** `scripts/python/20260428_random_forest_classification_analysis.py`  
Trains a Random Forest classifier model riverine parameters and EC features, to predict biotransforamtion. This script outputs a **Feature Importance** bar plot.

### 3. ASV based relative taxonomic abundances 
**File:** `scripts/R/20260826_asv_based_taxonomic_abundance_table.R`
Aggregates asv counts by taxonomic annotation 

### 4. Reformatting input tabel for FAPROTAX usage 
**File:** `scripts/R/20260826_prepare_input_table_faprotax.R`

### 5. FAPROTAX execution script 
**File:** `scripts/bash/20260826_execute_faprotax.sh`

### 6. Plotting estimated heterotroph/autotroph relative abundances 
**File:** `scripts/R/20260827_plotting_community_abundance_auto_hetero_ratios.R`

### 7. Plotting estimated heterotroph/autotroph relative abundances 
**File:** `scripts/R/20260827_plotting_faprotax_community_abundance_ratios.R`
