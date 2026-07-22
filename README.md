# Statistical Analysis of EC Abundances and Biotransformation Rates

This repository contains the Python scripts and datasets required to reproduce the statistical analyses described in the following study:

> **Assessing Field-Scale Effects of Treated Wastewater on Riverine Biofilms: Community Composition, Functional Genes, and Micropollutant Biotransformation** 
> > *Authors: 
Martina Kalt; Victoria Poltorak; Eleonora Mastrorilli; Elia Ceppi; Yaochun Yu; Shinichi Sunagawa; Michael Zimmermann; Serina L. Robinson;
> Corresponding Author: Dr. Kathrin Fenner]* > **DOI:** [Insert DOI Link Here]


## Repository Structure

| Directory | Description                                              |
| :--- |:---------------------------------------------------------|
| `scripts/python/` | Python scripts for data processing and modeling.         |
| `data/` | Input files (feature tables, EC abundances).             |
| `output/` | including correlation values, scripts output destination |

---

## Scripts

### 1. Feature Table Correlation
**File:** `scripts/python/20260424_mann_whitney_u_test_above_vs_below_lod.py`  
Mann-Whitney U-test for all the unique EC classes, comparing biotransformin samples vs non-biotransforming samples across samples, for individual compound.

### 2. Spearman Correlation (EC vs. Rates)
**File:** `scripts/python/20260428_random_forest_classification_analysis.py`  
Trains a Random Forest classifier model riverine parameters and EC features, to predict biotransforamtion. This script outputs a **Feature Importance** bar plot.



