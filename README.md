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
**File:** `20250610_feature_table_correlation_matrix.py`  
Calculates the correlation values between all EC classes, across samples, for each individual compound.

### 2. Spearman Correlation (EC vs. Rates)
**File:** `20260116_spearman_correlation_ec_abu_vs_rate.py`  
Computes and visualizes the correlation matrix comparing EC abundances against specific river parameters.

### 3. Random Forest Regression
**File:** `202560218_random_forest_regression_analysis.py`  
Trains a Random Forest model on the chosen set of parameters. This script outputs a **Feature Importance** bar plot.

