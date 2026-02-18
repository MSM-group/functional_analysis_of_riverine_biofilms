This repository contains the scripts and the input files for the statistical analysis of EC normalized relative abundances and the biotransforamtion rate constants, described in the paper (add title and doi)

In scripts/python you will find:

1. 20250610_feature_table_correlation_matrix.py
This script calculates the correlation values between all EC classes, across samples, for each compound.

2. 20260116_spearman_correlation_ec_abu_vs_rate.py
This script calculates and plots the correlation matrix of EC abundances and river parameters

3. 202560218_random_forest_regression_analysis.py
This script trains a random forest model on the chosen set of parameters. It plots the feature importance barplot

in the data folder you will find the input files for running the scripts above and the output folder contains the calculated correlation values.

