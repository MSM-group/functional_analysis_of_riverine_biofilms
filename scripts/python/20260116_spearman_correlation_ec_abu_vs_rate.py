import numpy as np
import pandas as pd
import json
from scipy import stats
from statsmodels.stats.multitest import multipletests


# Script to compute Spearman correlations between EC feature abundance and biotransformation rates

# === Feature Type Definition ===
feature_name = 'EC'

# === Load Abundances ===

df = pd.read_csv("data/20250604_EC_abundance_table.csv")

# === Load Biotransformation Rates ===

# Load the matrix of biotransformation rates (rows = samples, columns = compounds)
df_rates = pd.read_csv("data/20250212_bt_rates_matrix_without_AS.csv")

# Extract compound and EC lists for downstream correlation analysis
compound_list = df_rates.columns[2:]  # Skip the 'loc_batch' column
ec_list = df.columns[4:]

#Load mapping from accession (e.g., run ID) to sample name
with open("data/accession_lower_case_sample_name_dict.json", "r") as f:
    accession_sample_name_dict = json.load(f)

# === Align Sample Names Between DataFrames ===

# Map accessions in EC abundance table to sample names using loaded dictionary
df['sample_name'] = df['accession'].map(accession_sample_name_dict)

# Extract sample names from `loc_batch` strings in the biotransformation rates table
df_rates['sample_name'] = (
        df_rates['loc_batch'].str.split('_').str[0] + "_" +
        df_rates['loc_batch'].str.split('_').str[1] + "_" +
        df_rates['loc_batch'].str.split('_').str[-1].str[1]
)

# Ensure consistent ordering of samples in EC data to match df_rates
df['sample_name'] = pd.Categorical(
    df['sample_name'],
    categories=df_rates['sample_name'].tolist(),
    ordered=True
)
df = df.sort_values('sample_name')  # Sort for alignment

# === Spearman Correlation Analysis ===

# Initialize empty DataFrame to collect results
corr_df = pd.DataFrame(columns=['compound', 'EC', 'ro', 'p_value'])

# Loop through each compound and compute Spearman correlation with each EC
for compound in compound_list:

    #DEBUG:
    if compound == 'Saccharin':
        print(compound)

    corr_df_compound = pd.DataFrame(columns=['compound', 'EC', 'ro', 'p_value'])
    rate_vec = np.array(df_rates[compound], dtype=float)

    nan_ix = np.argwhere(np.isnan(rate_vec))
    nan_flag = False
    if nan_ix.any():
        print(f"{compound} has missing rate values")
        rate_vec = np.delete(rate_vec, nan_ix)
        nan_flag = True

    for index, ec in enumerate(ec_list):
        ec_vec = np.array(df[ec], dtype = float)
        if nan_flag:
           ec_vec = np.delete(ec_vec, nan_ix)

        res = stats.spearmanr(rate_vec, ec_vec)

        # Store correlation results
        corr_df_compound.loc[index, 'EC'] = ec
        corr_df_compound.loc[index, 'ro'] = res.statistic
        corr_df_compound.loc[index, 'p_value'] = res.pvalue

    corr_df_compound['compound'] = compound
    corr_df = pd.concat([corr_df, corr_df_compound])

# Save raw correlation results
corr_df.to_csv("output/20260120_spearman_corr_raw_results.csv", index=False)

# === Multiple Testing Correction and Filtering ===

# Adjust p-values for multiple testing using FDR (Benjamini-Hochberg)
corr_df["adjusted_p_value"] = multipletests(corr_df["p_value"], method="fdr_bh")[1]

# Filter for statistically significant positive correlations
corr_df_filtered = corr_df.loc[(corr_df['adjusted_p_value'] <= 0.05) & (corr_df['ro'] > 0)]

# Save filtered results
corr_df_filtered.to_csv("output/20260120_spearman_corr_filtered.csv", index=False)

print('done')