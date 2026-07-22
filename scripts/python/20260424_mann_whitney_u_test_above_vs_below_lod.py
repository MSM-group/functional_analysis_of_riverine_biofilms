import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json

def get_the_samples_for_sites_locations(df, annotated_df, include_location, include_sites):

    annotated_df['location'] = annotated_df['sample_name'].str.split('_').str[1]
    annotated_df['site'] = annotated_df['sample_name'].str.split('_').str[2]
    annotated_df.loc[annotated_df['site'].str[-1].str.isdigit(), 'site'] = annotated_df.loc[annotated_df['site'].str[
        -1].str.isdigit(), 'site'].str[:-1]
    include_location = ['Kol', 'Bir', 'Ehr', 'Alt', 'Air', 'Oli']
    include_sites = ['dw', 'up']  # , 'Eff', 'In', 'AS']
    annotated_df = annotated_df.loc[annotated_df['site'].isin(include_sites)]
    annotated_df = annotated_df.loc[annotated_df['location'].isin(include_location)]
    sample_location_dict = dict(zip(annotated_df['accession'], annotated_df['location']))
    sample_site_dict = dict(zip(annotated_df['accession'], annotated_df['site']))

    # parsing columns and keeping only the relevant by location and site
    df.drop(columns='length', inplace=True)
    column_list = list(df.columns)
    sample_list = [column.split('_')[1] for column in column_list]
    df.columns = sample_list
    df = df[list(sample_location_dict.keys())]

    return df, annotated_df, sample_site_dict, sample_location_dict


def label_group(annotated_df, label, df):
    """
    This function labels the group of interest in the gene abundance table

    :param annotated_df (df): the metadata df that contains the links between sample names and accessions
    :param label (list): list of sample names (not accessions) that belong to a group of interest
    :param df (df): abundances df
    :return: df (df): labeled abundances df
    """

    accession_list_label = annotated_df.loc[annotated_df['sample_name'].isin(label), 'accession'].to_list()
    df['label'] = 0
    df.loc[df['accession'].isin(accession_list_label), 'label'] = 1

    return df

def mann_whitney_u_test_features_labels(df, output_filename, feature_name):
    """
    performs individual test for every ec in the gene catalog abundance table

    :param df (df): abundances df
    :param feature_name (str) : analysis name/ label name (ec/ko)
    :param output_filename (str)
    :return: significant_features (df): table with significant ecs and their p values
    """

    #from chatgpt:

    # Separate ec columns and metadata
    feature_columns = df.columns[1:-2]  # Assuming first column is accession and last two are sample name and label

    # Define groups
    group_1 = df[df["label"] == 1]
    group_2 = df[df["label"] == 0]

    # 🔹 Step 2: Perform Wilcoxon Rank-Sum Test (Mann-Whitney U) on Each KO
    # ---------------------------------------------------------------------
    results = []
    for feature in feature_columns:

        values_1 = group_1[feature].dropna()  # Drop NaNs if present
        values_2 = group_2[feature].dropna()

        if len(values_1) > 1 and len(values_2) > 1:  # Ensure enough samples
            stat, p_value = stats.mannwhitneyu(values_1, values_2, alternative="two-sided")
        else:
            p_value = 1  # If not enough data, set p-value to 1 (no significance)

        results.append((feature, p_value))

    # Convert results to DataFrame
    results_df = pd.DataFrame(results, columns=[feature_name, "p_value"])

    # 🔹 Step 3: Adjust for Multiple Testing (FDR Correction)
    # -------------------------------------------------------
    results_df["adjusted_p_value"] = multipletests(results_df["p_value"], method="fdr_bh")[1]

    # 🔹 Step 4: Filter Significant KOs (FDR ≤ 0.05)
    # ----------------------------------------------
    significant_features = results_df[results_df["adjusted_p_value"] <= 0.05]

    # 🔹 Step 5: Save or Print Results
    # --------------------------------
    # significant_features.to_csv(f"output/{output_filename}.csv", index=False)
    print(significant_features)

    return significant_features

def plot_boxplots(significant_features, df, feature_name, top_n, plot_title, path_fig):
    """
    Plots boxplots for the top significant  features with adjusted p-values in the x-tick labels.

    Parameters:
        - significant_features: DataFrame containing 'ec' and 'adjusted_p_value' columns.
        - df: DataFrame with EC abundances and a 'label' column (grouping variable).
        - feature_name: String, name of the feature being analyzed (e.g., "ec"/"ko").
        - top_n: Number of top features to plot.
        - plot_title: Title for the plot.
        - path_fig: Path to save the figure.
    """

    # Select top significant ECs based on adjusted p-value
    top_features = significant_features.sort_values("adjusted_p_value").head(top_n)[[feature_name, "adjusted_p_value"]]


    # Log transform abundance values (excluding 'label' and 'sample_name' column)
    df[df.columns[1:-1]] = np.log2(df[df.columns[1:-1]] + 1)  # Adding 1 to avoid log(0)

    # Melt data for seaborn
    melted_data = df.melt(id_vars=["label"], value_vars=top_features[feature_name],
                          var_name=feature_name, value_name="cummulative_relative_abundance")

    # Create x-tick labels with adjusted p-values
    xtick_labels = [f"{feature} (p={pval:.3f})" for feature, pval in zip(top_features[feature_name], top_features["adjusted_p_value"])]

    # Plot
    plt.figure(figsize=(20, 12))
    sns.boxplot(x=feature_name, y="cummulative_relative_abundance", hue="label", data=melted_data)

    plt.xlabel(f"{feature_name.capitalize()} Labels", fontsize=14)
    plt.ylabel(f"log2({feature_name} Cumulative Relative Abundance +1)", fontsize=14)
    plt.title(f"Top {top_n} Significant {feature_name} Labels {plot_title}", fontsize=14)
    plt.legend()

    # Set xticks with labels including p-values
    plt.xticks(ticks=range(len(top_features)), labels=xtick_labels, rotation=35, ha="right", fontsize=14)

    plt.tight_layout()
    plt.savefig(path_fig)

    return


def get_significant_features_and_abundance_df(abundance_df, feature_name, feature_list, sample_site_dict,
                                              sample_location_dict):
    """

    :param abundance_df: This df contains the abundances of the features in the gene catalog
    :param significance_df: This df contains the significant features and their p-values
    :return: melted_df: This df contains melted info: the abundances of the significant features in different samples
    """

    # Melt data for seaborn
    melted_df = abundance_df.melt(id_vars=["label", 'accession'], value_vars=feature_list,
                                  var_name=feature_name, value_name="cummulative_relative_abundance")

    melted_df['sample_site'] = melted_df['accession'].map(lambda x: sample_site_dict[x])
    melted_df['sample_location'] = melted_df['accession'].map(lambda x: sample_location_dict[x])

    return melted_df


def get_EC_significance_df_with_abundances_info(significant_abu_df, significant_features_df):
    """
    This df provides the info about the enrichment/deficiency of the significant feature when comparing the
    two labeled groups

    :param significant_abu_df: This df contains the abundances of the features in the gene catalog
    :param significant_features_df: This df contains the significant features and their p-values
    :return: significant_features_df with the additional information about the mean, std dev,
    median in the two compared groups
    """

    g_1 = significant_abu_df[significant_abu_df["label"] == 1].groupby(by=feature_name)
    means_1 = g_1['cummulative_relative_abundance'].apply(lambda x: np.mean(x))
    std_dev_1 = g_1['cummulative_relative_abundance'].apply(lambda x: np.std(x))
    median_1 = g_1['cummulative_relative_abundance'].apply(lambda x: np.median(x))

    g_0 = significant_abu_df[significant_abu_df["label"] == 0].groupby(by=feature_name)
    means_0 = g_0['cummulative_relative_abundance'].apply(lambda x: np.mean(x))
    std_dev_0 = g_0['cummulative_relative_abundance'].apply(lambda x: np.std(x))
    median_0 = g_0['cummulative_relative_abundance'].apply(lambda x: np.median(x))

    significant_features_df['mean_1'] = means_1[significant_features_df[feature_name]].values
    significant_features_df['std_1'] = std_dev_1[significant_features_df[feature_name]].values
    significant_features_df['median_1'] = median_1[significant_features_df[feature_name]].values

    significant_features_df['mean_0'] = means_0[significant_features_df[feature_name]].values
    significant_features_df['std_0'] = std_dev_0[significant_features_df[feature_name]].values
    significant_features_df['median_0'] = median_0[significant_features_df[feature_name]].values

    significant_features_df['delta_mean'] = significant_features_df['mean_1'] - significant_features_df['mean_0']
    significant_features_df['delta_median'] = significant_features_df['median_1'] - significant_features_df['median_0']

    return significant_features_df

def check_EC(query_ec_list, target_ec_list):
    for ec in query_ec_list:
        num_list = ec.split('.')
        third_level_ec = f"{num_list[0]}.{num_list[1]}.{num_list[2]}"
        if third_level_ec in target_ec_list:
            return True
    return False

def reduce_EC(query_ec_list, target_ec_list):

    for ec in query_ec_list:
        num_list = ec.split('.')
        third_level_ec = f"{num_list[0]}.{num_list[1]}.{num_list[2]}"
        if third_level_ec in target_ec_list:
            return third_level_ec


def group_ec_in_significant_feature_df(significant_features_df, feature_list_path):

    feature_list = []
    with open(feature_list_path) as f:
        for line in f:
            feature_list.append(line.strip('\n'))

    significant_features_df['EC'] = significant_features_df['EC_list'].apply(lambda x: reduce_EC(x, feature_list))

def significant_feature_df_for_list_of_features(feature_name, feature_list, significant_features_df):

    significant_features_df[f"{feature_name}_list"] = significant_features_df[feature_name].str.split(',')
    if feature_name == "EC":
        significant_features_df['EC_target'] = significant_features_df['EC_list'].apply(lambda x: reduce_EC(x, feature_list))
    significant_features_df_in_list = significant_features_df.dropna(subset = 'EC_target')

    return significant_features_df_in_list



if __name__=="__main__":

    compound = 'Saccharin' #or "Cyclamate"
    date = 20260424
    target_EC = ['1.14.12'] #3.1.6 for Cyclamate
    # Load your dataset using Dask (replace with your file path)
    # Assuming 'df' is a Dask DataFrame where rows are genes and columns are samples
    df = pd.read_parquet("output/20260325_EC_annotated_counts.parquet")  # Adjust the path to your data file

    annotated_df = pd.read_csv("data/94_metadata_with_biosamples_final.csv")
    include_location = ['Kol', 'Bir', 'Ehr', 'Alt', 'Air', 'Oli']
    include_sites = ['dw', 'up']  # , 'Eff', 'In', 'AS']

    #get dfs and dicts only for the relevant samples
    df, annotated_df, sample_site_dict, sample_location_dict = get_the_samples_for_sites_locations(df, annotated_df, include_location, include_sites)

    # in metagenomics: ['alt_dw1_r3', 'alt_up2_r3'] but should be 'alt_dw1_r2' and 'alt_up2_r2'
    annotated_df['sample_name'] = annotated_df['sample_name'].str.split('_', n=1).str[1]
    df = df.T
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'accession'}, inplace=True)

    # labeling the samples that are above LOD == 0.005 for Saccharin biotransforamtion
    df_rates = pd.read_csv("data/20260325_bt_rates_matrix_without_AS.csv")

    # Load mapping from accession (e.g., run ID) to sample name
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

    df['label'] = 0
    label_samples = df_rates.loc[df_rates[compound] > 0.005, 'sample_name'].to_list()
    df.loc[df['sample_name'].isin(label_samples), 'label'] = 1


    feature_name = 'EC'
    significant_features_df = mann_whitney_u_test_features_labels(df, output_filename = f"{date}_mann_whitney_{feature_name}_{compound}_lod", feature_name = feature_name)
    top_n = 20
    significant_abu_df = get_significant_features_and_abundance_df(df, feature_name, significant_features_df[feature_name], sample_site_dict,
                                              sample_location_dict)

    significant_features_df_abu_info = get_EC_significance_df_with_abundances_info(significant_abu_df, significant_features_df)
    significant_features_df_abu_info.to_csv(f"output/mann_whitney_{compound}_above_vs_below_lod/{date}_mann_whitney_{feature_name}_{compound}_lod_with_abundance_info.csv")


    significant_feature_target = significant_feature_df_for_list_of_features(feature_name, target_EC, significant_features_df)
    significant_feature_target.to_csv(f"output/mann_whitney_{compound}_above_vs_below_lod/{date}_mann_whitney_{feature_name}_{compound}_lod_with_abundance_info_target_EC.csv")
    #

    df.drop(columns='sample_name', inplace=True)
    top_n_2 =10
    #significant & target

    plot_boxplots(significant_feature_target, df, feature_name, top_n_2, plot_title=f'Targeted ECs within Significant EC ({compound} above LOD vs rest)' ,
             path_fig=f"output/mann_whitney_{compound}_above_vs_below_lod/{date}_{compound}_mann_whitney_boxplots_{top_n_2}_lod_vs_rest_target_{feature_name}.png")
    #
    # significant & target & enriched
    plot_boxplots(significant_feature_target.loc[significant_feature_target['delta_mean'] > 0], df, feature_name, significant_feature_target.loc[significant_feature_target['delta_mean'] > 0].shape[0] ,
                  plot_title=f'Targeted and Enriched ECs within Significant EC ("{compound} above lod vs rest")',
                  path_fig=f"output/mann_whitney_{compound}_above_vs_below_lod/{date}_{compound}_mann_whitney_boxplots_lod_vs_rest_target_{feature_name}_only_enriched.png")

    # significant & enriched
    plot_boxplots(significant_features_df_abu_info.loc[significant_features_df_abu_info['delta_mean'] > 0], df, feature_name,
                  top_n,
                  plot_title=f'Enriched ECs within Significant EC ("{compound} above lod vs rest")',
                  path_fig=f"output/mann_whitney_{compound}_above_vs_below_lod/{date}_{compound}_mann_whitney_boxplots_lod_vs_rest_only_enriched.png")


    print('done')