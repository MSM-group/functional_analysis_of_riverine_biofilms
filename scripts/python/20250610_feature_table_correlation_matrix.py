import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

date = 20260206
compound = 'Cyclamate'

feature_df_org = pd.read_csv(f"data/20260121_{compound}_river_ec_rates_df.csv")

#dropping samples without rate constants
feature_df = feature_df_org.drop(columns = ['accession', 'sample_name'])
if pd.isna(feature_df['rate_constants']).any():
    feature_df = feature_df.dropna(axis = 0, subset = 'rate_constants')

#In my case the categorical values are nominal, therefor I will need to use dummy encoding
categorical_features = ['location', 'site', 'field_campaign', 'read_count']
categorical_df = pd.concat([feature_df.pop(x) for x in categorical_features], axis=1)

#Explore the correlations between the features and the constant rates
correlation_matrix_spearman = feature_df.corr(method = 'spearman')

corr_rate_spearman = correlation_matrix_spearman['rate_constants']


if compound == 'Cyclamate':
    fontsize = 10


#change column names for plotting
correlation_matrix_spearman.rename(columns = {'silicic_acid' : 'Salicylic acid'})
for column in correlation_matrix_spearman.columns:
    if column in correlation_matrix_spearman.columns:
        if column != 'pH':
            format_col = column.replace("_", " ")
            format_col = format_col.capitalize()
            correlation_matrix_spearman.rename(columns = {column : format_col}, inplace = True)



# Plotting the correlation matrix-spearman
plt.figure(figsize=(20, 20))
plt.rcParams['font.family'] = 'Arial'
sns.clustermap(correlation_matrix_spearman, annot = True, fmt = ".2f", cmap='coolwarm', annot_kws={"fontsize":fontsize}, )
plt.title('Spearman', fontsize = 8)
#plt.tight_layout()
plt.savefig(f"output/{date}_{compound}_spearman_cluster_correaltion_matrix.svg")

correlation_matrix_spearman.to_csv(f"output/{date}_{compound}_spearman_cluster_correaltion_matrix.csv")

print('done')