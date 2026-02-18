from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


compound = 'Saccharin' #or Cyclamate
date = 20260218
plot_top_n = 10
drop_genomic_features = False
drop_river_parameters = True


feature_df = pd.read_csv(f"data/20260121_{compound}_river_ec_rates_df.csv")
fig_path = f"output/{date}_{compound}_feature_importance_entire_set_top_{plot_top_n}.svg"


#dropping categotical features
categorical_features = ['location', 'site', 'field_campaign', 'accession', 'sample_name','read_count']
categorical_df = pd.concat([feature_df.pop(x) for x in categorical_features], axis=1)

if drop_river_parameters:
    non_gene_cols = [col for col in feature_df.columns if '.' not in col]
    non_gene_cols.remove('rate_constants')
    feature_df.drop(columns=non_gene_cols, inplace=True)
    fig_path = f"output/{date}_{compound}_feature_importance_entire_set_top_{plot_top_n}_only_genomic_features.svg"

if drop_genomic_features:
    ec_cols = [col for col in feature_df.columns if '.' in col]
    feature_df.drop(columns=ec_cols, inplace=True)
    fig_path = f"output/{date}_{compound}_feature_importance_entire_set_top_{plot_top_n}_only_river_parameters.svg"

#dropping rows with missing reaction rate constants
if pd.isna(feature_df['rate_constants']).any():
    feature_df = feature_df.dropna(axis = 0, subset = 'rate_constants')

#Formatting for nicer plotting

column_list = ['conductivity', 'alkalinity', 'hardness', 'sodium', 'magnesium',
       'calcium', 'potassium', 'chloride', 'nitrate', 'sulfate', 'ammonium',
       'phosphate', 'total_phosphorus', 'total_inorganic_carbon']


feature_df.rename(columns = {'silicic_acid' : 'Salicylic acid'})
for column in feature_df.columns:
    if column in column_list:
        format_col = column.replace("_", " ")
        format_col = format_col.capitalize()
        feature_df.rename(columns = {column : format_col}, inplace = True)

X = feature_df.drop(columns = 'rate_constants')
y = feature_df['rate_constants']

#Train Random Forest
rf = RandomForestRegressor(n_estimators=100, random_state=42)
rf.fit(X, y)

# Predict and evaluate
y_pred = rf.predict(X)
rmse = mean_squared_error(y, y_pred)
r2 = r2_score(y, y_pred)

print(f"RMSE (Trained on entire set): {rmse:.3f}")
print(f"R² (Trained on entire set): {r2:.3f}")

# Feature importance
importances = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)

# if plot_top_n:
#     importances = importances[:plot_top_n]
#     fig_path = f"output/{date}_{compound}_feature_importance_entire_set_top_{plot_top_n}.svg"
#
# else:
#     fig_path = f"output/{date}_{compound}_feature_importance_entire_set.svg"

if plot_top_n:
    plt.figure(figsize=(13, 8))
    plt.rcParams['font.family'] = 'Arial'


sns.barplot(x = importances, y = importances.index, color = "#7a9df8ff")

# sns.barplot(x = importances, y = importances.index, palette="coolwarm")
plt.title(f"{compound} Feature Importances from Random Forest, Entire Set", fontsize = 18)
plt.xlabel("Importance", fontsize = 22)
plt.ylabel("Feature", fontsize = 22)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)

# Add text box with RMSE and R²
metrics_text = f"RMSE = {rmse:.3f}\nR² = {r2:.3f}"
plt.text(0.95, 0.05, metrics_text, transform=plt.gca().transAxes,
         fontsize=24, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_path)

