from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score

from sklearn.metrics import accuracy_score, f1_score, classification_report
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.inspection import permutation_importance

compound = 'Saccharin'
date = 20260428
# with_dropping_low_correlating_features = False
plot_top_n = 10

feature_df = pd.read_csv(f"data/20260428_{compound}_river_ec_rates_df.csv")

#dropping categotical features
categorical_features = ['location', 'site', 'field_campaign', 'accession', 'sample_name','read_count']
categorical_df = pd.concat([feature_df.pop(x) for x in categorical_features], axis=1)

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

threshold = 0.005
X = feature_df.drop(columns = 'rate_constants')
y = (feature_df['rate_constants'] > threshold).astype(int)

#Train Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
rf.fit(X, y)

# Predict and evaluate
y_pred = rf.predict(X)
accuracy = accuracy_score(y, y_pred)
f1 = f1_score(y, y_pred)
print(f"Accuracy (Entire set): {accuracy:.3f}")


# Feature importance
importances = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)

if plot_top_n:
    importances = importances[:plot_top_n]
    fig_path = f"output/{date}_{compound}_feature_importance_entire_set_top_{plot_top_n}_wo_categorical_and_reads.svg"

else:
    fig_path = f"output/{date}_{compound}_feature_importance_entire_set_wo_categorical_and_reads.svg"

if plot_top_n:
    plt.figure(figsize=(13, 8))
    plt.rcParams['font.family'] = 'Arial'


sns.barplot(x = importances, y = importances.index, color = "#7a9df8ff")

plt.title(f"{compound} Feature Importances from Random Forest, Entire Set", fontsize = 18)
plt.xlabel("Importance", fontsize = 22)
plt.ylabel("Feature", fontsize = 22)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)

# Add text box with RMSE and R²
metrics_text = f"accuracy = {accuracy:.3f}\nF1 = {f1:.3f}"
plt.text(0.95, 0.05, metrics_text, transform=plt.gca().transAxes,
         fontsize=24, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_path)


#######################################################################
#Permutation feature importance
#######################################################################

# Stratified split is crucial to keep class proportions equal
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

clf_perm = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
clf_perm.fit(X_train, y_train)

# Calculate metrics for the test set
y_test_pred = clf_perm.predict(X_test)
test_f1 = f1_score(y_test, y_test_pred)
test_acc = accuracy_score(y_test, y_test_pred)

print(f"Test Accuracy: {test_acc:.3f}")
print(f"Test F1-Score: {test_f1:.3f}")

# We now use 'f1' as the scoring metric for permutation
perm_importance = permutation_importance(
    clf_perm, X_test, y_test, n_repeats=30, random_state=42, scoring='f1'
)

importances_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance_Mean': perm_importance.importances_mean,
    'Importance_Std': perm_importance.importances_std
}).sort_values(by='Importance_Mean', ascending=False)

if plot_top_n: importances_df = importances_df.head(plot_top_n)
importances_df['Feature'] = importances_df['Feature'].str.replace("_", " ")

plt.figure(figsize=(10, 6))
plt.barh(y=importances_df['Feature'], width=importances_df['Importance_Mean'],
         xerr=importances_df['Importance_Std'], color='skyblue', capsize=4)
plt.gca().invert_yaxis()
plt.xlabel("Mean Decrease in F1-Score", fontsize=14)
plt.title(f"{compound} Permutation Importance (Scored by F1)", fontsize=16)

# Add metrics text box
plt.text(0.95, 0.05, f"Test F1 = {test_f1:.3f}", transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(f"output/{date}_{compound}_permutation_f1.svg")

print('Done. See classification report below:')
print(classification_report(y_test, y_test_pred))

print('done')