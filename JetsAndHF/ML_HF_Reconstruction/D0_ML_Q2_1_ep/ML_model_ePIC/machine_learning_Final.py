#!/bin/bash
# Code for ML for ePIC using hipe4ml
#https://doi.org/10.5281/zenodo.5070131
#Install: pip install hipe4ml or brew install libomp
#script to run the machine learning differential in pT and y
# Shyam Kumar; INFN Bari, Italy;
# shyam.kumar@ba.infn.it; shyam055119@gmail.com
# Supported by The FAIR Spoke 6 Project, funded by the NextGenerationEU program in Italy
import numpy as np
import os
import shap
import argparse
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import plot_utils, analysis_utils
from ROOT import TFile, TH1F
from sklearn.metrics import log_loss, roc_curve, roc_auc_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

promptDmeson = 10000
def print_optuna_results(study):
    print(f"{'| Trial':<10}| {'Target':<9}| {'learning_rate':<14}| {'max_depth':<10}| {'n_estimators':<14}|")
    print("-" * 70)
    
    for trial in study.trials:
        target = round(trial.value, 4)
        lr = round(trial.params.get('learning_rate', 0), 5)
        depth = int(trial.params.get('max_depth', 0))
        estimators = int(trial.params.get('n_estimators', 0))

        print(f"| {trial.number:<9}| {target:<9}| {lr:<14}| {depth:<10}| {estimators:<14}|")

# Add an argument parser
parser = argparse.ArgumentParser(description="ML Script with pT and y ranges")
parser.add_argument('--ptmin',type=float,default=1.0, help='Pass the pTmin')
parser.add_argument('--ptmax',type=float,default=10.0, help='Pass the pTmax')
parser.add_argument('--ymin',type=float,default=-1.0, help='Pass the ymin')
parser.add_argument('--ymax',type=float,default=1.0, help='Pass the ymax')
args = parser.parse_args()
print(f"pT Range: {args.ptmin}-{args.ptmax} GeV/c, y Range: {args.ymin}-{args.ymax}")

# Create directory if it doesn't exist
output_dir = f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/massD0_ML_ePIC'
os.makedirs(output_dir, exist_ok=True)

sig = TreeHandler('Data_Preparation/SignalD0.root','treeMLSig')
bkg = TreeHandler('Data_Preparation/BkgD0.root','treeMLBkg')
print(f"Available Signal candidates: {sig.get_n_cand()} \t Available Background candidates: {bkg.get_n_cand()}")
selsig_cuts = sig.get_subset(f'(pt_D0 > {args.ptmin} and pt_D0 < {args.ptmax}) and (y_D0 > {args.ymin} and y_D0 < {args.ymax}) and (mass_D0 > 1.7 and mass_D0 < 2.1)') # ,size=promptDmeson
selbkg_cuts = bkg.get_subset(f'(pt_D0 > {args.ptmin} and pt_D0 < {args.ptmax}) and (y_D0 > {args.ymin} and y_D0 < {args.ymax}) and (1.6< mass_D0 < 1.7 or 2.1 < mass_D0 < 2.5)') # , size=selsig.get_n_cand()

# Get sizes
n_sig = selsig_cuts.get_n_cand()
n_bkg = selbkg_cuts.get_n_cand()
min_size = min(n_sig, n_bkg)

# Resize to match minimum
selsig = selsig_cuts.get_subset(size=min_size)
selbkg = selbkg_cuts.get_subset(size=min_size)

msg = f"Signal candidates for ML: {selsig.get_n_cand()} \t Background candidates for ML: {selbkg.get_n_cand()}"
print(msg)
# Construct the filename dynamically
filename = f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/sig_bkg_candidates.txt'
# Write to the dynamic file
with open(filename, "w") as f:
    f.write(msg + "\n")

print(f"Signal candidates for ML: {selsig.get_n_cand()} \t Background candidates for ML: {selbkg.get_n_cand()}")

#quit() # break from the code for debug
train_test_data = train_test_generator([selsig, selbkg], [1,0], test_size=0.2, random_state=42)
vars_to_draw = selsig.get_var_names()
#print("Variables",vars_to_draw)
#print(f"Training Data \n: {pd.DataFrame(train_test_data[0]),pd.DataFrame(train_test_data[1])}, \n Test Data \n: {pd.DataFrame(train_test_data[2]),pd.DataFrame(train_test_data[3])}")

leg_labels = ['Background', 'Signal']
plot_utils.plot_distr([selbkg, selsig], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/sig_bkg_distributions_plot.png') 
plt.close() 

features_for_train = [v for v in vars_to_draw if v not in ['mass_D0','pt_D0','d0xy_k','y_D0','dca_12','sum_d0xy']]
print("Variables used in Machine learning",features_for_train)

plot_utils.plot_corr([selbkg, selbkg], features_for_train, leg_labels[:1])
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bkg_correlation_plot.png', bbox_inches='tight') 
plt.close() 

plot_utils.plot_corr([selsig, selsig], features_for_train, leg_labels[:2])
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/signal_correlation_plot.png', bbox_inches='tight')  
plt.close()  

model_clf = xgb.XGBClassifier(verbosity=1)
model_hdl = ModelHandler(model_clf, features_for_train)
print("Features used for training (debug):", model_hdl.get_training_columns())

hyper_pars_ranges = {'n_estimators': (100, 500), 'max_depth': (1, 3), 'learning_rate': (0.01, 0.1)}
study = model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc', timeout=120, n_jobs=-1, n_trials=100, direction='maximize')

# Print all the information of the optimisation
print_optuna_results(study)

# Extract results from Optuna study
trials_data = []
for trial in study.trials:
    trials_data.append({
        "trial": trial.number,
        "n_estimators": trial.params["n_estimators"],
        "max_depth": trial.params["max_depth"],
        "learning_rate": trial.params["learning_rate"],
        "value": trial.value
    })

df = pd.DataFrame(trials_data)

# Plot ROC Value vs. n_estimators
plt.figure(figsize=(12, 10))
plt.scatter(df["n_estimators"], df["value"], c="blue", label="ROC AUC")
plt.xlabel("n_estimators")
plt.ylabel("Value (ROC AUC)")
plt.title("Value vs. n_estimators")
plt.legend()
plt.grid(True)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/trainvalue_vs_estimators.png')
plt.close()

# Plot ROC Value vs. Trial
plt.figure(figsize=(12, 10))
plt.scatter(df["trial"], df["value"], c="red", label="ROC AUC")
plt.xlabel("Trial")
plt.ylabel("Value (ROC AUC)")
plt.title("Value vs. Trial")
plt.legend()
plt.grid(True)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/trainvalue_vs_trail.png')
plt.close()

# Plot ROC Value vs. max_depth
plt.figure(figsize=(12, 10))
plt.scatter(df["max_depth"], df["value"], c="green", label="ROC AUC")
plt.xlabel("max_depth")
plt.ylabel("Value (ROC AUC)")
plt.title("Value vs. max_depth")
plt.legend()
plt.grid(True)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/trainvalue_vs_maxdepth.png')
plt.close()

# Plot ROC Value vs. learning_rate
plt.figure(figsize=(12, 10))
plt.scatter(df["learning_rate"], df["value"], c="purple", label="ROC AUC")
plt.xlabel("learning_rate")
plt.ylabel("Value (ROC AUC)")
plt.title("Value vs. learning_rate")
plt.legend()
plt.grid(True)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/trainvalue_vs_learningrate.png')
plt.close()
                               
model_hdl.train_test_model(train_test_data)
y_pred_train = model_hdl.predict(train_test_data[0], False)
y_pred_test = model_hdl.predict(train_test_data[2], False)
print("Train data pediction: ",y_pred_train)
print("Test data pediction: ",y_pred_test)

plt.rcParams["figure.figsize"] = (10, 7)
#plt.savefig('ROC_AUC.png')  # Save the second plot as a PNG
plt.close()  # Close the second figure

ml_out_fig = plot_utils.plot_output_train_test(model_hdl, train_test_data, 100, 
                                               False, leg_labels, True, density=True) 
roc_train_test_fig = plot_utils.plot_roc_train_test(train_test_data[3], y_pred_test,
                                                    train_test_data[1], y_pred_train, None, leg_labels)   


# Save the figures to files
ml_out_fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/ML_Output_Optuna_train_test.png')  # Saves ml_out_fig as a PNG file
roc_train_test_fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/roc_train_test.png')  # Saves roc_train_test_fig as a PNG file

plot_utils.plot_learning_curves(model_hdl,train_test_data,100)
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/learning_curves.png')
plt.close()  

# Create ROC using libraries
fpr, tpr, thresholds = roc_curve(train_test_data[1], y_pred_train) #ytrain and y_pred_train 

print(thresholds)

# Compute AUC
auc = roc_auc_score(train_test_data[1], y_pred_train)
print(f"ROC AUC: {auc:.2f}")

# Plot ROC Curve
plt.plot(fpr, tpr, marker='o', label=f"AUC = {auc:.2f}")
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')  # Diagonal line (random)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()
plt.grid()
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/ROC_AUC_library.png')  # Save the second plot as a PNG
plt.close()  # Close the second figure


# Ensure train_test_data[0] and train_test_data[2] are DataFrames and match model expectations
X_train, X_test = train_test_data[0], train_test_data[2]
y_train, y_test = train_test_data[1], train_test_data[3]


# Calculate confusion matrix
cm = confusion_matrix(y_test, (y_pred_test > 0.5).astype(int))

# Plot confusion matrix
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=leg_labels)
disp.plot(cmap='Blues')
plt.rcParams["figure.figsize"] = (10, 7)
plt.title('Confusion Matrix')
plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/confusion_matrix.png')
plt.close()

feature_names = X_test.columns.tolist()  # Extracts column names from DataFrame

# Compute log loss manually
train_logloss = log_loss(y_train, y_pred_train)
test_logloss = log_loss(y_test, y_pred_test)
print(f"Computed Training Log Loss: {train_logloss:.5f}")
print(f"Computed Validation Log Loss: {test_logloss:.5f}")


model = model_hdl.get_original_model()
explainer = shap.TreeExplainer(model)

# Shap is evaluated for selected feature for ML
X_train_sel = X_train[model.feature_names_in_]
X_test_sel = X_test[model.feature_names_in_]
y_train_sel = y_train[:len(X_train_sel)]  
y_test_sel = y_test[:len(X_test_sel)] 


print("X_train shape:", X_train_sel.shape)
print("X_test shape:", X_test_sel.shape)
print("Model input shape: \n", model.feature_importances_.shape)

shap_values = explainer.shap_values(X_test_sel)
shap_values_array = np.array(shap_values)  # Convert to NumPy array
print("Array Size, nfeatures, nclasses", shap_values_array.shape)  # Should be (6, len(X_test), 2)

base_prob= np.mean(model.predict(X_test_sel))  # Mean probability of class 1
print(f"Model Prediction from X_test (Mean): {base_prob} \t Base Value from SHAP: {explainer.expected_value}")

base_prob= np.mean(model.predict_proba(X_test_sel)[:, 1])  # Mean probability of class 1
print(f"Model Prediction from X_test (Mean probability): {base_prob} \t Base probability Value from SHAP = {1 / (1 + np.exp(-explainer.expected_value))}")

# Print SHAP values for the first 3 test samples
print("SHAP \n",shap_values, "\n Length", len(shap_values))
print("Test",X_test_sel, "\n Length", len(X_test_sel))
# Print SHAP values and corresponding predictions
for i in range(10):  # Analyze the first 3 test samples
    print(f"\n--- Test Sample {i} ---")
    print("X_test: \n", X_test_sel.iloc[i])
    print("y_test:", y_test_sel[i])
    print("SHAP Values:", shap_values[i])
    print("SHAP Sum:", shap_values[i].sum())
    print("Base Value:", explainer.expected_value)
    print("Actual Prediction:", model.predict(X_test_sel.iloc[i].values.reshape(1, -1))[0]) 
    print("SHAP Predicted Value:", explainer.expected_value + shap_values[i].sum())
    print("Predicted Probabilties:", model.predict_proba(X_test_sel.iloc[i].values.reshape(1, -1))[0]) 
     # Compute P(class 1) from SHAP predicted value using the sigmoid function
    shap_pred_prob = 1 / (1 + np.exp(-(explainer.expected_value + shap_values[i].sum())))
    print(f"SHAP Computed Probability: 1./[1+exp(-SHAP Predicted Value)] Background = {1.0-shap_pred_prob}, Signal = {shap_pred_prob}")


# Compute SHAP values for a single sample (X_test_sel must have exactly 1 row)
shap_values_sel = explainer.shap_values(X_test_sel.iloc[[0]])  # Note the double brackets to keep DataFrame format

# Use the result directly
shap_vals_sample = shap_values_sel[0]  # First and only row
#features_sample = X_test_sel.iloc[0].round(4) limit to 4 features
features_sample = X_test_sel.iloc[0].apply(lambda x: round(x, 2) if isinstance(x, (float, np.floating)) else x) # round decimals


# Confirm shape match
print("SHAP values length:", len(shap_vals_sample))
print("Feature length:", len(features_sample))

# Generate the force plot
fig = shap.plots.force(
    explainer.expected_value,   
    shap_vals_sample,
    features_sample,
    matplotlib=True,
    show=False
)

fig.set_size_inches(14, 5)  # wider and taller
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.4)  # more space at bottom for labels

for text in fig.findobj(match=plt.Text):
    text.set_fontsize(10)  # or try 12 for better readability

# Save the plot
fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/feature_shap_classifier_D0.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig)


shap_val_array = shap_values_sel[0]  # For the first sample

formatted_labels = pd.Series(
    [f"{v:.2f}" for v in shap_val_array],  # format to 2 decimal places
    index=X_test_sel.columns               # use feature names
)

# Plot with string-formatted SHAP values
fig = shap.plots.force(
    explainer.expected_value,
    shap_val_array,
    formatted_labels,         # this will now show 'costheta = 1.6500', etc.
    matplotlib=True,
    show=False
)

# Improve layout and readability
fig.set_size_inches(14, 5)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.4)

# Save the plot
fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/shap_classifier_D0.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig)

# Calculate the BDT efficiency as a function of the BDT score
EFFICIENCY, THRESHOLD = analysis_utils.bdt_efficiency_array(y_test, y_pred_test, 10000)
# Assuming EFFICIENCY and THRESHOLD are 1D NumPy arrays or lists
bdt_efficiency_signal = np.column_stack((THRESHOLD, EFFICIENCY))
# Save Threshold and Efficiency
np.savetxt(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_efficiency_signal_vs_threshold.txt', bdt_efficiency_signal, header="Threshold Efficiency", fmt='%.6f')

PRECISION_RECALL_PLOT = plot_utils.plot_precision_recall(y_test, y_pred_test)
PRECISION_RECALL_PLOT.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/precision.png')  # Saves ml_out_fig as a PNG file
BDT_EFFICIENCY_PLOT = plot_utils.plot_bdt_eff(THRESHOLD, EFFICIENCY)
BDT_EFFICIENCY_PLOT.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_efficiency.png')  # Saves ml_out_fig as a PNG file

# Calculate the BDT efficiency for BACKGROUND as a function of the BDT score
EFFICIENCY_BG, THRESHOLD_BG = analysis_utils.bdt_efficiency_array(1 - y_test, y_pred_test, 10000)

bdt_efficiency_bkg = np.column_stack((THRESHOLD_BG, EFFICIENCY_BG))
# Save Threshold and Efficiency
np.savetxt(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_efficiency_bkg_vs_threshold.txt', bdt_efficiency_bkg, header="Threshold Efficiency", fmt='%.6f')


# Plot BACKGROUND EFFICIENCY
BDT_EFFICIENCY_PLOT_BG = plot_utils.plot_bdt_eff(THRESHOLD_BG, EFFICIENCY_BG)
BDT_EFFICIENCY_PLOT_BG.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_efficiency_background.png')

# Calculate BACKGROUND REJECTION
REJECTION_BG = 1 - EFFICIENCY_BG

# Plot BACKGROUND REJECTION
BDT_REJECTION_PLOT_BG = plot_utils.plot_bdt_eff(THRESHOLD_BG, REJECTION_BG)
BDT_REJECTION_PLOT_BG.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_rejection_background.png')

#Superimposed
plt.figure(figsize=(12, 10))
plt.plot(THRESHOLD, EFFICIENCY, label='Signal Efficiency', color='red', linewidth=2)
plt.plot(THRESHOLD_BG, EFFICIENCY_BG, label='Background Efficiency', color='blue', linestyle='--', linewidth=2)
plt.plot(THRESHOLD_BG, REJECTION_BG, label='Background Rejection', color='green', linestyle='-.', linewidth=2)

plt.xlabel('BDT Score Threshold')
plt.ylabel('Efficiency / Rejection')
plt.title('BDT Efficiency and Rejection')
plt.legend(loc='best')
plt.grid(True)

plt.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/bdt_eff_rej_superimposed.png')
plt.close()


# Plot feature importance correctly using `plot_feature_imp` with the label vector passed as y_truth
train_feature_plots = plot_utils.plot_feature_imp(X_train, y_train, model_hdl,leg_labels)
test_feature_plots = plot_utils.plot_feature_imp(X_test, y_test, model_hdl,leg_labels)

# Save each plot individually (assuming the function returns a list of figures)
if isinstance(train_feature_plots, list):
    for i, fig in enumerate(train_feature_plots):
        fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/features_train_{i}.png')
else:
    train_feature_plots.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/features_train.png')

if isinstance(test_feature_plots, list):
    for i, fig in enumerate(test_feature_plots):
        fig.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/features_test_{i}.png')
else:
    test_feature_plots.savefig(f'ML_Output_Optuna_{args.ymin}_{args.ymax}_{args.ptmin}_{args.ptmax}/features_test.png')    

