import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# --- PARAMETERS ---
# Choose which absorbance column to use for calibration (e.g., 'r1-1')
calib_col = 'r1-1'  # Change as needed
pred_col = 'Y1-1'   # Change as needed for prediction  

# 1. Load training data
train = pd.read_excel('0611-train-new.xlsx')   
# 2. Find peak wavenumber for the chosen calibration column
peak_idx = train[calib_col].idxmax()
peak_wavenumber = train.loc[peak_idx, 'wave number']
print(f"Peak wavenumber for {calib_col}: {peak_wavenumber}")

# 3. Extract absorbance at the peak for all training samples
abs_cols = [col for col in train.columns if col.startswith('r') or col.startswith('RR')]
peak_absorbances = []
concentrations = []  
sample_names = []
for col in abs_cols: 
    idx = train[col].idxmax() 
    peak_abs = train.loc[idx, col]
    peak_wv = train.loc[idx, 'wave number']
    # Use K2CO3 as example
    conc = train.loc[abs_cols.index(col), 'c_K2CO3'] 
    peak_absorbances.append(peak_abs)
    concentrations.append(conc)
    sample_names.append(col)

# 4. Fit linear regression model
X_train = np.array(peak_absorbances).reshape(-1, 1)
y_train = np.array(concentrations)
model = LinearRegression()
model.fit(X_train, y_train)

# 5. Plot calibration curve
plt.figure(figsize=(8, 6))
plt.scatter(X_train, y_train, color='blue', label='Training Samples')
plt.plot(X_train, model.predict(X_train), color='red', label='Linear Fit')
for i, name in enumerate(sample_names):
    plt.text(X_train[i], y_train[i], name, fontsize=8, color='blue', ha='right')
plt.xlabel('Peak Absorbance')
plt.ylabel('K2CO3 Concentration')
plt.title(f'Calibration Curve at Peak (wavenumber {peak_wavenumber:.1f})')
plt.legend()

# 6. Load prediction data
pred = pd.read_excel('0611-verify & predict-new.xlsx')

# 7. Extract absorbance at the same wavenumber for each prediction sample
pred_cols = [col for col in pred.columns if col != 'wave number']
pred_absorbances = []
pred_sample_names = []
for col in pred_cols:
    # Find the row closest to the calibration peak wavenumber
    idx = (pred['wave number'] - peak_wavenumber).abs().idxmin()          
    pred_abs = pred.loc[idx, col]
    pred_absorbances.append(pred_abs)
    pred_sample_names.append(col)                

# 8. Predict concentrations  
X_pred = np.array(pred_absorbances).reshape(-1, 1)
pred_conc = model.predict(X_pred)

# 9. Plot predicted points        
plt.scatter(X_pred, pred_conc, color='green', marker='x', label='Predicted Samples')
for i, name in enumerate(pred_sample_names):
    plt.text(X_pred[i], pred_conc[i], name, fontsize=8, color='green', ha='left')
plt.legend()
plt.tight_layout()
plt.savefig('concentration_prediction_plot.png', dpi=300)
plt.close()

# 10. Save predicted concentrations
pred_results = pd.DataFrame({'Sample': pred_sample_names, 'Predicted_Concentration': pred_conc})
pred_results.to_excel('predicted_concentrations.xlsx', index=False)

print('Prediction complete. Results saved to predicted_concentrations.xlsx and concentration_prediction_plot.png.')
