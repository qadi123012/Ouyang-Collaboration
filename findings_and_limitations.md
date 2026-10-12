# Findings and limitations

All numbers come from `results/tables/` and can be regenerated with `bash run_all.sh`. "Validation" always means the unseen Y samples. RMSE values are in mol/L.

## Summary

| System | Calibration standards | Result on unseen samples |
|---|---|---|
| K2C2O4 + KHC2O4 | 11 | **Works.** RMSEP 0.009 and 0.029 against 0.118 and 0.104 for the mean predictor; R2 0.994 and 0.913 |
| Na3PO4 + Na2HPO4 + NaH2PO4 | 17 | **Partly.** Na3PO4 and NaH2PO4 beat the mean in all 24 pipelines (RMSEP 0.069 and 0.077); Na2HPO4 does not in any |
| K2CO3 + KHCO3 | 11 | **Unstable.** 16 of 24 pipelines beat the mean for K2CO3 and 23 of 24 for KHCO3, but the CV-selected model does not (K2CO3 RMSEP 0.498 vs 0.384) |
| (NH4)2SO4 + NH4HSO4 | 11 | **Fails in validation, with a clear signature.** CV R2 0.98 and 0.99, but 0 of 24 pipelines beat the mean; predictions are about 2x the labels |

## 1. Oxalate is a good calibration

Using the signal as recorded, full range, no derivative, PLS with 4 latent variables gives RMSECV 0.011 / 0.056 and RMSEP 0.009 / 0.029 mol/L over a 0 to 0.5 mol/L range. A single-band calibration gives 0.016 / 0.048, so PLS is clearly better for K2C2O4 and slightly better for KHC2O4. Every one of the 24 pipelines beats the mean predictor. The validation slope is 1.02 and 1.20 with r of 0.998 and 0.983.

## 2. Ammonium sulfate: right ranking, wrong scale

Cross-validation on the 11 calibration standards is excellent (RMSECV 0.089 and 0.057, R2 0.98 and 0.99), yet validation fails for every pipeline. The diagnostic is clear (`F5_sulfate_scale_diagnostic.png`):

* correlation between predicted and labelled concentration on Y is 0.999 and 0.996,
* the slope is 2.34 and 1.88, so predictions are roughly twice the labels,
* at the 50/50 mixture the model predicts about 1.05 mol/L of each salt, the label says 0.50,
* the validation spectra lie largely inside the calibration space (`F7`; only the two extreme mixtures, Y1 and Y9, sit slightly outside in PC2 and in residual), so this is not a case of wild extrapolation.

**Hypothesis (not established by the files):** the Y mixtures were made from stock solutions of about 2 mol/L and labelled with mixing fractions that sum to 1.0. Multiplying the Y labels by 2 reduces RMSEP from 0.672 and 0.584 to 0.119 and 0.066 mol/L (`results/tables/sulfate_scale_hypothesis.csv`). This is a diagnostic, not a corrected result. The preparation records for Y1 to Y9 should be checked before any of these numbers are quoted.

The MATLAB-style fit statistics for this system (R2 0.84 and 0.81) are lower than the honest CV values because the MATLAB scripts train on R and Y together, so the inconsistent Y labels contaminate the fit.

## 3. K2CO3 + KHCO3: too little data for PLS

With 11 standards and 9 distinct compositions, the number of latent variables and the pre-processing cannot be chosen reliably by cross-validation. The CV-selected pipeline (absorbance, fingerprint region, 6 and 3 latent variables) has RMSEC 0.020 for K2CO3, a near perfect fit of the training data, but RMSECV 0.338 and RMSEP 0.498. A univariate calibration reaches RMSEP 0.324 and 0.099 on the same split, better than PLS. The validation correlation is 0.46 for K2CO3 and 0.89 for KHCO3.

The picture is not hopeless, only unstable. Looking at `F4_pipeline_heatmap.png`, several pipelines with SNV reach roughly 0.2 to 0.3 of the mean-predictor error on the validation samples, but leave-one-composition-out CV on 9 distinct compositions did not rank them first. Choosing one of them now would be selection on the validation set, so it is mentioned here only as a pointer: with more standards, CV might identify such a pipeline in advance.

The two earlier workbooks cannot be used to improve on this: they differ from the corrected one in the spectra of R3 and R10 (see `data_catalogue.md`).

Carbonate and bicarbonate are linked by an acid-base equilibrium, so the spectrum depends on speciation (pH, ionic strength) and not only on the formal concentrations used as labels. That would reduce additivity. The data cannot prove it, but the same caution applies to the other three systems.

## 4. Phosphate: two of three species

Na3PO4 and NaH2PO4 are predicted better than the mean in every pipeline. Na2HPO4 is not predicted at all (validation correlation -0.40, CV R2 -0.20). The three species are linked by proton-transfer equilibria, which makes the formal concentrations of the three salts only partly independent in the spectrum.

## 5. The statistics printed by the MATLAB scripts are fit statistics

`pls_regression_basic.m` and its siblings print R2 and RMSE of the fitted data after fitting on all samples. They say how well the model describes the data it has seen, not how well it predicts new samples (`results/tables/metric_audit_matlab_style_vs_honest.csv`). For oxalate the two views agree (R2 0.997 and 0.976 versus 0.994 and 0.913 on unseen samples); for sulfate and carbonate they do not.

## 6. Pre-processing matters, and no single choice wins

Winning pipelines differ by system (oxalate: recorded signal, full range; sulfate: absorbance, fingerprint, SNV; phosphate: SNV and first derivative). `F4_pipeline_heatmap.png` shows validation error relative to the mean predictor for all 24 pipelines.

## What the results do not show

* They do not show that PLS beats a single-band calibration in general; it does for oxalate and not for K2CO3 + KHCO3.
* 5 to 9 validation samples per system give noisy R2 values, and one validation set per system is all that exists.
* No instrument information is available (technique, accessory, path length, temperature, pH), so the transmittance interpretation rests on the shape and range of the spectra.
* The ANN, peak-fitting and CNN results from the original repository are not reproduced: their code is not in the repository. See `legacy/README.md`.

## Suggested next steps

1. Verify the preparation records of the sulfate Y samples (stock concentration and mixing volumes).
2. Add calibration standards, especially off-axis mixtures, and measure replicates of every standard, not only the centre point.
3. Record pH for the carbonate, phosphate and sulfate systems and consider it as a predictor.
4. Hold out the Y samples in the original design (as here) and keep them out of every model-selection step.
5. Report RMSECV and RMSEP, not R2 of the fit.
