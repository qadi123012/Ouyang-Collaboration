# IR Chemometrics for Mixed Salt Solutions

![Python](https://img.shields.io/badge/python-3.9%2B-blue) ![MATLAB](https://img.shields.io/badge/MATLAB-legacy%20scripts-orange) ![License](https://img.shields.io/badge/license-MIT-green) ![Status](https://img.shields.io/badge/status-research%20code-lightgrey)

**Quantifying the concentration of each salt in a mixed aqueous solution from its infrared spectrum, with partial least squares (PLS) regression, and checking honestly how well it works.**

The project covers four mixed-salt systems: potassium carbonate and bicarbonate, potassium oxalate and hydrogen oxalate, ammonium sulfate and bisulfate, and the three sodium phosphates. It contains the spectra, a reproducible Python pipeline, the original MATLAB scripts, and a full validation study with figures and tables.

---

## Contents

1. [Why this project](#why-this-project) 
2. [Read this first](#read-this-first)
3. [Results at a glance](#results-at-a-glance)
4. [Repository structure](#repository-structure)
5. [Quick start](#quick-start)
6. [Data](#data)
7. [Method in brief](#method-in-brief)
8. [Figures](#figures)
9. [Using the code on your own spectra](#using-the-code-on-your-own-spectra)
10. [MATLAB scripts and legacy material](#matlab-scripts-and-legacy-material)
11. [Limitations](#limitations)
12. [Documentation](#documentation) | [Citation](#citation) | [License](#license) | [Contact](#contact)

## Why this project

Many industrial and environmental solutions contain a salt together with its protonated partner, for example carbonate and bicarbonate. Their infrared bands overlap strongly, so a single-peak reading is rarely enough. A multivariate calibration can in principle separate them, but with only a dozen standards it is easy to fool oneself: a model can fit its own training data perfectly and still predict new samples badly.

This repository therefore asks three practical questions for each system:

* Can the concentrations of the two (or three) salts be predicted from a spectrum?
* How large is the error on samples the model has **never seen**?
* Does PLS actually beat a simple single-band calibration, or even just guessing the average?

## Read this first

Three results matter more than any single number below.

1. **Potassium oxalate works well.** Errors of 0.009 and 0.029 mol/L on a 0 to 0.5 mol/L range, on unseen samples.
2. **The ammonium sulfate validation set appears to be on a different concentration scale.** The model ranks the validation samples almost perfectly (r = 0.999) but predicts values about **2 times** the labels. A likely explanation is a labelling issue in the validation samples (hypothesis, see [findings](docs/findings_and_limitations.md)); the preparation records should be checked before these numbers are quoted.
3. **The carbonate and phosphate data are too small for reliable PLS.** Cross-validation cannot pick the model settings reliably with 11 to 17 standards, and a single-band calibration is often as good or better.

The statistics printed by the original MATLAB scripts are R2 and RMSE of the **fitted** data. They are not measures of predictive ability, and this project does not use them as such.

## Results at a glance

PLS regression with one model per analyte, pre-processing and model size chosen by leave-one-composition-out cross-validation on the calibration standards only, then applied once to the unseen validation samples. All errors are in mol/L.

| System | Analyte | Calibration range (mol/L) | Pipeline | LV | RMSECV | RMSEP | Mean predictor | Single band | r (val) | Slope (val) |
|---|---|---|---|---|---|---|---|---|---|---|
| K2CO3 + KHCO3 | K2CO3 | 0 to 1.5 | `A\|fingerprint\|none\|d0` | 6 | 0.338 | 0.498 | 0.384 | 0.324 | 0.46 | 0.68 |
| K2CO3 + KHCO3 | KHCO3 | 0 to 1.5 | `A\|fingerprint\|none\|d0` | 3 | 0.449 | 0.189 | 0.374 | 0.099 | 0.89 | 0.78 |
| K2C2O4 + KHC2O4 | K2C2O4 | 0 to 0.5 | `T\|full\|none\|d0` | 4 | 0.011 | 0.009 | 0.118 | 0.016 | 1.00 | 1.02 |
| K2C2O4 + KHC2O4 | KHC2O4 | 0 to 0.5 | `T\|full\|none\|d0` | 4 | 0.056 | 0.029 | 0.104 | 0.048 | 0.98 | 1.20 |
| (NH4)2SO4 + NH4HSO4 | (NH4)2SO4 | 0 to 2 | `A\|fingerprint\|snv\|d0` | 3 | 0.089 | 0.672 | 0.564 | 0.724 | 1.00 | 2.34 |
| (NH4)2SO4 + NH4HSO4 | NH4HSO4 | 0 to 2 | `A\|fingerprint\|snv\|d0` | 7 | 0.057 | 0.584 | 0.564 | 0.709 | 1.00 | 1.88 |
| Na3PO4 + Na2HPO4 + NaH2PO4 | Na3PO4 | 0 to 0.5 | `T\|full\|snv\|d1` | 2 | 0.086 | 0.069 | 0.112 | 0.057 | 0.87 | 0.42 |
| Na3PO4 + Na2HPO4 + NaH2PO4 | Na2HPO4 | 0 to 0.5 | `T\|full\|snv\|d1` | 1 | 0.210 | 0.177 | 0.137 | 0.125 | -0.40 | -0.07 |
| Na3PO4 + Na2HPO4 + NaH2PO4 | NaH2PO4 | 0 to 0.5 | `T\|full\|snv\|d1` | 4 | 0.083 | 0.077 | 0.177 | 0.123 | 0.77 | 0.48 |

*LV = latent variables. RMSECV = cross-validation error. RMSEP = error on unseen samples. "Mean predictor" is the RMSEP of simply guessing the calibration average, and a model must beat it. "Single band" is a one-wavenumber calibration. Pipeline labels read signal, region, scatter correction, derivative; `A` is absorbance, `T` the signal as recorded. r and slope compare predicted with labelled concentrations on the validation samples.*

| System | Verdict |
|---|---|
| K2C2O4 + KHC2O4 | Good. Beats the mean predictor in all 24 pre-processing pipelines, and the single band in the selected one |
| Na3PO4, NaH2PO4 | Partly. Better than the mean in all 24 pipelines. Na2HPO4 is not predicted in any |
| K2CO3 + KHCO3 | Unstable. The CV-selected model does not beat the mean for K2CO3 (0.498 vs 0.384), although several SNV pipelines do far better on validation; CV on 9 compositions could not pick them |
| (NH4)2SO4 + NH4HSO4 | Cross-validation R2 0.98 and 0.99 but validation fails; scale mismatch, see above |

Full discussion: [`docs/findings_and_limitations.md`](docs/findings_and_limitations.md).

## Repository structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── run_all.sh                      regenerate every table and figure (about 75 s)
├── data/
│   ├── raw/                        11 spectral workbooks, grouped by chemical system
│   └── processed/                  SNV-preprocessed workbooks written by the MATLAB step
├── irchem/                         Python package
│   ├── io.py                       dataset registry and loader
│   ├── preprocessing.py            absorbance, region, SNV, Savitzky-Golay, baseline
│   ├── models.py                   PLS1, single-band and mean baselines
│   └── evaluation.py               grouped cross-validation and metrics
├── scripts/
│   ├── 01_describe_datasets.py     data catalogue and raw-spectra figure
│   ├── 02_run_benchmark.py         24 pipelines x every dataset
│   ├── 03_make_figures.py          parity, CV, coefficient and diagnostic figures
│   └── 04_metric_audit_vs_matlab.py  what the MATLAB statistics measure
├── tests/                          21 unit tests
├── results/
│   ├── tables/                     all numbers as CSV
│   └── figures/                    F1 to F7
├── matlab/                         original MATLAB scripts, renamed and documented
├── legacy/                         earlier outputs and scripts whose inputs or code are missing
└── docs/                           methods, data catalogue, findings and limitations
```

## Quick start

```bash
git clone https://github.com/qadi123012/IR-Chemometrics-Salt-Solution-Quantification.git
cd IR-Chemometrics-Salt-Solution-Quantification
pip install -r requirements.txt

python -m pytest -q tests        # 21 tests, a few seconds
bash run_all.sh                  # tables and figures, about 75 seconds on one CPU core
```

Each step also runs on its own, for example `python scripts/02_run_benchmark.py --datasets oxalate`. Results land in `results/tables/` and `results/figures/`.

## Data

Eleven Excel workbooks, one per measurement campaign. Each stores the sample table (concentrations in mol/L) and one spectrum per sample, 550 to 4000 cm-1 at about 1.9 cm-1 spacing.

* `R` samples are calibration standards; `Y` samples are unseen validation samples.
* The recorded signal looks like percent transmittance (baseline near 100, bands as dips), so the code can convert it to absorbance.
* Some workbooks still carry the template headers `c_K2CO3` and `c_KHCO3` although they hold other salts. The code reads concentrations by position and states the real chemistry.
* Two spectra (R3 and R10) of the K2CO3 + KHCO3 set were replaced in the corrected file, so earlier versions are not comparable.

The complete list, with original file names and caveats, is in [`docs/data_catalogue.md`](docs/data_catalogue.md).

## Method in brief

* **Pre-processing grid (24 pipelines):** recorded signal or absorbance, three wavenumber regions, optional SNV, optional Savitzky-Golay first derivative. The CO2 window 2220 to 2430 cm-1 is always excluded.
* **Models:** PLS regression with one model per analyte; baselines: best single wavenumber (chosen inside each fold) and the mean predictor.
* **Validation:** leave-one-**composition**-out cross-validation on the R standards (identical centre-point standards are held out together, otherwise near duplicates leak into training), parsimonious choice of latent variables, a single pass over the Y samples. The pipeline is selected by cross-validation only.

Details: [`docs/methods.md`](docs/methods.md).

## Figures

**Raw spectra of the four systems** (`F1`)

![Raw spectra](results/figures/F1_spectra_overview.png)

**Predicted versus labelled concentration** (`F3`): calibration points are cross-validated predictions, triangles are unseen samples.

![Parity plots](results/figures/F3_parity_plots.png)

**The ammonium sulfate scale problem** (`F5`)

![Sulfate diagnostic](results/figures/F5_sulfate_scale_diagnostic.png)

**Validation error of all 24 pipelines relative to the mean predictor** (`F4`)

![Pipeline heatmap](results/figures/F4_pipeline_heatmap.png)

Also in `results/figures/`: cross-validation curves (`F2`), PLS regression coefficients (`F6`), and PCA of calibration versus validation spectra (`F7`).

## Using the code on your own spectra

```python
import numpy as np
from irchem import preprocessing as P, evaluation as E
from irchem.models import PLS1

# X: (n_spectra, n_points) transmittance in %, wn: wavenumbers (ascending), y: concentrations, groups: composition ids
pipe = P.Pipeline(signal="A", scatter="snv", deriv=0, region="fingerprint")
w, Z = pipe.apply(wn, X)

curve = E.pls_cv_curve(Z, y, groups, kmax=8)       # leave-one-composition-out RMSECV for 1..8 LV
k = E.choose_k(curve)                              # parsimonious choice
model = PLS1(k).fit(Z, y)
y_new = model.predict(pipe.apply(wn, X_new)[1])
```

To add a new dataset, place the workbook under `data/raw/` in the same layout and add one line to the registry in `irchem/io.py`.

## MATLAB scripts and legacy material

* [`matlab/`](matlab/README.md) holds the nine original scripts (preprocessing, PLS, classical least squares, univariate calibration), unchanged but renamed, with a table of what each does and a list of known issues.
* [`legacy/`](legacy/README.md) holds figures and code from an earlier analysis (peak fitting, ANN). Their inputs or code are missing, so they are not reproduced here and should be read as exploratory.

## Limitations

* One validation set per system, with 5 to 9 samples. R2 on so few points is noisy.
* Calibration sets of 11 to 17 standards for about 1700 spectral variables.
* No instrument metadata (technique, accessory, path length, temperature, pH) was available; the transmittance interpretation rests on the spectra themselves.
* Formal concentrations are used as labels. Acid-base equilibria between each salt and its protonated partner mean the spectra may depend on speciation, which a linear model cannot capture.
* The ANN, peak-fitting and any CNN work mentioned in earlier versions of this repository are not part of the current pipeline.

## Documentation

| File | Content |
|---|---|
| [`docs/methods.md`](docs/methods.md) | every step, metric and the validation protocol |
| [`docs/data_catalogue.md`](docs/data_catalogue.md) | all workbooks, original names, caveats |
| [`docs/findings_and_limitations.md`](docs/findings_and_limitations.md) | results per system, diagnostics, next steps |
| [`matlab/README.md`](matlab/README.md) | MATLAB scripts and their known issues |
| [`legacy/README.md`](legacy/README.md) | earlier outputs and what is missing |

## Citation

If this repository helps your work, please cite it (see [`CITATION.cff`](CITATION.cff)).

## License

Code: MIT, see [`LICENSE`](LICENSE). The spectra in `data/` are provided for reproducibility of the analyses; please acknowledge the source when reusing them.

## Contact

Abdul Qadir. Please open a GitHub Issue for questions, bug reports or suggestions.

## Acknowledgements

The solutions were prepared and measured by the collaborating group whose name appears in the data file names (OUYANG). This repository was previously named *Ouyang-Collaboration*.
