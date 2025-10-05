# BaySIC-MATLAB

This repository provides MATLAB implementations of **BaySIC**: a Bayesian calibration for the Arctic sea-ice biomarker **IP<sub>25</sub>** and associated open-water phytoplankton-derived sterols (brassicasterol and dinosterol).  The models can be used to predict (forward model) the index **ln(PIP<sub>25</sub>)** from sea-ice concentration (**SIC**), or to invert measured biomarker concentrations into posterior wintertime (March-April-May) SIC estimates.

This code is adapted from the original Python implementation by C.Y. Fu, here: https://github.com/CrystalCYFu/PyBaySIC/. When using this model, please cite:

Fu, C. Y., Osman, M. B., & Aquino-López, M. A. (2025). Bayesian calibration for the Arctic sea ice biomarker IP<sub>25</sub>. _Paleoceanography and Paleoclimatology_, 40, e2024PA005048. https://doi.org/10.1029/2024PA005048

## Features

- **Nonlinearity**: BaySIC uses an inverse logistic function to characterise the nonlinear relationship between SIC and ln(PIP<sub>25</sub>), respecting the natural limit of SIC between 0 and 1.
- **Bi-directional uncertainty quantification**: Calibration uncertainties are quantified using highest density intervals (HDIs) in the outputs of both the forward and inverse models.
- **Non-stationary seasonality**: The forward model is based on a spatially varying calibration that correlates ln(PIP<sub>25</sub>) with the mean SIC for the three-month interval before the first SIC decrease, accounting for spatiotemporal variations in proxy seasonality.
- **Salinity as an additional environmental driver**: Thresholds have been identified for sea surface salinity below which SIC ceases to be the dominant driver of ln(PIP<sub>25</sub>); we caution against the use of BaySIC in such cases.
- **Direct alignment with PyBaySIC:** Uses identical calibration matrices for full cross-compatibility, which are automatically located (or prompted) by the code.  
- **MATLAB-native plotting capabilities:** Allows for optional PDF and HDI-series visualisations.

For more details, please refer to the source publication.

---

## Requirement

- MATLAB R2022a or newer  

## External dependencies

The following external packages are included in this repository:
- [**cbrewer**](https://uk.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2) - for built-in plotting capabilities.
- [**npy-matlab**](https://github.com/kwikteam/npy-matlab) - to read NumPy .npy posterior matrices.

Before first use, please ensure the calibration files are consistent with the PyBaySIC package.  This can be done automatically by running:

```matlab
download_calib_files
```

All .m functions (lnPIP25_forward.m, lnPIP25_predict.m, etc.) should reside in the same working directory as the calibration files.

## Usage

## Forward model — `lnPIP25_forward.m`

Predict ln(PIP<sub>25</sub>) given sea-ice concentration (SIC) using the forward calibration.

**Inputs**
1. `sic`: scalar or vector of fractional SIC values (0–1)
2. `index`: `'dino'` or `'bras'`
3. (optional) `bayes` — path name to posterior file (`.txt`); default uses built-in calibrations
4. (optional) `plotOpt` — boolean flag (if true, plots PDF or HDI series)

**Output**
`lnpip25` — N × 1000 matrix of posterior draws for ln(PIP<sub>25</sub>)

_Basic usage_:

```matlab
sic = [0.05, 0.50, 0.95];
lnpip25 = lnPIP25_forward(sic, 'dino');
```

_Plot single-posterior PDF (scalar input)_:

```matlab
lnpip25 = lnPIP25_forward(0.35, 'bras', [], true);
```

_Plot multi-point High Density Interval (HDI) series_:

```matlab
sic = linspace(0.05, 0.95, 20);
lnpip25 = lnPIP25_forward(sic, 'dino', [], true);
```

_Use a specific Bayesian posterior file_:

```matlab
lnpip25 = lnPIP25_forward([0.05, 0.50, 0.95], 'bras', 'server_bras_spavar_7900_2024-09-06_08-56-41.txt', true); 
```

## Inverse model — `lnPIP25_predict.m`

Predict MAM SIC from paired biomarker data using the inverse calibration.

**Inputs**
1. `ip25`: vector of IP<sub>25</sub> concentrations (≥0).
2. `sterol`: vector of brassicasterol or dinosterol (≥0) values, in the same units as ip25.
3. `index`: `'dino'` or `'bras'`.
4. `unit`: `'toc'` or `'sed'`
5. (optional) `bayes`: path to posterior matrix (`.npy` or `.mat`); auto-selected if empty
6. (optional) `plotOpt`: boolean flag (if true, plots PDF or HDI series)

**Output**
`sic`: N × 1000 matrix of posterior draws for MAM SIC (0–1).

_Basic usage_:

```matlab
ip25   = rand(20,1) * 0.09; % IP25 (0–0.09 mg/g TOC)
sterol = rand(20,1) * 9.0; % Sterol (0–9 mg/g TOC)
sic = lnPIP25_predict(ip25, sterol, 'dino', 'toc');
```

_Plot PDF for a single observation_:

```matlab
sic = lnPIP25_predict(ip25(1), sterol(1), 'bras', 'toc', [], true); % note scalar inputs for ip25 and sterol
```

_Plot the full HDI time series_:

```matlab
sic = lnPIP25_predict(ip25, sterol, 'bras', 'toc', [], true);
```

_Use a specific Bayesian posterior file_:

```matlab
sic = lnPIP25_predict(ip25, sterol, 'dino', 'sed', 'server_MAM_bras_2024-08-05_09-18-15.npy');
```

## Get in touch

If issues arise when using this Matlab repository, please reach out directly to Matt Osman (mo549@cam.ac.uk).

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/).

**Copyright (c) 2025**
