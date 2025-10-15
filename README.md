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
1. `sic`: scalar or vector of fractional SIC values (0–1) -- see also `calc_meanSIC.m`
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

## Sea-ice climatology utility — `calc_meanSIC.m`

Compute the mean SIC for the three months ending at the first post-maximum decrease in sea ice concentration. _This mirrors the seasonal logic used by the BaySIC forward model, and can be used to generate the `sic` input for `lnPIP25_forward.m`._

**Inputs**
1. `sic_climo`:
*Vector mode*: a 12×1 or 1×12 climatology of SIC values for a site location (fractions in [0,1]). *Grid mode*: a 3-D array of SIC climatological values with one dimension having length 12
2. (Optional, used for *Grid Mode* only) `site_lat`, `site_lon`, `all_lat`, `all_lon`: scalars (first two) and arrays (last two) of coordinates, each provided in degrees. Longitude convention (0--360˚ or −180--180˚) is handled automatically; inputs for `all_lat` and `all_lon` can be provided either as vectors or matrices, but must match the dimensionality of `sic_climo`.

**Outputs**
`meanSIC`: scalar mean SIC value over [m−2, m−1, m], where m is the month of first SIC decrease
`monthsUsed`: a 1 × 3 list of of integer months (1–12; Jan=1, Dec=12)
`meta`: a diagnostic structure with fields that include computation flags, grid indices, lat/lon, distance to nearest decreasing SIC value, and applied SIC series)

_Usage_:

_Basic 12-month SIC climatology (fractions 0–1)_:

```matlab
sic = [0.6, 0.7, 0.8, 0.9, 1, 0.9, 0.4, 0.2, 0.1, 0.3, 0.4, 0.5]';
[meanSIC, monthsUsed, meta] = calc_meanSIC(sic);
```

_Grid example using 3D (curvilinear gridded) climatology provided in `example_sic_data`_

```matlab
load example_sic_data/curvilinear_sic_climo.mat  % provides SIC, lat (1D), lon (1D)
site_lat = lat(166); site_lon = lon(36);
[meanSIC, monthsUsed, meta] = calc_meanSIC(SIC./100, site_lat, site_lon, lat, lon); % SIC converted from units of percent to fraction
% visualise it!
figure; hold on
plot(1:12, squeeze(SIC(36,166,:))./100, '-o', 'DisplayName','Input site');
scatter(monthsUsed, squeeze(SIC(36,166,monthsUsed))./100, 60, 'filled', 'DisplayName','Months used'); % could alternatively use meta.sic(monthsUsed) for input y-variable in this example!
xlabel('Month'); ylabel('SIC (fraction)');
```

_Grid example using 3D (tripolar gridded) climatology provided in `example_sic_data`, with search required for nearest decreasing SIC location_

```matlab
load example_sic_data/tripolar_sic_climo.mat  % provides SIC, lat (2D), lon (2D)
site_lat = lat(74,359); site_lon = lon(74,359); % location with near-stationary SIC values
[meanSIC, monthsUsed, meta] = calc_meanSIC(SIC./100, site_lat, site_lon, lat, lon);
% visualise it!
figure; hold on
plot(1:12, squeeze(SIC(74,359,:))./100, '-o', 'DisplayName','Input site');
plot(1:12, meta.sic, '-s', 'DisplayName', sprintf('Nearest variable (%.1f km)', meta.distance));
scatter(monthsUsed, meta.sic(monthsUsed), 60, 'filled', 'DisplayName','Months used');
xlabel('Month'); ylabel('SIC (fraction)');
```


## Get in touch

If issues arise when using this Matlab repository, please reach out directly to Matt Osman (mo549@cam.ac.uk).

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/).

**Copyright (c) 2025**
