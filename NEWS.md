# arimatools 0.3.4

Fixes from the v0.3.3 audit.

* Outlier regressors are recognised only by the exact names that
  `outlier_regressors()` creates (type plus period or observation number, e.g.
  `AO2020Q2`, `LS2008M09`, `TC1990`, `AO151`). Ordinary regressors whose names
  start with AO, LS or TC (e.g. `LSE`, `TCU`) are no longer plotted as outlier
  effects by `plot_outliers()`, and the forecasting functions give the correct
  "needs future regressor values" message for them.
* `outlier_regressors()`: `types` must have length 1 or one entry per date.
  Other lengths (e.g. two types for four dates) are an error instead of being
  recycled.
* `forecast_error_variance()`, `dynamic_forecast()`, `plot_forecast()` and
  `plot_fan_chart()` warn when the MA polynomial has a root on or very near
  the unit circle (inverse-root modulus >= 0.99, typically after
  overdifferencing). The psi-weight intervals then understate the uncertainty
  relative to the Kalman filter. `plot_fan_chart()` warns once, not once per
  band.
* `recursive_estimation()` rejects series with missing values inside the
  sample with a clear message. Previously a plain vector with interior gaps
  was silently compressed, joining non-adjacent observations.
* DESCRIPTION gains `URL` and `BugReports` (the GitHub repository).
* API guide: `self-contained-math: true`.

# arimatools 0.3.3

* `compare_models()`'s `sample_id` is now a fingerprint of the complete series
  (an xxhash64 digest from the 'digest' package of every observation in its
  original order, with missing-value positions and time attributes), instead
  of a checksum of summary statistics. Different series with identical sums,
  squared sums and weighted sums, or with different missing-value positions,
  no longer receive the same identifier, so `select_model()` rejects them.
* New dependency: 'digest'.

# arimatools 0.3.2

Fixes the two remaining defects from the v0.3.1 final verification.

* `detect_outliers()` accepts a `drift` regressor only if it is exactly the
  built-in drift of `forecast::Arima(..., include.drift = TRUE)`; shifted or
  rescaled trends named `drift` are rejected instead of being replaced.
* `compare_models()` records each model's estimation sample (`start`, `end`)
  and a fingerprint of its data (`sample_id`). `select_model()` refuses to
  choose among models fitted to different windows (even with equal n) or to
  different data over the same dates (e.g. `y` and `log(y)`).
* README and API guide describe the ATSE software suite, including the
  forthcoming ADL-Tools.

# arimatools 0.3.1

Bug-fix release following the v0.3.0 final validation review. All changes are
guards and documentation; no new statistical methods.

* **Forecasting (F1, F7):** `forecast_error_variance()`, `dynamic_forecast()`,
  `plot_forecast()` and `plot_fan_chart()` reject models with a Box-Cox
  transformation (`lambda`), whose point forecasts and variances would be on
  different scales, and series with missing values, for which the psi-weight
  variances understate the uncertainty. `dynamic_effects()` and
  `plot_dynamic_effects()` also reject `lambda` models.
* **Outlier detection (F2, F8):** argument names in `...` are matched as
  `tso()` matches them, so abbreviations such as `xr` (for `xreg`) are caught;
  unknown arguments are errors. A regressor named `drift` is accepted only if
  it is the built-in linear trend. Models with fixed coefficients are rejected.
  Infinite values are rejected like missing ones.
* **Model selection (F3, F4):** `order_search()` estimates and reports; it
  groups rows by d and sorts by the criterion within each d, with a note when
  several d are present. `select_model()` is the only function that chooses,
  and only among comparable models: it requires `d =` when the table has
  several d, and stops if candidates have different numbers of observations.
  CSS-only fits have no likelihood, so their information criteria are `NA`.
* **Diagnostics (F5, F6, F10):** ACF/PACF, Ljung-Box and the residual
  diagnostic plot reject residuals with gaps inside the sample.
  `recursive_estimation()` rejects `xreg` and `lambda`, and keeps the error
  messages of failed endpoints in `attr(result, "failures")`.
  `ljung_box_test()` rejects fractional lags instead of truncating them.
* **Manual regressors (F9):** `outlier_regressors()` rejects multi-column
  series and empty or missing positions/types.
* **Documentation:** help pages are generated with Markdown enabled (no more
  literal backticks or brackets); forecasts are described as single-origin
  forecasts with intervals conditional on the estimated model; the
  `forecast::Arima` innovation-variance denominator is stated precisely;
  per-function support restrictions are documented; README licence link fixed.

# arimatools 0.3.0

* New outlier functions (R only for now): `outlier_regressors()`,
  `detect_outliers()` (Chen-Liu procedure via `tsoutliers::tso()`) and
  `plot_outliers()`.
* `detect_outliers()` supports univariate, nonseasonal ARIMA models of the
  series as given, with an optional mean or drift. Models with explanatory
  variables (`xreg`) or a Box-Cox transformation (`lambda`) are rejected.
* Warnings from `tso()` and non-convergence of the final model are reported
  (`$warnings`, `$converged`) instead of being silenced.
* Outlier dates are checked against the series frequency (e.g. `"2000Q5"` and
  `"2008M09"` on a quarterly series are errors); weekly and daily series need
  `positions`.
* Clear errors for missing values and samples too short for the model.
* `plot_outliers()` draws one marker per date when several outliers share it.
* Forecasting models with outlier regressors gives a "not supported yet" error.

# arimatools 0.2.0

* Initial R port of the Python package ARIMA-Tools (26 functions).
