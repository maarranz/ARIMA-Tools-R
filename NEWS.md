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
