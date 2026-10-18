# arimatools (R)

Teaching-oriented R tools for empirical ARIMA modelling, developed for **Applied Time Series Econometrics (ATSE)** at **Universidad Carlos III de Madrid**. This is the R port of the Python package [ARIMA-Tools](https://github.com/maarranz/ARIMA-Tools). It complements `forecast` and `stats` with tables, summaries, numerical forecasts, and `ggplot2` figures (interactive via `plotly`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("maarranz/ARIMA-Tools-R")
library(arimatools)
```

Requires R ≥ 4.1. Imports `digest`, `forecast`, `ggplot2`, `patchwork` and `zoo`; `plotly` is optional (`interactive = TRUE`), and `tsoutliers` is needed only for `detect_outliers()`.

## Quick start

```r
library(forecast)
library(arimatools)

set.seed(12345)
y <- ts(arima.sim(list(ar = 0.7, ma = 0.4), n = 240), start = c(2000, 1), frequency = 12)
fit <- Arima(y, order = c(1, 0, 1), include.mean = FALSE)

coefficient_diagnostics(fit)
root_summary(fit)
ljung_box_test(fit, lags = c(5, 10, 20))
plot_residual_diagnostics(fit)
plot_fan_chart(fit, steps = 20)
```

## Public API

The same 26 functions as the Python package, with the same names, plus three outlier functions (R only for now), exported flat (`arimatools::function()` or `library(arimatools)`).

- **Estimation and model properties:** `coefficient_diagnostics`, `parameter_correlation`, `joint_significance`, `root_diagnostics`, `is_causal`, `is_invertible`, `common_roots`, `root_summary`, `plot_inverse_roots`.
- **Diagnostics:** `ljung_box_test`, `residual_correlations`, `jarque_bera_test`, `plot_residual_diagnostics`, `dynamic_effects`, `plot_dynamic_effects`, `recursive_estimation`, `plot_parameter_stability`.
- **Model comparison and selection:** `compare_models`, `order_search`, `select_model`, `auto_arima`, `auto_arima_summary`. `order_search` and `compare_models` estimate and report; only `select_model` chooses, and only among comparable models (same d, same sample, a genuine likelihood).
- **Forecasting:** `forecast_error_variance`, `dynamic_forecast`, `plot_forecast`, `plot_fan_chart` (default fan coverages `c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99)`).
- **Outliers (AO, LS, TC; nonseasonal models):** `outlier_regressors` (dummies for dates you choose, e.g. `"2020Q2"`), `detect_outliers` (Chen–Liu procedure via `tsoutliers::tso()`, for a fitted model, a fixed order or an auto.arima choice), `plot_outliers` (observed vs adjusted series and the combined outlier effect, for detected or hand-picked outliers). Models with explanatory variables or a Box–Cox transformation are rejected; warnings and non-convergence are reported.

[API documentation — rendered HTML](docs/index.html) · [Quarto source](docs/arimatools_Documentation.qmd) · in R: `?recursive_estimation`, etc.

## Differences from the Python version

The port is idiomatic R rather than a line-by-line translation.

| Topic | Python (ARIMA-Tools) | R (arimatools) |
|:--|:--|:--|
| Fitted model | Statsmodels `ARIMAResults` | `forecast::Arima()`, `auto.arima()` or `stats::arima()` |
| `auto_arima()` result | StatsForecast object, not accepted by other functions | ordinary `Arima` object, accepted everywhere |
| Coefficient names | `ar.L1`, `ma.L1`, `const`, `sigma2` | `ar1`, `ma1`, `intercept`, `drift` |
| Deterministic terms | `trend=` | `include.constant=` (and `...` to `forecast::Arima`) |
| Table columns | `Estimate`, `Std. Error`, `P-value`, ... | `estimate`, `std.error`, `p.value`, ... (broom style) |
| Single-result outputs | pandas `Series` | named list that prints as two columns; `as.data.frame()` |
| Figures | Plotly `Figure` | `ggplot`/`patchwork`; `interactive = TRUE` → plotly |
| Figure size arguments | `width`, `height` | chunk options (`fig-width`, `fig-height`) |
| `parameter_correlation(include_variance=)` | available | dropped: R has no covariance for `sigma2` |
| Standard errors | OPG (Statsmodels default) | Hessian (= Statsmodels `cov_type="approx"`) |
| `sigma2` | ML estimate, SSR/n | `forecast::Arima`: SSR/(n′ − #coefficients); `stats::arima`: SSR/n′ (n′ = observations after differencing) |
| Ljung–Box df | h − p − q | h − p − q − P − Q (identical for nonseasonal models) |
| `recursive_estimation()` | | also stores the recursive `sigma2`; retries with exact ML if CSS-ML fails |

Point estimates, log likelihoods, information criteria, roots, psi weights and point forecasts agree with the Python package to about 1e-5 on the same data (exact ML). Standard errors, `sigma2` and residual-based statistics differ slightly for the reasons in the table.

## Tests

```r
testthat::test_local()   # from the package root
```

The suite ports `tests/test_forecasting.py` (FEV against Kalman variances and closed forms for AR(1), ARMA(1,1) and the random walk) and adds checks for every other function against base R, `forecast`, `lmtest` and `tsoutliers`.

## The ATSE software suite

arimatools is one of a set of companion packages for the ATSE course:

- **[ARIMA-Tools](https://github.com/maarranz/ARIMA-Tools)** (Python) — the original package; arimatools is its R version.
- **[ARMA-Theory](https://github.com/maarranz/ARMA-theory)** — theoretical ARMA properties, roots, ACF/PACF and impulse responses.
- **[ADL-Theory](https://github.com/maarranz/adl-theory)** — theoretical autoregressive distributed lag (ADL) models.
- **[Forecasting-Tools](https://github.com/maarranz/Forecasting-Tools)** — forecast evaluation and comparison (Python; R version planned).
- **ADL-Tools** — empirical ADL modelling (forthcoming).

## License

[MIT License](LICENSE.md).

## Author

**Miguel A. Arranz, PhD**\
Department of Economics\
Universidad Carlos III de Madrid (UC3M)\
Email: [maarranz@eco.uc3m.es](mailto:maarranz@eco.uc3m.es)\
ORCID: [https://orcid.org/0000-0002-5951-2284](https://orcid.org/0000-0002-5951-2284)
