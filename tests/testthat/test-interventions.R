# v0.5.0: forecasting with regressors, extend_interventions(),
# intervention_effects(), plot_interventions()

make_case <- function() {
  set.seed(1)
  y <- ts(cumsum(0.2 + rnorm(120)), start = c(1990, 1), frequency = 4)
  y[60:120] <- y[60:120] + 4
  y[100] <- y[100] + 3
  X <- outlier_regressors(y, dates = c("2004Q4", "2014Q4", "2012Q1"),
                          types = c("LS", "AO", "TC"), delta = 0.6)
  fit <- forecast::Arima(y, order = c(1, 1, 0), xreg = X, include.drift = TRUE)
  list(y = y, X = X, fit = fit)
}

test_that("extend_interventions continues AO, LS and TC from the event dates", {
  cs <- make_case()
  F <- extend_interventions(cs$X, 8)
  expect_equal(dim(F), c(8L, 3L))
  expect_equal(colnames(F), colnames(cs$X))
  expect_true(stats::is.ts(F))
  expect_equal(stats::tsp(F)[1], 2020)
  expect_equal(as.numeric(F[, "AO2014Q4"]), rep(0, 8))
  expect_equal(as.numeric(F[, "LS2004Q4"]), rep(1, 8))
  # TC at overall observation 121..128, event at 89: 0.6^(32..39)
  expect_equal(as.numeric(F[, "TC2012Q1"]), 0.6^(32:39))
  # Chaining keeps the event ages
  expect_equal(unclass(extend_interventions(F, 4))[, 3],
               unclass(extend_interventions(cs$X, 12))[9:12, 3],
               ignore_attr = TRUE)
  # Plain vectors give plain matrices
  Xn <- outlier_regressors(as.numeric(cs$y), positions = c(10, 20), types = c("AO", "LS"))
  Fn <- extend_interventions(Xn, 3)
  expect_false(stats::is.ts(Fn))
  expect_equal(unname(Fn[, 2]), rep(1, 3))
})

test_that("extend_interventions rejects modified or foreign matrices", {
  cs <- make_case()
  expect_error(extend_interventions(cs$X[, 1:2], 4), "created by outlier_regressors")
  Z <- cs$X; Z[5, 1] <- 1
  expect_error(extend_interventions(Z, 4), "no longer match")
  expect_error(extend_interventions(matrix(0, 5, 1), 4), "created by outlier_regressors")
  expect_error(extend_interventions(cs$X, 0), "positive integer")
})

test_that("dynamic_forecast with xreg equals forecast::forecast and Kalman intervals", {
  cs <- make_case()
  F <- extend_interventions(cs$X, 8)
  d <- dynamic_forecast(cs$fit, 8, xreg = F)
  k <- forecast::forecast(cs$fit, h = 8, xreg = F, level = 95)
  expect_equal(d$forecast, as.numeric(k$mean))
  expect_equal(d$upper, as.numeric(k$upper[, 1]), tolerance = 1e-8)
  # Column order does not matter
  expect_equal(dynamic_forecast(cs$fit, 8, xreg = F[, c(3, 1, 2)])$forecast, d$forecast)
  # Plots accept xreg
  expect_s3_class(plot_forecast(cs$fit, 8, xreg = F), "ggplot")
  expect_s3_class(plot_fan_chart(cs$fit, 8, xreg = F), "ggplot")
})

test_that("ordinary regressors are supported too", {
  set.seed(5)
  x <- rnorm(150)
  y <- 2 * x + as.numeric(arima.sim(list(ar = 0.5), n = 150))
  fit <- forecast::Arima(y, order = c(1, 0, 0), xreg = cbind(TCU = x))
  newx <- cbind(TCU = c(0.5, -1, 2))
  expect_equal(dynamic_forecast(fit, 3, xreg = newx)$forecast,
               as.numeric(forecast::forecast(fit, h = 3, xreg = newx)$mean))
  # A single unnamed vector works for a single regressor
  expect_equal(dynamic_forecast(fit, 3, xreg = c(0.5, -1, 2))$forecast,
               dynamic_forecast(fit, 3, xreg = newx)$forecast)
  expect_error(dynamic_forecast(fit, 3), "future values with `xreg =`")
})

test_that("invalid future regressors give clear errors", {
  cs <- make_case()
  F <- extend_interventions(cs$X, 8)
  expect_error(dynamic_forecast(cs$fit, 8), "extend_interventions")
  expect_error(dynamic_forecast(cs$fit, 6, xreg = F), "8 rows but 6")
  G <- unclass(F)[, 1:2]
  expect_error(dynamic_forecast(cs$fit, 8, xreg = G), "missing: TC2012Q1")
  H <- F; H[, "LS2004Q4"] <- 0
  expect_error(dynamic_forecast(cs$fit, 8, xreg = H), "level shift is permanent")
  H <- F; H[, "AO2014Q4"] <- 1
  expect_error(dynamic_forecast(cs$fit, 8, xreg = H), "happened once")
  H <- F; H[, "TC2012Q1"] <- 0.5
  expect_error(dynamic_forecast(cs$fit, 8, xreg = H), "geometric decay")
  late <- stats::ts(unclass(F), start = c(2030, 1), frequency = 4)
  expect_error(dynamic_forecast(cs$fit, 8, xreg = late), "after the end of the sample")
  plain <- forecast::Arima(cs$y, order = c(1, 1, 0))
  expect_error(dynamic_forecast(plain, 8, xreg = F), "not needed")
  M <- unclass(F); colnames(M) <- NULL
  expect_error(dynamic_forecast(cs$fit, 8, xreg = M), "column names")
})

test_that("intervention_effects reproduces beta * X and the adjusted series", {
  cs <- make_case()
  eff <- intervention_effects(cs$fit, cs$X)
  b <- cs$fit$coef[colnames(cs$X)]
  expect_named(eff, c("time", "observed", colnames(cs$X), "combined", "adjusted"))
  expect_equal(eff$LS2004Q4, as.numeric(cs$X[, 1]) * b[[1]])
  expect_equal(eff$combined, as.numeric(unclass(cs$X) %*% b))
  expect_equal(eff$adjusted, as.numeric(cs$y) - eff$combined)
  expect_equal(attr(eff, "coefficients"), b)
  expect_true(inherits(eff$time, "Date"))
  # Errors
  expect_error(intervention_effects(cs$fit, extend_interventions(cs$X, 4)),
               "original sample")
  plain <- forecast::Arima(cs$y, order = c(1, 1, 0))
  expect_error(intervention_effects(plain, cs$X), "not fitted with")
  expect_error(intervention_effects(cs$fit), "Supply `X`")
})

test_that("intervention_effects and forecasting work with detect_outliers output", {
  skip_if_not_installed("tsoutliers")
  cs <- make_case()
  found <- suppressWarnings(detect_outliers(cs$y, order = c(1, 1, 0), include.constant = TRUE))
  skip_if(is.null(found$regressors))
  eff <- intervention_effects(found)
  expect_equal(eff$combined, as.numeric(found$effects))
  h <- 4
  d <- dynamic_forecast(found$fit, h, xreg = extend_interventions(found$regressors, h))
  expect_equal(nrow(d), h)
})

test_that("plot_interventions with and without a fit; plot_outliers unchanged", {
  cs <- make_case()
  p0 <- plot_interventions(cs$y, cs$X)
  expect_s3_class(p0, "ggplot")
  expect_silent(ggplot2::ggplot_build(p0))
  p1 <- plot_interventions(cs$y, cs$X, cs$fit, individual = TRUE)
  expect_s3_class(p1, "patchwork")
  # Combined effect line equals intervention_effects()$combined
  eff <- intervention_effects(cs$fit, cs$X)
  expect_equal(ggplot2::layer_data(p1[[2]], 2)$y, eff$combined)
  ld3 <- ggplot2::layer_data(p1[[2]], 3)          # individual effects (dashed)
  expect_equal(length(unique(ld3$group)), 3L)
  expect_warning(plot_interventions(cs$y, cs$X, individual = TRUE), "ignored")
  expect_error(plot_interventions(cs$y + 1, cs$X, cs$fit), "differs")
  expect_error(plot_interventions(cs$y[1:50], cs$X), "same number")
  expect_s3_class(plot_outliers(cs$fit), "patchwork")
})
