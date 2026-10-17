skip_tso <- function() skip_if_not_installed("tsoutliers")

set.seed(3)
y <- ts(as.numeric(arima.sim(list(ar = 0.6), 120)), start = c(1990, 1), frequency = 4)
y[40] <- y[40] + 6          # AO at 1999Q4
y[70:120] <- y[70:120] + 5  # LS at 2007Q2

test_that("outlier_regressors builds AO, LS and TC columns", {
  X <- outlier_regressors(y, dates = c("1999Q4", "2007:2", "2010-08-15"),
                          types = c("AO", "LS", "TC"), delta = 0.7)
  expect_equal(colnames(X), c("AO1999Q4", "LS2007Q2", "TC2010Q3"))
  expect_equal(tsp(X), tsp(y))
  expect_equal(which(X[, 1] != 0), 40L)
  expect_equal(which(X[, 2] != 0), 70:120)
  expect_equal(X[83:86, 3], c(1, 0.7, 0.49, 0.343))
  expect_equal(sum(X[1:82, 3]), 0)
  expect_equal(attr(X, "outliers")$index, c(40L, 70L, 83L))
})

test_that("regressors match tsoutliers::outliers.effects", {
  skip_tso()
  X <- outlier_regressors(y, positions = c(10, 20, 30), types = c("TC", "LS", "AO"))
  ref <- tsoutliers::outliers.effects(
    data.frame(type = c("TC", "LS", "AO"), ind = c(10, 20, 30)), length(y), delta = 0.7)
  expect_equal(unclass(X)[, 1:3], unname(ref), ignore_attr = TRUE)
})

test_that("dates are parsed in several formats and checked", {
  expect_equal(colnames(outlier_regressors(y, dates = as.Date("2007-05-20"))), "AO2007Q2")
  expect_equal(colnames(outlier_regressors(y, dates = 2007.25)), "AO2007Q2")
  m <- ts(rnorm(60), start = c(2005, 1), frequency = 12)
  expect_equal(colnames(outlier_regressors(m, dates = "2008M9", types = "LS")), "LS2008M09")
  a <- ts(rnorm(30), start = 1990)
  expect_equal(colnames(outlier_regressors(a, dates = 2001)), "AO2001")
  expect_equal(colnames(outlier_regressors(rnorm(50), positions = 12)), "AO12")
  expect_error(outlier_regressors(y, dates = "2030Q1"), "outside the sample")
  expect_error(outlier_regressors(y, dates = "1999Q4", types = "XX"), "types")
  expect_error(outlier_regressors(y), "exactly one")
  expect_error(outlier_regressors(y, dates = c("1999Q4", "1999Q4")), "twice")
})

test_that("manual dummies work with Arima, diagnostics and plot_outliers", {
  X <- outlier_regressors(y, dates = c("1999Q4", "2007Q2"), types = c("AO", "LS"))
  fit <- forecast::Arima(y, order = c(1, 0, 0), xreg = X)
  cd <- coefficient_diagnostics(fit)
  expect_true(all(colnames(X) %in% cd$term))
  expect_equal(cd$estimate[cd$term == "LS2007Q2"], 5, tolerance = 0.15)
  expect_s3_class(plot_outliers(fit), "patchwork")
  expect_error(plot_outliers(forecast::Arima(y, order = c(1, 0, 0))), "no outlier")
  expect_error(dynamic_forecast(fit), "not supported yet")
})

test_that("detect_outliers reproduces tso with a fixed order", {
  skip_tso()
  fit0 <- forecast::Arima(y, order = c(1, 0, 0))
  res <- detect_outliers(fit0)
  ref <- tsoutliers::tso(y, types = c("AO", "LS", "TC"), tsmethod = "arima",
                         args.tsmethod = list(order = c(1, 0, 0)), check.rank = TRUE)
  expect_s3_class(res, "arimatools_outliers")
  expect_equal(res$outliers$type, as.character(ref$outliers$type))
  expect_equal(res$outliers$index, ref$outliers$ind)
  expect_equal(res$outliers$date, c("1999Q4", "2007Q2"))
  expect_equal(res$outliers$estimate, ref$outliers$coefhat, tolerance = 1e-4)
  expect_equal(res$outliers$statistic, ref$outliers$tstat, tolerance = 1e-3)
  expect_equal(as.numeric(res$adjusted + res$effects), as.numeric(y))
  expect_s3_class(res$fit, "Arima")
  expect_equal(nrow(ljung_box_test(res$fit, 8)), 1L)
  expect_s3_class(plot_outliers(res), "patchwork")
  expect_output(print(res), "Outliers found:  2")
  # same result from y + order
  res2 <- detect_outliers(y, order = c(1, 0, 0))
  expect_equal(res2$outliers$index, res$outliers$index)
})

test_that("detect_outliers: auto.arima path, drift, no outliers, errors", {
  skip_tso()
  res <- detect_outliers(y)
  expect_true(res$auto)
  expect_true(all(c(40L, 70L) %in% res$outliers$index))
  set.seed(5)
  yd <- ts(cumsum(0.3 + rnorm(150)), start = c(1980, 1), frequency = 4)
  yd[100:150] <- yd[100:150] - 6
  rd <- detect_outliers(forecast::Arima(yd, order = c(0, 1, 1), include.drift = TRUE))
  expect_true("drift" %in% names(rd$fit$coef))
  expect_true(100L %in% rd$outliers$index)
  set.seed(9)
  r0 <- detect_outliers(ts(rnorm(100)), order = c(0, 0, 0), cval = 5)
  expect_equal(nrow(r0$outliers), 0L)
  expect_true(all(r0$effects == 0))
  expect_s3_class(plot_outliers(r0), "patchwork")
  seas <- forecast::Arima(log(AirPassengers), order = c(0, 1, 1), seasonal = c(0, 1, 1))
  expect_error(detect_outliers(seas), "nonseasonal")
  expect_error(detect_outliers(y, types = "SLS"), "types")
})

# ---- Regression tests for the numerical-validation review (D1-D7) ----------

test_that("D1: explanatory variables are rejected", {
  skip_tso()
  set.seed(1701)
  Z <- cbind(economic_driver = rnorm(120))
  yz <- ts(3 * Z[, 1] + rnorm(120, sd = 0.3), start = c(1980, 1), frequency = 4)
  expect_error(detect_outliers(yz, order = c(0, 0, 0), xreg = Z), "not supported")
  expect_error(detect_outliers(stats::arima(yz, order = c(0, 0, 0), xreg = Z)),
               "explanatory variables")
  expect_error(detect_outliers(forecast::Arima(yz, order = c(0, 0, 0), xreg = Z)),
               "explanatory variables")
  expect_error(detect_outliers(yz, order = c(0, 0, 0), tsmethod = "arima"),
               "cannot be passed")
})

test_that("D2: Box-Cox models are rejected", {
  skip_tso()
  set.seed(1702)
  yb <- ts(exp(1 + rnorm(120, sd = 0.15)), start = c(1980, 1), frequency = 4)
  f <- forecast::Arima(yb, order = c(0, 0, 0), lambda = 0)
  expect_error(detect_outliers(f), "Box-Cox")
  X <- outlier_regressors(yb, positions = 60)
  fm <- forecast::Arima(yb, order = c(0, 0, 0), lambda = 0, xreg = X)
  expect_error(plot_outliers(fm), "Box-Cox")
  expect_error(detect_outliers(yb, order = c(0, 0, 0), lambda = 0), "not supported")
})

test_that("D3: non-convergence and tso warnings are reported", {
  skip_tso()
  set.seed(1710)
  yn <- ts(rnorm(40))
  w <- testthat::capture_warnings(r <- detect_outliers(yn, order = c(4, 0, 4),
                                                       cval = 100))
  expect_false(r$converged)
  expect_true(any(grepl("did not converge", w)))
  expect_type(r$warnings, "character")
  expect_output(print(r), "Converged:       NO")
})

test_that("D4: invalid periods and unsupported calendars are rejected", {
  z <- ts(rep(0, 100), start = 2000, frequency = 4)
  expect_error(outlier_regressors(z, dates = "2000Q5"), "between 1 and 4")
  expect_error(outlier_regressors(z, dates = "2000:5"), "between 1 and 4")
  expect_error(outlier_regressors(z, dates = "2000M09"), "monthly date")
  expect_error(outlier_regressors(z, dates = 2005), "ambiguous")
  expect_error(outlier_regressors(z, dates = 2005.1), "not matching")
  m <- ts(rep(0, 60), start = c(2005, 1), frequency = 12)
  expect_error(outlier_regressors(m, dates = "2006M13"), "between 1 and 12")
  expect_error(outlier_regressors(m, dates = "2006Q1"), "quarterly date")
  expect_error(outlier_regressors(m, dates = "2006-02-30"), "Invalid date")
  daily <- ts(rep(0, 366), start = 2020, frequency = 366)
  expect_error(outlier_regressors(daily, dates = as.Date("2020-01-15")), "positions")
  expect_equal(attr(outlier_regressors(daily, positions = 15), "outliers")$index, 15L)
})

test_that("D5: same-date interior LS and TC are distinct; last-date ones coincide", {
  set.seed(1703)
  z <- ts(rnorm(100), start = 2000, frequency = 4)
  X <- outlier_regressors(z, positions = c(50, 50), types = c("LS", "TC"))
  expect_equal(qr(X)$rank, 2L)
  Xend <- outlier_regressors(z, positions = c(100, 100), types = c("LS", "TC"))
  expect_equal(qr(Xend)$rank, 1L)
})

test_that("D6: missing values and short samples give clear errors", {
  skip_tso()
  set.seed(1703)
  z <- ts(rnorm(100), start = 2000, frequency = 4)
  z[30] <- NA
  expect_error(detect_outliers(z, order = c(1, 0, 0)), "missing values")
  expect_error(detect_outliers(forecast::Arima(z, order = c(1, 0, 0))), "missing values")
  expect_error(detect_outliers(ts(rnorm(8)), order = c(0, 0, 0)), "Too few")
  expect_error(detect_outliers(ts(rnorm(15))), "Too few")
  expect_error(detect_outliers(cbind(a = rnorm(50), b = rnorm(50))), "univariate")
})

test_that("D7: outliers on the same date share one marker", {
  set.seed(1)
  z <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
  X <- outlier_regressors(z, positions = c(40, 40, 60), types = c("AO", "LS", "TC"))
  fit <- forecast::Arima(z, order = c(1, 0, 0), xreg = X)
  p <- plot_outliers(fit)
  pts <- p[[1]]$layers[[2]]$data
  expect_equal(nrow(pts), 2L)
  expect_equal(pts$type, c("Several", "TC"))
  expect_match(pts$text[1], "AO2009Q4.*LS2009Q4")
})
