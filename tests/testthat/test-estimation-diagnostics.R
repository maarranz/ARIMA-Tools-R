set.seed(12345)
y <- arima.sim(list(ar = 0.7, ma = 0.4), n = 600)
fit <- forecast::Arima(y, order = c(1, 0, 1), include.mean = FALSE)

test_that("roots solve the AR and MA polynomials", {
  rd <- root_diagnostics(fit)
  expect_equal(rd$type, c("AR", "MA"))
  expect_equal(rd$modulus[rd$type == "AR"], 1 / abs(fit$coef[["ar1"]]))
  expect_equal(rd$modulus[rd$type == "MA"], 1 / abs(fit$coef[["ma1"]]))
  set.seed(3); y2 <- arima.sim(list(ar = c(1.2, -0.5)), n = 500)
  f2 <- forecast::Arima(y2, order = c(2, 0, 0), include.mean = FALSE)
  z <- root_diagnostics(f2)$value
  phi <- f2$coef
  expect_equal(Mod(1 - phi[1] * z - phi[2] * z^2), c(0, 0), tolerance = 1e-8)
  expect_equal(rd$inverse_modulus, 1 / rd$modulus)
  expect_true(is_causal(fit))
  expect_true(is_invertible(fit))
})

test_that("common roots are detected and classified", {
  set.seed(2468)
  yn <- arima.sim(list(ar = 0.70, ma = -0.68), n = 2000)
  fn <- forecast::Arima(yn, order = c(1, 0, 1), include.mean = FALSE)
  cr <- common_roots(fn, tol = 2)
  expect_equal(nrow(cr), 1L)
  expect_equal(cr$distance, Mod(cr$ar_inverse_root - cr$ma_inverse_root))
  expect_equal(nrow(common_roots(fit)), 0L)
  rs <- root_summary(fn, tol = 2)
  expect_s3_class(rs, "arimatools_summary")
  expect_equal(rs[["Near-common root pairs"]] + rs[["Common root pairs"]], 1)
  expect_s3_class(plot_inverse_roots(fn, tol = 2), "ggplot")
})

test_that("pure AR and white-noise models have empty MA root sets", {
  f0 <- forecast::Arima(y, order = c(0, 0, 0))
  expect_equal(nrow(root_diagnostics(f0)), 0L)
  expect_true(is_causal(f0) && is_invertible(f0))
  expect_equal(dynamic_effects(f0, 3)$effect, c(1, 0, 0, 0))
})

test_that("coefficient inference matches lmtest::coeftest and confint", {
  cd <- coefficient_diagnostics(fit, confidence = 0.90)
  ct <- lmtest::coeftest(fit)
  expect_equal(cd$estimate, unname(ct[, 1]))
  expect_equal(cd$std.error, unname(ct[, 2]))
  expect_equal(cd$p.value, unname(ct[, 4]))
  expect_equal(cd$conf.low, unname(confint(fit, level = 0.90)[, 1]))
  expect_equal(parameter_correlation(fit), cov2cor(vcov(fit)))
})

test_that("Wald test of one coefficient equals squared z", {
  js <- joint_significance(fit, "ma1")
  z <- coefficient_diagnostics(fit)$statistic[2]
  expect_equal(js[["Wald statistic"]], z^2)
  expect_error(joint_significance(fit, "ar9"), "Unknown")
  expect_error(joint_significance(fit, c("ar1", "ar1")), "repeated")
})

test_that("Ljung-Box uses df = h - p - q and rejects small lags", {
  lb <- ljung_box_test(fit, c(5, 10))
  expect_equal(lb$df, c(3L, 8L))
  bt <- Box.test(residuals(fit), lag = 10, type = "Ljung-Box", fitdf = 2)
  expect_equal(lb$statistic[2], unname(bt$statistic))
  expect_error(ljung_box_test(fit, 2), "greater than")
})

test_that("residual correlations and Jarque-Bera", {
  rc <- residual_correlations(fit, lags = 12)
  expect_equal(nrow(rc), 12L)
  expect_equal(rc$acf, acf(residuals(fit), 12, plot = FALSE)$acf[-1])
  expect_equal(rc$upper, rep(qnorm(0.975) / sqrt(600), 12))
  jb <- jarque_bera_test(fit)
  r <- as.numeric(residuals(fit)); n <- length(r); m <- r - mean(r)
  s <- mean(m^3) / mean(m^2)^1.5; k <- mean(m^4) / mean(m^2)^2
  expect_equal(jb[["Jarque-Bera statistic"]], n / 6 * (s^2 + (k - 3)^2 / 4))
})

test_that("dynamic effects include integration", {
  expect_equal(dynamic_effects(fit, 5)$effect, c(1, ARMAtoMA(0.7, 0.4, 5)),
               tolerance = 0.2)
  expect_equal(dynamic_effects(fit, 5)$effect,
               c(1, ARMAtoMA(fit$coef[1], fit$coef[2], 5)))
  rw <- forecast::Arima(cumsum(rnorm(100)), order = c(0, 1, 0))
  expect_equal(dynamic_effects(rw, 4)$effect, rep(1, 5))
  expect_equal(dynamic_effects(rw, 4)$cumulative, 1:5)
})

test_that("diagnostic plots build", {
  expect_s3_class(plot_residual_diagnostics(fit), "patchwork")
  expect_s3_class(plot_dynamic_effects(fit), "patchwork")
})

test_that("recursive estimation: windows, columns and plot", {
  rec <- recursive_estimation(y, c(1, 0, 1), initial = 500)
  expect_named(rec, c("endpoint", "n", "term", "estimate", "std.error",
                      "conf.low", "conf.high", "converged"))
  expect_equal(range(rec$n), c(500, 600))
  last <- rec[rec$n == 600 & rec$term == "ar1", ]
  full <- forecast::Arima(y, order = c(1, 0, 1))
  expect_equal(last$estimate, unname(full$coef["ar1"]), tolerance = 1e-6)
  expect_equal(nrow(recursive_estimation(y, c(1, 0, 1), initial = 0.9)),
               4L * 61L)
  expect_error(recursive_estimation(y, c(1, 0, 1), initial = 600), "too large")
  expect_error(recursive_estimation(y, c(1, 0, 1), initial = 4), "too small")
  p <- plot_parameter_stability(rec, breakpoint = 550)
  expect_s3_class(p, "ggplot")
  expect_error(plot_parameter_stability(rec, "ar7"), "Unknown")
})


test_that("F1: Box-Cox fits are rejected by dynamic effects", {
  set.seed(301)
  yb <- ts(exp(3 + as.numeric(arima.sim(list(ar = 0.6), n = 120, sd = 0.3))),
           start = 1990, frequency = 4)
  fb <- forecast::Arima(yb, order = c(1, 0, 0), lambda = 0)
  expect_error(dynamic_effects(fb), "Box-Cox")
  expect_error(plot_dynamic_effects(fb), "Box-Cox")
})

test_that("F5: residual gaps are rejected by autocorrelation diagnostics", {
  set.seed(305)
  yg <- ts(as.numeric(arima.sim(list(ar = 0.6), n = 100)))
  yg[seq(10, 90, 10)] <- NA
  fg <- forecast::Arima(yg, order = c(1, 0, 0))
  expect_error(residual_correlations(fg, lags = 5), "inside the sample")
  expect_error(ljung_box_test(fg, 5), "inside the sample")
  expect_error(plot_residual_diagnostics(fg), "inside the sample")
  expect_s3_class(jarque_bera_test(fg), "arimatools_summary")
})

test_that("F6: recursive estimation rejects xreg and keeps failure reasons", {
  set.seed(304)
  Z <- cbind(z = rnorm(80))
  yz <- ts(0.6 * Z[, 1] + rnorm(80))
  expect_error(recursive_estimation(yz, c(0, 0, 0), initial = 60, xreg = Z), "xreg")
  expect_error(recursive_estimation(yz, c(0, 0, 0), initial = 60, xr = Z), "xreg")
  expect_error(recursive_estimation(yz, c(0, 0, 0), initial = 60, lambda = 0), "lambda")
  rec <- recursive_estimation(y, c(1, 0, 1), initial = 590)
  expect_s3_class(attr(rec, "failures"), "data.frame")
  expect_equal(nrow(attr(rec, "failures")), 0L)
})

test_that("F10: fractional or invalid Ljung-Box lags are rejected", {
  expect_error(ljung_box_test(fit, 10.9), "positive integers")
  expect_error(ljung_box_test(fit, c(5, NA)), "positive integers")
  expect_error(ljung_box_test(fit, "10"), "positive integers")
})
