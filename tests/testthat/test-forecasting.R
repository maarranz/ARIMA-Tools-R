# Controlled numerical checks for the forecasting functions (port of
# tests/test_forecasting.py).

set.seed(20261005)
models <- list(
  "AR(1)" = forecast::Arima(arima.sim(list(ar = 0.65), n = 1200, n.start = 300),
                            order = c(1, 0, 0), include.mean = FALSE),
  "ARMA(1,1)" = forecast::Arima(arima.sim(list(ar = 0.55, ma = 0.35), n = 1200,
                                          n.start = 300),
                                order = c(1, 0, 1), include.mean = FALSE),
  "Random walk" = forecast::Arima(cumsum(rnorm(1200)), order = c(0, 1, 0))
)

test_that("FEV matches Kalman variances and closed forms", {
  for (nm in names(models)) {
    fit <- models[[nm]]
    s2 <- fit$sigma2
    h <- 1:1000
    fev <- forecast_error_variance(fit, steps = 1000)
    kalman <- predict(fit, n.ahead = 1000)$se^2
    expect_equal(fev, as.numeric(kalman), tolerance = 1e-8, label = nm)
    if (nm == "Random walk") {
      expect_equal(fev, h * s2, tolerance = 1e-12)
    } else {
      phi <- fit$coef[["ar1"]]
      if (nm == "AR(1)") {
        expect_equal(fev, s2 * (1 - phi^(2 * h)) / (1 - phi^2), tolerance = 1e-12)
        limit <- s2 / (1 - phi^2)
      } else {
        theta <- fit$coef[["ma1"]]
        limit <- s2 * (1 + (phi + theta)^2 / (1 - phi^2))
      }
      expect_equal(fev[1000], limit, tolerance = 1e-12)
      expect_gte(fev[1000], fev[1])
    }
  }
})

test_that("forecast table and intervals are coherent", {
  for (fit in models) {
    tab <- dynamic_forecast(fit, steps = 12, level = 0.90)
    expect_named(tab, c("horizon", "forecast", "fev", "se", "lower", "upper"))
    expect_equal(tab$horizon, 1:12)
    expect_equal(tab$forecast, as.numeric(forecast::forecast(fit, h = 12)$mean))
    expect_equal(tab$fev, forecast_error_variance(fit, 12))
    expect_equal(tab$se, sqrt(tab$fev))
    expect_equal(tab$lower, tab$forecast - qnorm(0.95) * tab$se)
    expect_equal(tab$upper, tab$forecast + qnorm(0.95) * tab$se)
  }
})

test_that("seasonal and drift models are handled", {
  y <- log(AirPassengers)
  fit <- forecast::Arima(y, order = c(0, 1, 1), seasonal = c(0, 1, 1))
  expect_equal(forecast_error_variance(fit, 36),
               as.numeric(predict(fit, n.ahead = 36)$se^2), tolerance = 1e-4)
  set.seed(1)
  fd <- forecast::Arima(cumsum(0.3 + rnorm(200)), order = c(1, 1, 0),
                        include.drift = TRUE)
  tab <- dynamic_forecast(fd, 10)
  expect_equal(tab$forecast, as.numeric(forecast::forecast(fd, h = 10)$mean))
})

test_that("invalid arguments are rejected", {
  fit <- models[["AR(1)"]]
  expect_error(forecast_error_variance(fit, 0), "positive integer")
  expect_error(forecast_error_variance(fit, 2.5), "positive integer")
  expect_error(dynamic_forecast(fit, 5, level = 1), "strictly between")
  expect_error(plot_fan_chart(fit, levels = c(0.5, 0.5)), "duplicates")
  expect_error(plot_fan_chart(fit, levels = c(0.5, 1.2)), "strictly between")
  expect_error(plot_forecast(fit, history = 0), "positive integer")
  expect_error(forecast_error_variance(lm(dist ~ speed, cars)), "fitted ARIMA")
})

test_that("forecast plots build, with dates for quarterly ts", {
  y <- ts(as.numeric(arima.sim(list(ar = 0.8), n = 120)), start = c(1990, 1),
          frequency = 4)
  fit <- forecast::Arima(y, order = c(1, 0, 0))
  p1 <- plot_forecast(fit, steps = 8)
  p2 <- plot_fan_chart(fit, steps = 8)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_silent(ggplot2::ggplot_build(p2))
  expect_true(inherits(attr(dynamic_forecast(fit, 8), "time"), "Date"))
})
