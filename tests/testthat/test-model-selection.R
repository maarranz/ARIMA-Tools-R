set.seed(1357)
y <- arima.sim(list(ar = c(1.2, -0.5)), n = 400)
f1 <- forecast::Arima(y, order = c(2, 0, 0))
f2 <- forecast::Arima(y, order = c(1, 0, 1))

test_that("information criteria agree with forecast::Arima", {
  cm <- compare_models(f1, f2)
  expect_equal(cm$AIC, c(f1$aic, f2$aic))
  expect_equal(cm$AICc, c(f1$aicc, f2$aicc))
  expect_equal(cm$BIC, c(f1$bic, f2$bic))
  expect_equal(cm$k, c(4L, 4L))
  expect_equal(rownames(cm), c("ARIMA(2,0,0)", "ARIMA(1,0,1)"))
  expect_equal(compare_models(a = f1, b = f2)$model, c("a", "b"))
  expect_error(compare_models(f1), "At least two")
  expect_error(compare_models(f1, f1), "unique")
})

test_that("order search covers the grid and is sorted", {
  expect_message(
    s <- order_search(y, p_max = 2, q_max = 2, d_min = 0, d_max = 1, criterion = "aic"),
    "same d")
  expect_equal(nrow(s), 18L)
  expect_false(is.unsorted(s$d))                     # grouped by d
  for (dd in 0:1) {
    ok <- s[s$converged & s$d == dd, ]
    expect_false(is.unsorted(ok$AIC))                # sorted within d
    best <- select_model(s, "AIC", d = dd)
    expect_equal(best$AIC, min(ok$AIC))
    expect_equal(best$d, dd)
  }
  expect_error(select_model(s, "AIC"), "different orders of differencing")
  expect_error(select_model(s, "AIC", d = 2), "No models with d = 2")
  expect_equal(attr(s, "criterion"), "AIC")
  expect_error(order_search(y, criterion = "foo"), "criterion")
  expect_error(order_search(y, d_min = 2, d_max = 1), "d_min")
})

test_that("auto_arima is nonseasonal, BIC-based and summarized", {
  a <- auto_arima(y)
  expect_s3_class(a, "Arima")
  expect_equal(attr(a, "ic"), "BIC")
  s <- auto_arima_summary(a)
  expect_equal(s[["Criterion value"]], a$bic)
  expect_match(s[["Selected model"]], "^ARIMA\\(")
  set.seed(24680)
  yd <- cumsum(0.4 + rnorm(300))
  sd <- auto_arima_summary(auto_arima(yd))
  expect_equal(sd[["Deterministic component"]], "Drift")
})


test_that("F3: CSS fits have no information criteria and are not selected", {
  set.seed(303)
  y3 <- ts(as.numeric(arima.sim(list(ar = 0.6, ma = 0.2), n = 100)))
  css <- forecast::Arima(y3, order = c(1, 0, 1), method = "CSS")
  ml <- forecast::Arima(y3, order = c(1, 0, 1), method = "ML")
  cm <- compare_models(CSS = css, ML = ml)
  expect_true(is.na(cm["CSS", "AIC"]) && is.na(cm["CSS", "BIC"]))
  expect_equal(select_model(cm, "AIC")$model, "ML")
  s <- order_search(y3, p_max = 1, q_max = 1, method = "CSS")
  expect_true(all(is.na(s$BIC)))
  expect_true(all(grepl("^CSS", s$status)))
  expect_error(select_model(s, "BIC"), "No eligible")
})

test_that("F4: select_model checks d and sample size in compare_models tables", {
  f0 <- forecast::Arima(y, order = c(1, 0, 0))
  f1 <- forecast::Arima(y, order = c(1, 1, 0))
  cm <- compare_models(f0, f1)
  expect_error(select_model(cm, "BIC"), "different orders of differencing")
  expect_equal(select_model(cm, "BIC", d = 1)$model, "ARIMA(1,1,0)")
  fs <- forecast::Arima(window(y, end = 300), order = c(1, 0, 0))
  expect_error(select_model(compare_models(f0, fs, names = c("all", "part")), "BIC"),
               "different numbers of observations")
})
