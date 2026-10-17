# Residual diagnostics, dynamic effects, recursive parameter stability -------

#' Ljung-Box tests on the residuals
#'
#' Ljung-Box statistics at one or several lags, with degrees of freedom
#' adjusted for the number of estimated ARMA coefficients,
#' `df = h - (p + q + P + Q)`. Lags with non-positive degrees of freedom are
#' rejected. Models whose residuals have gaps inside the sample are rejected,
#' because the test assumes consecutive residuals.
#'
#' Related: `stats::Box.test(residuals(fit), lag = h, type = "Ljung-Box",
#' fitdf = p + q)` for one lag; `forecast::checkresiduals(fit)`.
#'
#' @inheritParams root_diagnostics
#' @param lags Integer vector of lags.
#' @return A data frame with columns `lag`, `statistic`, `df` and `p.value`.
#' @export
ljung_box_test <- function(fit, lags) {
  .check_fit(fit)
  if (!is.numeric(lags) || length(lags) == 0L || anyNA(lags) ||
      any(!is.finite(lags)) || any(lags != round(lags)) || any(lags <= 0)) {
    stop("All lags must be positive integers.", call. = FALSE)
  }
  lags <- as.integer(lags)
  o <- .orders(fit)
  model_df <- o$p + o$q + o$P + o$Q
  bad <- lags[lags <= model_df]
  if (length(bad)) {
    stop("Ljung-Box lags must be greater than the number of ARMA ",
         "coefficients (", model_df, "). Invalid lag(s): ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  r <- .resid_consecutive(fit)
  res <- lapply(lags, function(h) {
    bt <- stats::Box.test(r, lag = h, type = "Ljung-Box", fitdf = model_df)
    c(unname(bt$statistic), unname(bt$parameter), bt$p.value)
  })
  res <- do.call(rbind, res)
  data.frame(lag = lags, statistic = res[, 1], df = as.integer(res[, 2]),
             p.value = res[, 3])
}

#' Residual ACF and PACF with white-noise bounds
#'
#' Lag zero is omitted. The bounds are the individual approximate white-noise
#' bounds +/- z / sqrt(T); they are not a joint test (use [ljung_box_test()]).
#' Residuals with gaps inside the sample are rejected, because the lags would
#' no longer correspond to calendar distances.
#'
#' @inheritParams root_diagnostics
#' @param lags Maximum lag.
#' @param confidence Confidence level for the bounds.
#' @return A data frame with columns `lag`, `acf`, `pacf`, `lower` and `upper`.
#' @export
residual_correlations <- function(fit, lags = 20, confidence = 0.95) {
  .check_fit(fit)
  lags <- .check_count(lags, "lags", min = 1)
  .check_level(confidence)
  r <- .resid_consecutive(fit)
  n <- length(r)
  if (lags >= n) {
    stop("`lags` must be smaller than the number of residuals (", n, ").",
         call. = FALSE)
  }
  a <- stats::acf(r, lag.max = lags, plot = FALSE, demean = TRUE)$acf[-1]
  pa <- stats::pacf(r, lag.max = lags, plot = FALSE)$acf[, 1, 1]
  bound <- stats::qnorm(1 - (1 - confidence) / 2) / sqrt(n)
  data.frame(lag = seq_len(lags), acf = a, pacf = pa,
             lower = -bound, upper = bound)
}

# One correlogram panel (stems, markers, shaded white-noise band).
.correlogram_panel <- function(lag, value, bound, ylab, title, lags) {
  df <- data.frame(lag = lag, value = value,
                   text = sprintf("Lag: %d<br>%s: %.4f", lag, ylab, value))
  ymax <- min(1, max(0.25, 1.15 * max(abs(value), bound)))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$lag, y = .data$value)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -bound,
                      ymax = bound, alpha = 0.12, fill = "#1F5A96") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                          linewidth = 0.8, colour = "#1F5A96") +
    .geom_point_text(colour = "#1F5A96", size = 2.2) +
    ggplot2::scale_x_continuous(breaks = .lag_breaks(lags),
                                limits = c(0.5, lags + 0.5)) +
    ggplot2::coord_cartesian(ylim = c(-ymax, ymax)) +
    ggplot2::labs(title = title, x = "Lag", y = "Correlation") +
    .at_theme()
}

.lag_breaks <- function(lags) {
  step <- if (lags <= 24) 1 else if (lags <= 60) 5 else 10
  unique(c(1, seq(step, lags, by = step)))
}

# geom_point carrying a `text` aesthetic for plotly tooltips. ggplot2 warns
# about the unknown aesthetic; the warning is harmless and silenced here.
.geom_point_text <- function(...) {
  suppressWarnings(ggplot2::geom_point(ggplot2::aes(text = .data$text), ...))
}

#' Residual diagnostic plot
#'
#' Four stacked panels: residuals over time, residual ACF, residual PACF, and
#' Ljung-Box p-values by lag (only lags with positive degrees of freedom), with
#' a reference line at 0.05.
#'
#' Related: `forecast::checkresiduals(fit)` (residuals, ACF, histogram, one
#' Ljung-Box test) and `stats::tsdiag(fit)` (whose Ljung-Box p-values are
#' *not* adjusted for the estimated parameters).
#'
#' @inheritParams residual_correlations
#' @param title Overall title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A patchwork object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_residual_diagnostics <- function(fit, lags = 20, confidence = 0.95,
                                      title = "Residual Diagnostics",
                                      interactive = FALSE) {
  cors <- residual_correlations(fit, lags = lags, confidence = confidence)
  bound <- cors$upper[1]

  res_ts <- stats::residuals(fit)
  tt <- .time_values(res_ts)
  rdf <- data.frame(t = tt, r = as.numeric(res_ts))
  rdf <- rdf[is.finite(rdf$r), ]
  rdf$text <- sprintf("Observation: %s<br>Residual: %.4f", format(rdf$t), rdf$r)

  o <- .orders(fit)
  first <- o$p + o$q + o$P + o$Q + 1
  lb <- if (first <= lags) ljung_box_test(fit, first:lags) else
    data.frame(lag = integer(0), p.value = numeric(0))
  lb$text <- sprintf("Lag: %d<br>p-value: %.4f", lb$lag, lb$p.value)

  p1 <- ggplot2::ggplot(rdf, ggplot2::aes(x = .data$t, y = .data$r)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(colour = "#1F5A96", linewidth = 0.5) +
    suppressWarnings(ggplot2::geom_point(ggplot2::aes(text = .data$text),
                                         alpha = 0, size = 0.1)) +
    ggplot2::labs(title = "Residuals", x = "Observation", y = "Residual") +
    .at_theme()

  p2 <- .correlogram_panel(cors$lag, cors$acf, bound, "ACF",
                           "Autocorrelation Function", lags)
  p3 <- .correlogram_panel(cors$lag, cors$pacf, bound, "PACF",
                           "Partial Autocorrelation Function", lags)

  p4 <- ggplot2::ggplot(lb, ggplot2::aes(x = .data$lag, y = .data$p.value)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
    ggplot2::annotate("text", x = lags + 0.5, y = 0.05, label = "0.05",
                      hjust = 1, vjust = -0.5, size = 3) +
    ggplot2::geom_line(colour = "#1F5A96") +
    .geom_point_text(colour = "#1F5A96", size = 2) +
    ggplot2::scale_x_continuous(breaks = .lag_breaks(lags),
                                limits = c(0.5, lags + 0.5)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Ljung-Box Test", x = "Lag", y = "p-value") +
    .at_theme()

  .stack_panels(list(p1, p2, p3, p4), heights = c(0.30, 0.24, 0.24, 0.22),
                title = title, interactive = interactive)
}

#' Jarque-Bera normality test on the residuals
#'
#' The statistic uses the (biased) sample skewness and kurtosis, as in
#' `tseries::jarque.bera.test()` and SciPy. The reported skewness and kurtosis
#' are the bias-corrected versions; kurtosis is on the scale where the Normal
#' distribution has kurtosis 3.
#'
#' @inheritParams root_diagnostics
#' @return A named list of class `arimatools_summary`.
#' @export
jarque_bera_test <- function(fit) {
  .check_fit(fit)
  x <- .resid(fit)
  n <- length(x)
  if (n < 3) stop("At least three finite residuals are required.", call. = FALSE)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m3 <- mean((x - m)^3)
  m4 <- mean((x - m)^4)
  s <- m3 / m2^1.5
  k <- m4 / m2^2
  jb <- n / 6 * (s^2 + (k - 3)^2 / 4)
  skew_adj <- sqrt(n * (n - 1)) / (n - 2) * s
  kurt_adj <- if (n > 3) {
    ((n^2 - 1) * k - 3 * (n - 1)^2) / ((n - 2) * (n - 3)) + 3
  } else NA_real_
  .at_summary(
    list("Jarque-Bera statistic" = jb, "Degrees of freedom" = 2L,
         "P-value" = stats::pchisq(jb, 2, lower.tail = FALSE),
         "Skewness" = skew_adj, "Kurtosis" = kurt_adj),
    "Jarque-Bera test"
  )
}

#' Dynamic effects of a one-unit innovation
#'
#' Psi weights of the fitted model for horizons 0 through `horizon`, and their
#' cumulative sum. For integrated models (ordinary or seasonal), the weights
#' incorporate the differencing polynomial, so they are responses of the
#' original (level) series and need not decay to zero.
#'
#' Models with a Box-Cox transformation (`lambda`) are rejected: their psi
#' weights would describe the transformed series, not the original one.
#'
#' Related: [stats::ARMAtoMA()] for the stationary ARMA part only.
#'
#' @inheritParams root_diagnostics
#' @param horizon Maximum horizon (non-negative integer).
#' @return A data frame with columns `horizon`, `effect` and `cumulative`.
#' @export
dynamic_effects <- function(fit, horizon = 20) {
  .check_fit(fit)
  horizon <- .check_count(horizon, "horizon", min = 0)
  .reject_lambda(fit, "dynamic_effects()")
  psi <- .psi_weights(fit, horizon + 1L)
  data.frame(horizon = 0:horizon, effect = psi, cumulative = cumsum(psi))
}

#' Plot dynamic and cumulative effects
#'
#' Upper panel: psi weights as stems. Lower panel: their cumulative sum.
#'
#' @inheritParams dynamic_effects
#' @param title Overall title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A patchwork object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_dynamic_effects <- function(fit, horizon = 20,
                                 title = "Dynamic Effects of an Innovation",
                                 interactive = FALSE) {
  tab <- dynamic_effects(fit, horizon = horizon)
  tab$text1 <- sprintf("Horizon: %d<br>Dynamic effect: %.5f", tab$horizon, tab$effect)
  tab$text2 <- sprintf("Horizon: %d<br>Cumulative effect: %.5f", tab$horizon, tab$cumulative)
  brks <- if (horizon <= 30) 0:horizon else pretty(0:horizon)

  p1 <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$horizon, y = .data$effect)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$horizon, yend = 0),
                          colour = "#1F5A96", linewidth = 0.8) +
    suppressWarnings(ggplot2::geom_point(ggplot2::aes(text = .data$text1),
                                         colour = "#1F5A96", size = 2.2)) +
    ggplot2::scale_x_continuous(breaks = brks) +
    ggplot2::labs(title = "Dynamic Effects", x = NULL, y = "Effect") +
    .at_theme()

  p2 <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$horizon, y = .data$cumulative)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(colour = "#1F5A96") +
    suppressWarnings(ggplot2::geom_point(ggplot2::aes(text = .data$text2),
                                         colour = "#1F5A96", size = 1.8)) +
    ggplot2::scale_x_continuous(breaks = brks) +
    ggplot2::labs(title = "Cumulative Effects", x = "Horizon",
                  y = "Cumulative effect") +
    .at_theme()

  .stack_panels(list(p1, p2), heights = c(0.5, 0.5), title = title,
                interactive = interactive)
}

#' Recursive (expanding-window) ARIMA estimation
#'
#' Re-estimates the model on samples `1:initial`, `1:(initial + 1)`, ..., `1:n`
#' with [forecast::Arima()] and stores the estimates, standard errors and
#' normal confidence intervals of every coefficient, plus the innovation
#' variance (`sigma2`, without standard error).
#'
#' If `initial` is `NULL`, the initial window is the larger of 10% of the
#' sample and the window that leaves at least `min_df` residual degrees of
#' freedom. Routine estimation warnings are silenced; convergence status is
#' kept in the output. If the default CSS-ML estimation fails at an endpoint,
#' exact ML is tried before recording the endpoint as failed.
#'
#' @param y A numeric vector or `ts` object. Leading and trailing missing
#'   values of a `ts` are dropped.
#' @param order Integer vector `c(p, d, q)`.
#' @param initial `NULL`, a fraction in (0, 1), or an integer number of
#'   observations for the first window.
#' @param min_df Minimum residual degrees of freedom for the automatic rule.
#' @param confidence Confidence level for the intervals.
#' @param include.constant Passed to [forecast::Arima()]. `NULL` (default)
#'   keeps forecast's default: a mean when d = 0, no constant when d > 0.
#'   `TRUE` gives a mean when d = 0 and a drift when d = 1; `FALSE` none.
#' @param ... Further arguments for [forecast::Arima()] (e.g. `method`).
#'   External regressors (`xreg`) and Box-Cox transformations (`lambda`) are
#'   not supported.
#' @return A long data frame with columns `endpoint`, `n`, `term`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high` and `converged`. `endpoint` is a
#'   `Date` for monthly, quarterly and annual `ts` objects, decimal time for
#'   other `ts`, and the observation number otherwise. Endpoints where
#'   estimation failed have one row with missing values; the error messages are
#'   kept in `attr(result, "failures")`.
#' @export
recursive_estimation <- function(y, order, initial = NULL, min_df = 20,
                                 confidence = 0.95, include.constant = NULL,
                                 ...) {
  dot_names <- .resolve_dot_names(list(...), forecast::Arima)
  if (any(c("xreg", "lambda", "y", "order") %in% dot_names)) {
    stop("recursive_estimation() supports univariate models without Box-Cox ",
         "transformation: `xreg` and `lambda` cannot be passed through `...` ",
         "(and `y`, `order` are set by the function).", call. = FALSE)
  }
  if (!is.numeric(order) || length(order) != 3L ||
      any(!vapply(order, .is_count, logical(1)))) {
    stop("`order` must be c(p, d, q) with non-negative integers.", call. = FALSE)
  }
  order <- as.integer(order)
  p <- order[1]; d <- order[2]; q <- order[3]
  min_df <- .check_count(min_df, "min_df", min = 1)
  .check_level(confidence)

  if (stats::is.ts(y)) {
    y <- stats::na.omit(y)
  } else {
    if (!is.null(dim(y)) && NCOL(y) > 1) stop("`y` must be univariate.", call. = FALSE)
    y <- as.numeric(y)
    y <- y[!is.na(y)]
  }
  if (!all(is.finite(y))) stop("`y` must contain only finite values.", call. = FALSE)
  n <- length(y)
  times <- .time_values(y)

  k <- p + q + 1L
  has_const <- if (is.null(include.constant)) d == 0 else isTRUE(include.constant) && d <= 1
  k <- k + as.integer(has_const)
  df_window <- k + min_df + d

  if (is.null(initial)) {
    initial_n <- max(floor(0.10 * n), df_window)
  } else if (is.numeric(initial) && length(initial) == 1L && initial > 0 && initial < 1) {
    initial_n <- floor(initial * n)
  } else if (.is_count(initial, min = 1)) {
    initial_n <- as.integer(initial)
  } else {
    stop("`initial` must be NULL, a fraction in (0, 1), or a positive integer.",
         call. = FALSE)
  }
  if (initial_n >= n) {
    stop("The initial window is too large: no recursive estimation period remains.",
         call. = FALSE)
  }
  approx_df <- initial_n - d - k
  if (approx_df <= 0) {
    stop("The initial window is too small for ARIMA(", paste(order, collapse = ","),
         "). Approximate residual degrees of freedom: ", approx_df, ".",
         call. = FALSE)
  }

  crit <- stats::qnorm(1 - (1 - confidence) / 2)
  is_ts <- stats::is.ts(y)

  rows <- lapply(seq.int(initial_n, n), function(end) {
    ys <- if (is_ts) {
      stats::ts(as.numeric(y)[seq_len(end)], start = stats::start(y),
                frequency = stats::frequency(y))
    } else y[seq_len(end)]
    fit <- tryCatch(.fit_arima(ys, order, include.constant, ...),
                    error = function(e) conditionMessage(e))
    if (is.character(fit)) {
      return(structure(
        data.frame(endpoint = times[end], n = end, term = NA_character_,
                   estimate = NA_real_, std.error = NA_real_,
                   conf.low = NA_real_, conf.high = NA_real_,
                   converged = FALSE, stringsAsFactors = FALSE),
        failure = fit))
    }
    est <- .free_coef(fit)
    se <- if (length(est)) sqrt(diag(fit$var.coef)) else numeric(0)
    data.frame(
      endpoint = times[end], n = end,
      term = c(names(est), "sigma2"),
      estimate = c(unname(est), fit$sigma2),
      std.error = c(unname(se), NA_real_),
      conf.low = c(unname(est - crit * se), NA_real_),
      conf.high = c(unname(est + crit * se), NA_real_),
      converged = .converged(fit),
      stringsAsFactors = FALSE
    )
  })
  failed <- vapply(rows, function(r) !is.null(attr(r, "failure")), logical(1))
  failures <- if (any(failed)) {
    data.frame(n = vapply(rows[failed], function(r) r$n[1], numeric(1)),
               message = vapply(rows[failed], function(r) attr(r, "failure"),
                                character(1)),
               stringsAsFactors = FALSE)
  } else data.frame(n = numeric(0), message = character(0))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "failures") <- failures
  attr(out, "order") <- order
  attr(out, "confidence") <- confidence
  out
}

#' Plot recursive parameter estimates
#'
#' One panel per parameter with the recursive estimate, its confidence band,
#' an optional dashed line at the full-sample estimate, and an optional
#' vertical line at a known or hypothesized break point. The break point is a
#' graphical reference only; it is not a formal test.
#'
#' @param recursive Output of [recursive_estimation()].
#' @param parameters Character vector of terms to plot. `NULL` plots every
#'   coefficient except `sigma2`.
#' @param reference Show the full-sample estimate?
#' @param breakpoint Optional value on the `endpoint` scale (e.g. a `Date`, or
#'   an observation number for plain vectors).
#' @param title Plot title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A ggplot object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_parameter_stability <- function(recursive, parameters = NULL,
                                     reference = TRUE, breakpoint = NULL,
                                     title = "Recursive Parameter Stability",
                                     interactive = FALSE) {
  need <- c("endpoint", "n", "term", "estimate", "conf.low", "conf.high", "converged")
  miss <- setdiff(need, names(recursive))
  if (length(miss)) {
    stop("`recursive` is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  available <- unique(stats::na.omit(recursive$term))
  if (!length(available)) stop("No parameter estimates are available.", call. = FALSE)
  selected <- if (is.null(parameters)) setdiff(available, "sigma2") else unique(parameters)
  if (!length(selected)) stop("No parameters have been selected.", call. = FALSE)
  unknown <- setdiff(selected, available)
  if (length(unknown)) {
    stop("Unknown parameter(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "), call. = FALSE)
  }

  dat <- recursive[!is.na(recursive$term) & recursive$term %in% selected, ]
  dat <- dat[order(dat$term, dat$n), ]
  dat$term <- factor(dat$term, levels = selected)
  dat$text <- sprintf("Endpoint: %s<br>Estimate: %.5f", format(dat$endpoint),
                      dat$estimate)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$endpoint)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
                         fill = "grey50", alpha = 0.25, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate, colour = "Recursive estimate"),
                       linewidth = 0.8, na.rm = TRUE) +
    suppressWarnings(ggplot2::geom_point(
      ggplot2::aes(y = .data$estimate, text = .data$text),
      alpha = 0, size = 0.1, na.rm = TRUE))

  if (reference) {
    ok <- dat[!is.na(dat$estimate), ]
    ref <- do.call(rbind, lapply(split(ok, ok$term, drop = TRUE), function(d) d[which.max(d$n), ]))
    ref$label <- sprintf("Full sample: %.4f", ref$estimate)
    p <- p +
      ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = .data$estimate),
                          linetype = "dashed", colour = "grey20") +
      ggplot2::geom_text(data = ref,
                         ggplot2::aes(x = .data$endpoint, y = .data$estimate,
                                      label = .data$label),
                         hjust = 1, vjust = -0.6, size = 3)
  }
  if (!is.null(breakpoint)) {
    p <- p + ggplot2::geom_vline(xintercept = breakpoint, linetype = "dashed",
                                 colour = "#B2182B")
  }

  p <- p +
    ggplot2::facet_wrap(~term, ncol = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c("Recursive estimate" = "#1F5A96"),
                                 name = NULL) +
    ggplot2::labs(title = title, x = "Recursive sample endpoint", y = "Estimate") +
    .at_theme()

  if (interactive) .as_interactive(p) else p
}
