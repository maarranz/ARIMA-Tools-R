# Numerical forecasts and innovation-based forecast uncertainty --------------

#' Forecast-error variances from the psi weights
#'
#' Expands the full ARIMA transfer function, including ordinary and seasonal
#' differencing, to obtain psi_0, ..., psi_(steps-1), and returns
#' `sigma2 * cumsum(psi^2)` for horizons 1 through `steps`.
#'
#' Parameters are treated as fixed and past innovations as known, so
#' parameter-estimation uncertainty is not included. The Kalman-filter
#' variances behind `predict()` / [forecast::forecast()] are exact in finite
#' samples and can differ slightly (e.g. at the start of short samples or with
#' missing values); for typical samples the two coincide to many decimals.
#'
#' The innovation variance is `fit$sigma2`. Note that [forecast::Arima()]
#' reports a degrees-of-freedom-adjusted variance, SSR / (n - number of
#' coefficients), whereas [stats::arima()] and Statsmodels report the ML
#' estimate SSR / n. Using `fit$sigma2` keeps these intervals consistent with
#' `forecast::forecast(fit)` for the same object.
#'
#' @inheritParams root_diagnostics
#' @param steps Number of horizons (positive integer).
#' @return A numeric vector of length `steps`.
#' @export
forecast_error_variance <- function(fit, steps = 20) {
  .check_fit(fit)
  steps <- .check_count(steps, "steps", min = 1)
  sigma2 <- fit$sigma2
  if (!is.finite(sigma2) || sigma2 < 0) {
    stop("The estimated innovation variance must be finite and non-negative.",
         call. = FALSE)
  }
  sigma2 * cumsum(.psi_weights(fit, steps)^2)
}

#' Recursive (dynamic) forecasts with normal-theory intervals
#'
#' Point forecasts from [forecast::forecast()] (which handles means and drift
#' automatically); uncertainty from [forecast_error_variance()], not from the
#' Kalman filter.
#'
#' Models with external regressors other than drift cannot be forecast through
#' this interface.
#'
#' @inheritParams forecast_error_variance
#' @param level Coverage of the central prediction interval.
#' @return A data frame with columns `horizon`, `forecast`, `fev`, `se`,
#'   `lower` and `upper`.
#' @export
dynamic_forecast <- function(fit, steps = 20, level = 0.95) {
  .check_fit(fit)
  .check_level(level, "level")
  fev <- forecast_error_variance(fit, steps)
  fc <- .point_forecast(fit, length(fev))
  se <- sqrt(fev)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(horizon = seq_along(fev), forecast = as.numeric(fc),
                    fev = fev, se = se,
                    lower = as.numeric(fc) - crit * se,
                    upper = as.numeric(fc) + crit * se)
  attr(out, "level") <- level
  attr(out, "time") <- .time_values(fc)
  out
}

.point_forecast <- function(fit, steps) {
  xreg_names <- setdiff(colnames(fit$xreg), "drift")
  if (length(xreg_names)) {
    stop("Models with external regressors (", paste(xreg_names, collapse = ", "),
         ") need future regressor values; use forecast::forecast(fit, xreg = ...).",
         call. = FALSE)
  }
  forecast::forecast(fit, h = steps, level = 95)$mean
}

# Observed history and future time values on a common axis.
.forecast_axis <- function(fit, steps, history) {
  if (!is.null(history)) history <- .check_count(history, "history", min = 1)
  y <- .response(fit)
  obs <- as.numeric(y)
  tt <- .time_values(y)
  fc <- .point_forecast(fit, steps)
  ft <- .time_values(fc)
  if (inherits(tt, "Date") != inherits(ft, "Date") || is.numeric(tt) && !stats::is.ts(y)) {
    tt <- seq_along(obs)
    ft <- length(obs) + seq_len(steps)
  }
  start <- if (is.null(history)) 1L else max(1L, length(obs) - history + 1L)
  keep <- start:length(obs)
  axis_title <- if (inherits(tt, "Date") ||
                     (stats::is.ts(y) && stats::frequency(y) != 1)) "Time" else "Observation"
  list(t = tt[keep], y = obs[keep], future = ft, axis_title = axis_title)
}

.forecast_base <- function(fit, ax, title, subtitle = NULL) {
  hist_df <- data.frame(t = ax$t, y = ax$y, series = "Observed")
  t0 <- utils::tail(ax$t, 1)
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = t0, linetype = "dotted", colour = "grey55") +
    ggplot2::annotate("text", x = t0, y = Inf, label = "Forecast origin",
                      vjust = 1.5, hjust = 1.05, size = 3, colour = "grey40") +
    ggplot2::geom_line(data = hist_df,
                       ggplot2::aes(x = .data$t, y = .data$y, colour = .data$series),
                       linewidth = 0.6) +
    ggplot2::labs(title = title, subtitle = subtitle, x = ax$axis_title,
                  y = .series_name(fit), colour = NULL, fill = NULL,
                  linetype = NULL) +
    .at_theme()
}

#' Plot recursive forecasts with prediction-interval lines
#'
#' Observed history, the point forecast (joined to the last observation), and
#' two dashed interval boundaries; a dotted line marks the forecast origin.
#' All numbers come from [dynamic_forecast()].
#'
#' Related: `autoplot(forecast::forecast(fit))` draws Kalman-filter intervals.
#'
#' @inheritParams dynamic_forecast
#' @param history Number of final observations to display (`NULL` = all).
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A ggplot object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_forecast <- function(fit, steps = 20, level = 0.95, history = 50,
                          interactive = FALSE) {
  tab <- dynamic_forecast(fit, steps = steps, level = level)
  ax <- .forecast_axis(fit, nrow(tab), history)
  lab <- sprintf("%g%% prediction interval", 100 * level)
  t0 <- utils::tail(ax$t, 1); y0 <- utils::tail(ax$y, 1)
  fc_df <- data.frame(t = c(t0, ax$future), y = c(y0, tab$forecast),
                      series = "Forecast")
  pi_df <- rbind(
    data.frame(t = ax$future, y = tab$lower, bound = "Lower"),
    data.frame(t = ax$future, y = tab$upper, bound = "Upper")
  )
  pi_df$band <- lab

  p <- .forecast_base(fit, ax, "Dynamic Forecast") +
    ggplot2::geom_line(data = pi_df,
                       ggplot2::aes(x = .data$t, y = .data$y, group = .data$bound,
                                    linetype = .data$band),
                       colour = "#6484A6", linewidth = 0.6) +
    ggplot2::geom_line(data = fc_df,
                       ggplot2::aes(x = .data$t, y = .data$y, colour = .data$series),
                       linewidth = 1) +
    ggplot2::scale_colour_manual(values = c(Observed = "#444444", Forecast = "#1F5A96"),
                                 breaks = c("Observed", "Forecast")) +
    ggplot2::scale_linetype_manual(values = stats::setNames("dashed", lab))

  if (interactive) .as_interactive(p, tooltip = c("x", "y")) else p
}

#' Fan chart of recursive forecasts
#'
#' Nested central prediction-interval bands (widest drawn first, so the centre
#' is darkest), joined to the last observation, with the observed series and
#' the point forecast drawn on top. Each band comes from [dynamic_forecast()].
#'
#' Related: the 'fanplot' package draws fan charts in base graphics from
#' simulated paths.
#'
#' @inheritParams plot_forecast
#' @param levels Distinct coverages in (0, 1); order does not matter.
#' @return A ggplot object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_fan_chart <- function(fit, steps = 20,
                           levels = c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99),
                           history = 50, interactive = FALSE) {
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) ||
      any(levels <= 0 | levels >= 1)) {
    stop("Every level must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (anyDuplicated(levels)) stop("`levels` must not contain duplicates.", call. = FALSE)
  levels <- sort(levels, decreasing = TRUE)
  tabs <- lapply(levels, function(l) dynamic_forecast(fit, steps = steps, level = l))
  ax <- .forecast_axis(fit, nrow(tabs[[1]]), history)
  t0 <- utils::tail(ax$t, 1); y0 <- utils::tail(ax$y, 1)
  bx <- c(t0, ax$future)
  cover <- paste0(sprintf("%g%%", 100 * rev(levels)), collapse = ", ")

  p <- .forecast_base(fit, ax, "Forecast Fan Chart",
                      paste("Central prediction intervals:", cover))
  L <- length(levels)
  for (i in seq_len(L)) {
    band <- data.frame(t = bx, lo = c(y0, tabs[[i]]$lower), hi = c(y0, tabs[[i]]$upper),
                       band = "Prediction intervals")
    p <- p + ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi, fill = .data$band),
      alpha = 0.08 + 0.12 * (i - 1) / max(1, L - 1)
    )
  }
  fc_df <- data.frame(t = bx, y = c(y0, tabs[[1]]$forecast), series = "Forecast")
  p <- p +
    ggplot2::geom_line(data = fc_df,
                       ggplot2::aes(x = .data$t, y = .data$y, colour = .data$series),
                       linewidth = 1) +
    ggplot2::scale_colour_manual(values = c(Observed = "#444444", Forecast = "#1F5A96"),
                                 breaks = c("Observed", "Forecast")) +
    ggplot2::scale_fill_manual(values = c("Prediction intervals" = "#1F5A96"))

  if (interactive) .as_interactive(p, tooltip = c("x", "y", "ymin", "ymax")) else p
}
