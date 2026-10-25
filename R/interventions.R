# Interventions: future values, fitted effects and plots ----------------------
#
# Works with the regressor matrices created by outlier_regressors(), which store
# the definition of every event: attr "outliers" (type, date, index), attr
# "delta" (TC decay), and, from version 0.5.0, attr "row_start" (0 for the
# original sample, n for a block that continues it) and attr "origin_length".
# Matches extend_interventions(), intervention_effects() and plot_interventions()
# in the Python package ARIMA-Tools.

# Plain numeric matrix (drops ts and other attributes) for comparisons.
.plain <- function(X) matrix(as.numeric(X), nrow = NROW(X))

# Values of the intervention regressors at overall observations
# start + 1, ..., start + len (1-based positions counted from the start of the
# original sample).
.intervention_values <- function(types, idx, delta, start, len) {
  tt <- start + seq_len(len)
  X <- vapply(seq_along(idx), function(j) {
    age <- tt - idx[j]
    switch(types[j],
           AO = as.numeric(age == 0),
           LS = as.numeric(age >= 0),
           TC = ifelse(age >= 0, delta^pmax(age, 0), 0))
  }, numeric(len))
  matrix(X, nrow = len)
}

# Validate a matrix from outlier_regressors() (or extend_interventions()) and
# return its definition. Rejects matrices whose values no longer match the
# stored events (edited, subset by rows, or combined with other columns).
.intervention_meta <- function(X, arg = "X") {
  ev <- attr(X, "outliers")
  if (is.null(dim(X)) || !is.data.frame(ev) ||
      !all(c("type", "date", "index") %in% names(ev)) || is.null(attr(X, "delta"))) {
    stop("`", arg, "` must be a matrix created by outlier_regressors() (or ",
         "found$regressors from detect_outliers()), which stores the definition of ",
         "each intervention. Keep that matrix: subsetting or combining it with other ",
         "columns removes the definitions.", call. = FALSE)
  }
  types <- as.character(ev$type); idx <- as.integer(ev$index)
  delta <- attr(X, "delta")
  start <- attr(X, "row_start"); if (is.null(start)) start <- 0L
  origin <- attr(X, "origin_length"); if (is.null(origin)) origin <- nrow(X)
  nm <- paste0(types, ev$date)
  if (ncol(X) != length(types) || !identical(colnames(X), nm)) {
    stop("The columns of `", arg, "` do not match its stored interventions (",
         paste(nm, collapse = ", "), ").", call. = FALSE)
  }
  ref <- .intervention_values(types, idx, delta, start, nrow(X))
  if (!isTRUE(all.equal(.plain(X), ref, tolerance = 1e-12))) {
    stop("The values of `", arg, "` no longer match its stored interventions. ",
         "Use the matrix exactly as returned by outlier_regressors().", call. = FALSE)
  }
  list(types = types, idx = idx, dates = as.character(ev$date), delta = delta,
       start = as.integer(start), origin = as.integer(origin), names = nm)
}

#' Future values of intervention regressors
#'
#' Continues the AO, LS and TC regressors created by [outlier_regressors()]
#' beyond the end of the sample, for forecasting a model that contains them
#' (see the `xreg` argument of [dynamic_forecast()]).
#'
#' An additive outlier happened once, so its future values are 0. A level shift
#' is permanent, so its future values are 1. A temporary change keeps decaying
#' from its original date: at observation t its value is `delta^(t - t0)`, so
#' the decay does not restart at the forecast origin.
#'
#' The function continues the interventions you specified; it does not add new
#' future events. A model with other regressors needs their future values too,
#' supplied by you as additional columns.
#'
#' @param X Matrix of intervention regressors from [outlier_regressors()], or
#'   `found$regressors` from [detect_outliers()], unchanged. A block returned by
#'   `extend_interventions()` can itself be extended (the event ages continue).
#' @param steps Number of future observations (positive integer).
#' @return A matrix with `steps` rows and the same columns as `X`: a `ts`
#'   matrix continuing the time axis when `X` is a `ts`, an ordinary matrix
#'   otherwise. It keeps the definitions of the interventions.
#' @examples
#' y <- ts(rnorm(40), start = c(2015, 1), frequency = 4)
#' X <- outlier_regressors(y, dates = c("2020Q2", "2020Q3", "2021Q1"),
#'                         types = c("AO", "TC", "LS"), delta = 0.5)
#' extend_interventions(X, steps = 4)   # AO stays 0, LS stays 1, TC keeps decaying
#'
#' # Forecasting a model with an intervention
#' set.seed(2)
#' z <- ts(cumsum(rnorm(100)), start = c(2000, 1), frequency = 4)
#' z[60:100] <- z[60:100] + 5                       # level shift in 2014Q4
#' L <- outlier_regressors(z, dates = "2014Q4", types = "LS")
#' fit <- forecast::Arima(z, order = c(1, 1, 0), xreg = L)
#' dynamic_forecast(fit, steps = 4, xreg = extend_interventions(L, 4))
#' @export
extend_interventions <- function(X, steps = 20) {
  m <- .intervention_meta(X)
  steps <- .check_count(steps, "steps", min = 1)
  start <- m$start + nrow(X)
  out <- .intervention_values(m$types, m$idx, m$delta, start, steps)
  colnames(out) <- m$names
  if (stats::is.ts(X)) {
    f <- stats::frequency(X)
    out <- stats::ts(out, start = stats::tsp(X)[2] + 1 / f, frequency = f)
  }
  attr(out, "outliers") <- attr(X, "outliers")
  attr(out, "delta") <- m$delta
  attr(out, "row_start") <- start
  attr(out, "origin_length") <- m$origin
  out
}

#' Estimated intervention effects
#'
#' The contribution of each intervention to the observed series, its estimated
#' coefficient times its regressor, together with their sum and the adjusted
#' series (observed minus the combined effect).
#'
#' The adjusted series is not a residual: the ARIMA dynamics, the mean or drift
#' and any other regressors are still in it. The contributions are on the
#' scale of the modelled series (for a model of `log(y)`, on the log scale).
#' They treat the estimated coefficients, dates, types and `delta` as known.
#'
#' A warning is issued if the model did not converge or if the standard errors
#' of the intervention coefficients are not available, which suggests that the
#' interventions are not identified (for example, an LS on the first
#' observation together with a mean, or an LS and a TC on the last
#' observation).
#'
#' @param fit A model fitted with [forecast::Arima()] whose `xreg` contains the
#'   columns of `X`, or the output of [detect_outliers()] (then `X` is taken
#'   from it).
#' @param X The matrix from [outlier_regressors()] used in estimation,
#'   unchanged and covering the original sample.
#' @return A data frame with columns `time`, `observed`, one column per
#'   intervention (named as in `X`), `combined` and `adjusted`. Attribute
#'   `coefficients` holds the estimated intervention coefficients.
#' @examples
#' set.seed(440)
#' y <- ts(arima.sim(list(ar = 0.65), n = 120), start = c(1990, 1), frequency = 4)
#' y[60] <- y[60] + 6                   # additive outlier in 2004Q4
#' y[90:120] <- y[90:120] + 3           # level shift in 2012Q2
#' X <- outlier_regressors(y, dates = c("2004Q4", "2012Q2"), types = c("AO", "LS"))
#' fit <- forecast::Arima(y, order = c(1, 0, 0), xreg = X)
#' eff <- intervention_effects(fit, X)
#' attr(eff, "coefficients")
#' eff[58:62, ]
#' @export
intervention_effects <- function(fit, X = NULL) {
  if (inherits(fit, "arimatools_outliers")) {
    if (is.null(X)) X <- fit$regressors
    if (is.null(X)) stop("No outliers were detected, so there are no effects.",
                         call. = FALSE)
    fit <- fit$fit
  }
  .check_fit(fit)
  if (is.null(X)) stop("Supply `X`, the matrix from outlier_regressors() used in ",
                       "estimation.", call. = FALSE)
  .reject_lambda(fit, "intervention_effects()")
  m <- .intervention_meta(X)
  if (m$start != 0L) {
    stop("`X` must cover the original sample, not future observations.", call. = FALSE)
  }
  y <- .response(fit)
  if (length(y) != nrow(X)) {
    stop("`X` has ", nrow(X), " rows but the model was fitted to ", length(y),
         " observations.", call. = FALSE)
  }
  xr <- fit$xreg
  missing <- setdiff(m$names, colnames(xr))
  if (length(missing)) {
    stop("The model was not fitted with these interventions: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!isTRUE(all.equal(.plain(xr[, m$names, drop = FALSE]), .plain(X)))) {
    stop("The interventions in `X` differ from the regressors used in estimation.",
         call. = FALSE)
  }
  beta <- fit$coef[m$names]
  if (any(!is.finite(beta)) || any(!.free(fit)[match(m$names, names(fit$coef))])) {
    stop("Some intervention coefficients are fixed or not finite.", call. = FALSE)
  }
  if (!.converged(fit)) {
    warning("The model did not converge (optimizer code ", fit$code, "); the effects ",
            "use its reported coefficients.", call. = FALSE)
  }
  vc <- fit$var.coef
  se <- if (!is.null(vc) && all(m$names %in% rownames(vc))) {
    suppressWarnings(sqrt(diag(vc)[m$names]))
  } else rep(NA_real_, length(m$names))
  if (any(!is.finite(se))) {
    warning("Standard errors are not available for: ",
            paste(m$names[!is.finite(se)], collapse = ", "),
            ". The interventions may not be identified (e.g. collinear with the ",
            "mean, the drift or each other).", call. = FALSE)
  }
  contrib <- sweep(.plain(X), 2, beta, `*`)
  colnames(contrib) <- m$names
  comb <- rowSums(contrib)
  out <- data.frame(time = .time_values(y), observed = as.numeric(y),
                    contrib, combined = comb, adjusted = as.numeric(y) - comb,
                    check.names = FALSE)
  attr(out, "coefficients") <- beta
  out
}

#' Plot interventions: chosen dates, estimated effects
#'
#' Without a fitted model, shows the series with the dates of the
#' interventions you specified, marked by type: a check of the specification
#' before estimation. With a fitted model, adds the estimated effects: the
#' observed and adjusted series in the upper panel and the combined
#' intervention effect (optionally each intervention separately) in the lower
#' panel. Several interventions on the same date share one marker.
#'
#' For outliers found by [detect_outliers()], or to plot any fitted model
#' directly, see [plot_outliers()]; numbers are given by
#' [intervention_effects()].
#'
#' @param y The observed series (the one used in estimation).
#' @param X The matrix from [outlier_regressors()] for `y`.
#' @param fit Optional fitted model containing the columns of `X` as `xreg`.
#' @param individual With `fit`: also draw the effect of each intervention.
#' @param title Overall title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A ggplot object without `fit`, a patchwork object with `fit`, or a
#'   plotly widget when `interactive = TRUE`.
#' @examples
#' set.seed(440)
#' y <- ts(arima.sim(list(ar = 0.65), n = 120), start = c(1990, 1), frequency = 4)
#' y[60] <- y[60] + 6
#' y[90:120] <- y[90:120] + 3
#' X <- outlier_regressors(y, dates = c("2004Q4", "2012Q2"), types = c("AO", "LS"))
#' plot_interventions(y, X)                       # the specification
#' fit <- forecast::Arima(y, order = c(1, 0, 0), xreg = X)
#' plot_interventions(y, X, fit, individual = TRUE)   # the estimated effects
#' @export
plot_interventions <- function(y, X, fit = NULL, individual = FALSE,
                               title = "Intervention Effects", interactive = FALSE) {
  if (!is.logical(individual) || length(individual) != 1L || is.na(individual)) {
    stop("`individual` must be TRUE or FALSE.", call. = FALSE)
  }
  m <- .intervention_meta(X)
  if (m$start != 0L) stop("`X` must cover the original sample.", call. = FALSE)
  if (!is.null(dim(y)) && NCOL(y) > 1) stop("`y` must be univariate.", call. = FALSE)
  if (length(y) != nrow(X)) {
    stop("`y` and `X` must have the same number of observations.", call. = FALSE)
  }
  if (stats::is.ts(y) && stats::is.ts(X) &&
      !isTRUE(all.equal(stats::tsp(y), stats::tsp(X)))) {
    stop("`y` and `X` cover different periods.", call. = FALSE)
  }
  ylab <- deparse1(substitute(y))
  if (is.null(fit)) {
    if (individual) warning("`individual` needs a fitted model; ignored.", call. = FALSE)
    return(.effects_figure(y, X, beta = NULL, word = "Intervention", ylab = ylab,
                           title = title, interactive = interactive))
  }
  tab <- intervention_effects(fit, X)
  if (!isTRUE(all.equal(as.numeric(y), tab$observed))) {
    stop("`y` differs from the series used to fit the model.", call. = FALSE)
  }
  .effects_figure(y, X, beta = attr(tab, "coefficients"), word = "Intervention",
                  ylab = ylab, title = title, interactive = interactive,
                  individual = individual)
}

# Shared figure for plot_outliers() and plot_interventions(). With beta = NULL,
# one panel with the series and the event dates; otherwise observed and
# adjusted series plus the combined (and optionally individual) effects.
.effects_figure <- function(y, X, beta, word, ylab, title, interactive,
                            individual = FALSE, none_label = NULL) {
  n <- length(y)
  tt <- .time_values(y)
  have <- !is.null(X) && ncol(X) > 0
  if (have) {
    Xm <- .plain(X); colnames(Xm) <- colnames(X)
    first <- apply(Xm, 2, function(col) which(col != 0)[1])
    one <- if (is.null(beta)) colnames(Xm) else sprintf("%s: %.4f", colnames(Xm), beta)
    pos <- sort(unique(first))
    grp <- lapply(pos, function(i) which(first == i))
    eff_all <- if (is.null(beta)) numeric(n) else as.numeric(Xm %*% beta)
    pts <- data.frame(
      t = tt[pos], y = as.numeric(y)[pos], eff = eff_all[pos],
      type = vapply(grp, function(g) {
        tp <- unique(substr(colnames(Xm)[g], 1, 2))
        if (length(tp) > 1) "Several" else tp
      }, character(1)),
      text = vapply(grp, function(g) paste(one[g], collapse = "<br>"), character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    eff_all <- numeric(n)
    pts <- NULL
  }
  shapes <- c(AO = 21, LS = 22, TC = 24, Several = 23)
  add_points <- function(p, yvar, size, legend) {
    if (is.null(pts)) return(p)
    p + suppressWarnings(ggplot2::geom_point(
      data = pts, ggplot2::aes(x = .data$t, y = .data[[yvar]], shape = .data$type,
                               text = .data$text),
      fill = "#B2182B", colour = "#B2182B", size = size)) +
      ggplot2::scale_shape_manual(values = shapes, name = NULL,
                                  guide = if (legend) "legend" else "none")
  }

  if (is.null(beta)) {
    df <- data.frame(t = tt, value = as.numeric(y),
                     text = sprintf("%s<br>%s: %.4f", format(tt), ylab, as.numeric(y)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$value)) +
      suppressWarnings(ggplot2::geom_line(ggplot2::aes(text = .data$text, group = 1),
                                          colour = "grey35", linewidth = 0.6)) +
      ggplot2::labs(title = title, subtitle = "Specified intervention dates",
                    x = NULL, y = ylab) +
      .at_theme()
    p <- add_points(p, "y", 2.6, TRUE)
    return(if (interactive) .as_interactive(p) else p)
  }

  series <- rbind(
    data.frame(t = tt, value = as.numeric(y), series = "Observed"),
    data.frame(t = tt, value = as.numeric(y) - eff_all, series = "Adjusted")
  )
  series$series <- factor(series$series, levels = c("Observed", "Adjusted"))
  series$text <- sprintf("%s<br>%s: %.4f", format(series$t), series$series, series$value)
  p1 <- ggplot2::ggplot(series, ggplot2::aes(x = .data$t, y = .data$value)) +
    suppressWarnings(ggplot2::geom_line(
      ggplot2::aes(colour = .data$series, group = .data$series, text = .data$text),
      linewidth = 0.6)) +
    ggplot2::scale_colour_manual(values = c(Observed = "grey55", Adjusted = "#1F5A96"),
                                 name = NULL) +
    ggplot2::labs(title = paste0("Observed and ", word, "-Adjusted Series"),
                  subtitle = if (!have) none_label, x = NULL, y = ylab) +
    .at_theme()
  p1 <- add_points(p1, "y", 2.6, TRUE)

  edf <- data.frame(t = tt, eff = eff_all, what = "Combined",
                    text = sprintf("%s<br>Combined effect: %.4f", format(tt), eff_all))
  ind <- NULL
  if (individual && have && ncol(Xm) > 1) {
    ind <- do.call(rbind, lapply(seq_len(ncol(Xm)), function(j) {
      e <- Xm[, j] * beta[j]
      data.frame(t = tt, eff = e, what = colnames(Xm)[j],
                 text = sprintf("%s<br>%s: %.4f", format(tt), colnames(Xm)[j], e))
    }))
  }
  labs_all <- c("Combined", if (!is.null(ind)) colnames(Xm))
  pal <- c("#1F5A96", c("#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                        "#999999")[(seq_len(length(labs_all) - 1) - 1) %% 6 + 1])
  names(pal) <- labs_all
  edf$what <- factor(edf$what, levels = labs_all)
  # Combined effect as a solid line; individual effects (optional) as thin
  # dashed lines on top, so that both remain visible where they coincide.
  p2 <- ggplot2::ggplot(edf, ggplot2::aes(x = .data$t, y = .data$eff)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    suppressWarnings(ggplot2::geom_line(
      ggplot2::aes(text = .data$text, group = .data$what, colour = .data$what),
      linewidth = 0.9))
  if (!is.null(ind)) {
    ind$what <- factor(ind$what, levels = labs_all)
    p2 <- p2 + suppressWarnings(ggplot2::geom_line(
      data = ind, ggplot2::aes(text = .data$text, group = .data$what, colour = .data$what),
      linewidth = 0.5, linetype = "dashed"))
  }
  p2 <- p2 +
    ggplot2::scale_colour_manual(values = pal, name = NULL, drop = FALSE,
                                 guide = if (is.null(ind)) "none" else "legend") +
    ggplot2::labs(title = paste0("Combined ", word, " Effect"), x = NULL, y = "Effect") +
    .at_theme()
  p2 <- add_points(p2, "eff", 2.2, FALSE)

  .stack_panels(list(p1, p2), heights = c(0.6, 0.4), title = title,
                interactive = interactive)
}
