# Outliers: regressors, detection (via tsoutliers) and effect plots ---------
#
# Three outlier types, as in Chen and Liu (1993) and tsoutliers, for an
# outlier at observation t0:
#   AO (additive outlier):   1 at t0, 0 elsewhere
#   LS (level shift):        1 from t0 onwards
#   TC (temporary change):   delta^(t - t0) from t0 onwards, 0 before
#
# detect_outliers() currently wraps tsoutliers::tso(). Its arguments and
# output are defined independently of tsoutliers so that the detection
# algorithm can later be replaced by our own implementation (shared with the
# Python package) without changing code that uses it.

.outlier_types <- c("AO", "LS", "TC")
.outlier_pattern <- "^(AO|LS|TC)"

# ---- Periods: parsing and labels ---------------------------------------------

# Labels like 2020Q2 (quarterly), 2008M09 (monthly), 1990 (annual), or the
# observation number for other series.
.period_label <- function(y, idx) {
  if (!stats::is.ts(y)) return(as.character(idx))
  f <- stats::frequency(y)
  tt <- as.numeric(stats::time(y))[idx]
  year <- floor(tt + 1e-8)
  per <- round((tt - year) * f) + 1
  if (f == 4 && min(year) >= 1000) return(sprintf("%dQ%d", year, per))
  if (f == 12 && min(year) >= 1000) return(sprintf("%dM%02d", year, per))
  if (f == 1 && min(year) >= 1000) return(sprintf("%d", year))
  as.character(idx)
}

# Convert user-supplied dates to observation numbers of y. Calendar formats
# are accepted only where they are unambiguous:
#   quarterly ts (frequency 4):  "2020Q2", "2020:2", Date / "2020-04-01"
#   monthly ts (frequency 12):   "2008M09", "2008:9", Date / "2008-09-01"
#   annual ts (frequency 1):     1990, "1990", Date / "1990-01-01"
# Decimal time (e.g. 2020.25) works for any ts if it matches an observation
# exactly. For other frequencies, use `positions` (observation numbers).
.periods_to_index <- function(y, dates) {
  n <- length(y)
  if (!stats::is.ts(y)) {
    if (!is.numeric(dates) || any(dates != round(dates))) {
      stop("For a series without time attributes, give `positions` ",
           "(observation numbers) instead of `dates`.", call. = FALSE)
    }
    idx <- as.integer(dates)
  } else {
    f <- stats::frequency(y)
    tt <- as.numeric(stats::time(y))
    calendar <- f %in% c(1, 4, 12)
    no_calendar <- function(s) {
      stop("Calendar dates such as '", s, "' are supported for annual, ",
           "quarterly and monthly series only. Use `positions` (observation ",
           "numbers) or decimal time for this series (frequency ", f, ").",
           call. = FALSE)
    }
    from_period <- function(yr, p, s) {
      if (p < 1 || p > f) {
        stop("Invalid period in '", s, "': must lie between 1 and ", f, ".",
             call. = FALSE)
      }
      yr + (p - 1) / f
    }
    from_date <- function(d, s) {
      if (!calendar) no_calendar(s)
      yr <- as.numeric(format(d, "%Y")); m <- as.numeric(format(d, "%m"))
      yr + floor((m - 1) * f / 12) / f
    }
    to_time <- function(d) {
      if (inherits(d, "Date")) return(from_date(d, format(d)))
      if (is.numeric(d)) {
        if (f == 1 || d != round(d)) return(d)
        stop("A bare year such as ", d, " is ambiguous for a series with ",
             "frequency ", f, ". Use e.g. \"", d, "Q1\" or \"", d, "M01\".",
             call. = FALSE)
      }
      s <- trimws(as.character(d))
      if (grepl("^[0-9]{4}[Qq][0-9]{1,2}$", s)) {
        if (f != 4) stop("'", s, "' is a quarterly date but the series has ",
                         "frequency ", f, ".", call. = FALSE)
        return(from_period(as.numeric(substr(s, 1, 4)), as.numeric(substring(s, 6)), s))
      }
      if (grepl("^[0-9]{4}[Mm][0-9]{1,2}$", s)) {
        if (f != 12) stop("'", s, "' is a monthly date but the series has ",
                          "frequency ", f, ".", call. = FALSE)
        return(from_period(as.numeric(substr(s, 1, 4)), as.numeric(substring(s, 6)), s))
      }
      if (grepl("^[0-9]{4}:[0-9]{1,3}$", s)) {
        if (!calendar) no_calendar(s)
        return(from_period(as.numeric(substr(s, 1, 4)), as.numeric(substring(s, 6)), s))
      }
      if (grepl("^[0-9]{4}-[0-9]{1,2}(-[0-9]{1,2})?$", s)) {
        s2 <- if (grepl("^[0-9]{4}-[0-9]{1,2}$", s)) paste0(s, "-01") else s
        d2 <- tryCatch(as.Date(s2), error = function(e) NA)
        if (is.na(d2)) stop("Invalid date '", s, "'.", call. = FALSE)
        return(from_date(d2, s))
      }
      if (grepl("^[0-9]{4}$", s)) return(to_time(as.numeric(s)))
      stop("Cannot interpret the date '", s, "'. Use e.g. \"2020Q2\", ",
           "\"2008M09\", \"2008-09-01\" or 1990.", call. = FALSE)
    }
    tq <- if (inherits(dates, "Date")) {
      vapply(seq_along(dates), function(i) to_time(dates[i]), numeric(1))
    } else {
      vapply(as.list(dates), to_time, numeric(1))
    }
    idx <- round((tq - tt[1]) * f) + 1
    bad <- is.na(idx) | idx < 1 | idx > n
    bad[!bad] <- abs(tt[idx[!bad]] - tq[!bad]) > 1e-6
    if (any(bad)) {
      stop("Date(s) outside the sample or not matching an observation: ",
           paste(format(dates[bad]), collapse = ", "), call. = FALSE)
    }
  }
  if (any(idx < 1 | idx > n)) {
    stop("Positions must lie between 1 and ", n, ".", call. = FALSE)
  }
  as.integer(idx)
}

.outlier_column <- function(type, idx, n, delta) {
  x <- numeric(n)
  if (type == "AO") x[idx] <- 1
  if (type == "LS") x[idx:n] <- 1
  if (type == "TC") x[idx:n] <- delta^(0:(n - idx))
  x
}

#' Outlier regressors for chosen dates
#'
#' Builds additive-outlier (AO), level-shift (LS) and temporary-change (TC)
#' regressors at dates you choose, for example from known historical events,
#' to be passed as `xreg` to [forecast::Arima()] or [forecast::auto.arima()].
#' For an outlier at observation t0, AO is 1 at t0 and 0 elsewhere; LS is 1
#' from t0 onwards; TC is `delta^(t - t0)` from t0 onwards and 0 before.
#'
#' Columns are named by type and period (`AO2020Q2`, `LS2008M09`, `TC1990`,
#' or the observation number for series without calendar time), so that
#' [coefficient_diagnostics()], [joint_significance()] and [plot_outliers()]
#' identify them.
#'
#' @param y The series (preferably a `ts` object).
#' @param dates Outlier dates, checked against the frequency of `y`.
#'   Quarterly series: `"2020Q2"` or `"2020:2"`; monthly: `"2008M09"` or
#'   `"2008:9"`; annual: `1990`. For these three frequencies, `Date` objects or
#'   strings like `"2008-09-01"` also work (the period containing the date).
#'   Decimal time (`2020.25`) works for any `ts` if it matches an observation.
#'   For other frequencies (weekly, daily), use `positions`.
#' @param types Outlier types, `"AO"`, `"LS"` or `"TC"`, recycled to the
#'   number of dates.
#' @param delta Decay rate of temporary changes, in (0, 1).
#' @param positions Alternative to `dates`: observation numbers.
#' @return A matrix with one column per outlier: a `ts` matrix with the time
#'   attributes of `y` when `y` is a `ts`, an ordinary matrix otherwise.
#'   Attribute `outliers` holds a table of type, period and position.
#' @examples
#' set.seed(1)
#' y <- ts(cumsum(rnorm(100)), start = c(2000, 1), frequency = 4)
#' X <- outlier_regressors(y, dates = c("2008Q4", "2020Q2"), types = c("LS", "AO"))
#' fit <- forecast::Arima(y, order = c(1, 1, 0), xreg = X)
#' coefficient_diagnostics(fit)
#' @export
outlier_regressors <- function(y, dates = NULL, types = "AO", delta = 0.7,
                               positions = NULL) {
  if (is.null(dates) == is.null(positions)) {
    stop("Give exactly one of `dates` or `positions`.", call. = FALSE)
  }
  if (!is.null(dim(y)) && NCOL(y) > 1) {
    stop("`y` must be a univariate series.", call. = FALSE)
  }
  if (!length(y)) stop("`y` is empty.", call. = FALSE)
  .check_level(delta, "delta")
  given <- if (is.null(positions)) dates else positions
  if (!length(given) || anyNA(given)) {
    stop("Give at least one date or position, without missing values.", call. = FALSE)
  }
  if (!length(types) || anyNA(types)) {
    stop("`types` must contain \"AO\", \"LS\" or \"TC\".", call. = FALSE)
  }
  idx <- if (!is.null(positions)) {
    if (!is.numeric(positions) || any(!is.finite(positions)) ||
        any(positions != round(positions))) {
      stop("`positions` must be observation numbers.", call. = FALSE)
    }
    .periods_to_index(as.numeric(y), positions)
  } else {
    .periods_to_index(y, dates)
  }
  types <- toupper(as.character(types))
  if (!all(types %in% .outlier_types)) {
    stop("`types` must be \"AO\", \"LS\" or \"TC\".", call. = FALSE)
  }
  if (length(types) != length(idx)) {
    if (length(idx) %% length(types) != 0) {
      stop("`types` must have length 1 or the number of dates.", call. = FALSE)
    }
    types <- rep_len(types, length(idx))
  }
  n <- length(y)
  labels <- .period_label(y, idx)
  X <- vapply(seq_along(idx), function(i) .outlier_column(types[i], idx[i], n, delta),
              numeric(n))
  X <- matrix(X, nrow = n)
  colnames(X) <- paste0(types, labels)
  if (anyDuplicated(colnames(X))) {
    stop("The same outlier (type and date) is given twice.", call. = FALSE)
  }
  if (stats::is.ts(y)) {
    X <- stats::ts(X, start = stats::start(y), frequency = stats::frequency(y))
  }
  attr(X, "outliers") <- data.frame(type = types, date = labels, index = idx,
                                    stringsAsFactors = FALSE)
  attr(X, "delta") <- delta
  X
}

#' Automatic outlier detection (Chen-Liu procedure via tsoutliers)
#'
#' Detects additive outliers (AO), level shifts (LS) and temporary changes (TC)
#' in a nonseasonal ARIMA model with the iterative procedure of Chen and Liu
#' (1993), as implemented in [tsoutliers::tso()], and re-estimates the final
#' model with [forecast::Arima()] and outlier regressors built by
#' [outlier_regressors()].
#'
#' Three ways to specify the model:
#' * `detect_outliers(fit)`: a fitted nonseasonal model; its series, order and
#'   constant (mean or drift) are kept fixed. Models with fixed (constrained)
#'   coefficients are rejected, because the constraints would not be kept.
#' * `detect_outliers(y, order = c(1, 0, 1))`: the same without fitting first.
#' * `detect_outliers(y)`: the order is chosen by `auto.arima` with the course
#'   defaults (BIC, nonseasonal). The order selected in the presence of
#'   outliers can differ from the one chosen for the raw series.
#'
#' **Supported models.** Univariate, nonseasonal ARIMA models of the series as
#' given, with an optional mean (d = 0) or drift (d = 1). Drift means the
#' built-in trend of `forecast::Arima(..., include.drift = TRUE)`; a
#' user-supplied regressor named `drift` (even a shifted trend) is rejected. Transform the series
#' yourself (logs, growth rates) before calling the function. Models with other
#' explanatory variables (`xreg`) or a Box-Cox transformation (`lambda`) are
#' rejected, because the outlier effects and the adjusted series would then be
#' defined on a different scale or conditional on other variables.
#'
#' **Interpretation.** AO, LS and TC are estimated intervention components of
#' the modelled series, not proof of a cause. The reported standard errors and
#' z-statistics treat the selected dates, types and model order as known; they
#' are not adjusted for the search that selected them. The table comes from the
#' final fit, so it matches [coefficient_diagnostics()] and can differ in the
#' last decimals from the numbers printed by `tso()`.
#'
#' **Warnings and convergence.** Warnings from `tso()` (for example, iteration
#' limits reached during the search) are passed on as a single warning and kept
#' in `$warnings`. If the final model did not converge, a warning is issued and
#' `$converged` is `FALSE`.
#'
#' **Collinearity.** Rank checking is switched on in `tso()`, so a candidate
#' outlier that is numerically collinear with others is dropped instead of
#' stopping the search. Collinear cases include an LS and a TC on the last
#' observation (identical regressors), and an LS on the first observation
#' together with a mean (or with differencing). An LS and a TC on the same
#' interior date are distinct regressors, but become nearly collinear as
#' `delta` approaches 1.
#'
#' Requires the 'tsoutliers' package.
#'
#' @param x A fitted model ([forecast::Arima()], [stats::arima()]) or a series.
#' @param order `c(p, d, q)` when `x` is a series; `NULL` lets auto.arima
#'   choose it.
#' @param types Outlier types to search for.
#' @param cval Critical value for the t-statistics. `NULL` uses the tsoutliers
#'   rule: 3 for n <= 50, 4 for n >= 450, linear in between.
#' @param delta Decay rate of temporary changes.
#' @param include.constant With a fixed order: `NULL` (mean if d = 0), `TRUE`
#'   (mean if d = 0, drift if d = 1) or `FALSE`. Ignored when auto.arima
#'   chooses the model.
#' @param ... Further arguments for [tsoutliers::tso()] that control the search
#'   (e.g. `maxit.iloop`, `discard.method`). Names are matched as `tso()`
#'   would match them (so `xr` means `xreg`); `xreg`, `lambda` and the
#'   arguments that specify the model are not allowed here.
#' @return An object of class `arimatools_outliers`: a list with `outliers`
#'   (table with `type`, `date`, `index`, `time`, `estimate`, `std.error`,
#'   `statistic`), `fit` (final `Arima` model), `regressors`, `effects`
#'   (combined outlier effect, a `ts`), `adjusted` (series minus effects),
#'   `y`, `order`, `cval`, `delta`, `auto` (was the order selected?),
#'   `converged` (did the final fit converge?) and `warnings`.
#' @export
detect_outliers <- function(x, order = NULL, types = c("AO", "LS", "TC"),
                            cval = NULL, delta = 0.7, include.constant = NULL,
                            ...) {
  if (!requireNamespace("tsoutliers", quietly = TRUE)) {
    stop("detect_outliers() needs the 'tsoutliers' package: ",
         "install.packages(\"tsoutliers\").", call. = FALSE)
  }
  types <- toupper(as.character(types))
  if (!length(types) || !all(types %in% .outlier_types)) {
    stop("`types` must be a subset of \"AO\", \"LS\", \"TC\".", call. = FALSE)
  }
  .check_level(delta, "delta")
  dots <- list(...)
  # Resolve partial names (e.g. `xr` -> `xreg`) as tso() itself would.
  dot_names <- .resolve_dot_names(dots, tsoutliers::tso)
  if ("lambda" %in% names(dots)) dot_names[names(dots) == "lambda"] <- "lambda"
  unknown <- setdiff(dot_names, c(names(formals(tsoutliers::tso)), "lambda"))
  if (length(unknown)) {
    stop("Unknown argument(s) in `...`: ", paste(unknown, collapse = ", "),
         ". See ?tsoutliers::tso for the search options.", call. = FALSE)
  }
  names(dots) <- dot_names
  reserved <- c("y", "xreg", "lambda", "types", "cval", "delta", "tsmethod",
                "args.tsmethod", "check.rank")
  bad_dots <- intersect(dot_names, reserved)
  if (length(bad_dots)) {
    if (any(c("xreg", "lambda") %in% bad_dots)) {
      stop("detect_outliers() supports univariate models of the series as given: ",
           "external regressors (`xreg`) and Box-Cox transformations (`lambda`) ",
           "are not supported. Transform the series yourself before calling it.",
           call. = FALSE)
    }
    stop("These arguments cannot be passed through `...`: ",
         paste(bad_dots, collapse = ", "), ".", call. = FALSE)
  }

  # ---- Model specification -----------------------------------------------------
  if (inherits(x, "Arima")) {
    if (!is.null(order)) stop("Give either a fitted model or `order`, not both.",
                              call. = FALSE)
    o <- .orders(x)
    if (o$P + o$D + o$Q > 0) {
      stop("detect_outliers() supports nonseasonal models only.", call. = FALSE)
    }
    if (!is.null(x$lambda)) {
      stop("The model uses a Box-Cox transformation (lambda), which ",
           "detect_outliers() does not support. Transform the series yourself ",
           "(e.g. log(y)) and fit the model without `lambda`.", call. = FALSE)
    }
    if (any(!.free(x))) {
      stop("The model has fixed (constrained) coefficients, which ",
           "detect_outliers() does not preserve. Fit the model without `fixed`.",
           call. = FALSE)
    }
    cn <- names(x$coef)
    other <- cn[!grepl("^(ar|ma)[0-9]+$", cn) & !cn %in% c("intercept", "drift")]
    if (length(other)) {
      stop("The model has explanatory variables (", paste(other, collapse = ", "),
           "), which detect_outliers() does not support. Fit the model without ",
           "them.", call. = FALSE)
    }
    y <- .response(x)
    if (!all(is.finite(y))) {
      stop("The series must not contain missing or infinite values.", call. = FALSE)
    }
    if ("drift" %in% cn && !.is_builtin_drift(x)) {
      stop("The model has a regressor named 'drift' that is not the built-in ",
           "drift. Only the trend created by forecast::Arima(..., include.drift = ",
           "TRUE) is supported; other regressors, including shifted or rescaled ",
           "trends, are not.", call. = FALSE)
    }
    ord <- c(o$p, o$d, o$q)
    use_mean <- "intercept" %in% cn
    use_drift <- "drift" %in% cn
    auto <- FALSE
  } else {
    if (!is.null(dim(x)) && NCOL(x) > 1) {
      stop("`x` must be a univariate series.", call. = FALSE)
    }
    y <- x
    if (!is.numeric(y)) stop("`x` must be a numeric series or a fitted model.",
                             call. = FALSE)
    if (!stats::is.ts(y)) y <- stats::ts(as.numeric(y))
    if (!all(is.finite(y))) {
      stop("The series must not contain missing or infinite values.", call. = FALSE)
    }
    auto <- is.null(order)
    if (!auto) {
      if (!is.numeric(order) || length(order) != 3L ||
          any(!vapply(order, .is_count, logical(1)))) {
        stop("`order` must be c(p, d, q) with non-negative integers.", call. = FALSE)
      }
      ord <- as.integer(order)
      d <- ord[2]
      use_mean <- if (is.null(include.constant)) d == 0 else isTRUE(include.constant) && d == 0
      use_drift <- isTRUE(include.constant) && d == 1
    } else if (!is.null(include.constant)) {
      warning("`include.constant` is ignored when auto.arima chooses the model.",
              call. = FALSE)
    }
  }
  n <- length(y)

  # ---- Sample size -------------------------------------------------------------
  if (auto) {
    if (n < 20) {
      stop("Too few observations (", n, ") for automatic model selection with ",
           "outlier detection; at least 20 are needed.", call. = FALSE)
    }
  } else {
    k <- ord[1] + ord[3] + use_mean + use_drift
    if (n - ord[2] - k < 10) {
      stop("Too few observations (", n, ") for ARIMA(", paste(ord, collapse = ","),
           ") with outlier detection: after differencing, at least 10 more ",
           "observations than estimated coefficients are needed.", call. = FALSE)
    }
  }

  # ---- Detection with tsoutliers::tso -----------------------------------------
  tso_args <- c(list(y = y, types = types, cval = cval, delta = delta,
                     check.rank = TRUE), dots)
  if (auto) {
    tso_args$tsmethod <- "auto.arima"
    tso_args$args.tsmethod <- list(ic = "bic", seasonal = FALSE,
                                   approximation = FALSE)
  } else {
    tso_args$tsmethod <- "arima"
    tso_args$args.tsmethod <- list(order = ord, include.mean = use_mean)
    if (use_drift) tso_args$xreg <- cbind(drift = seq_len(n))
  }
  det <- .collect_warnings(tryCatch(
    do.call(tsoutliers::tso, tso_args),
    error = function(e) {
      stop("tsoutliers::tso() failed: ", conditionMessage(e), call. = FALSE)
    }
  ))
  res <- det$value
  if (auto) {
    a <- res$fit$arma
    ord <- c(a[1], a[6], a[2])
    cn <- names(res$fit$coef)
    use_mean <- "intercept" %in% cn
    use_drift <- "drift" %in% cn
  }

  # ---- Final model with our own regressors ------------------------------------
  out <- res$outliers
  if (!is.null(out) && nrow(out)) {
    out <- out[order(out$ind), , drop = FALSE]
    X <- outlier_regressors(y, positions = out$ind, types = as.character(out$type),
                            delta = delta)
  } else {
    X <- NULL
  }
  final <- .collect_warnings(forecast::Arima(
    y, order = ord, xreg = X, include.mean = use_mean, include.drift = use_drift
  ))
  fit <- final$value
  fit$series <- if (inherits(x, "Arima")) .series_name(x) else deparse1(substitute(x))
  converged <- .converged(fit)

  warns <- unique(c(det$warnings, final$warnings))
  if (length(warns)) {
    warning("Warnings during outlier detection:\n  ",
            paste(warns, collapse = "\n  "),
            "\nThe results may still be valid; see `$warnings` and check the fit.",
            call. = FALSE)
  }
  if (!converged) {
    warning("The final model did not converge (optimizer code ", fit$code, "). ",
            "Treat the estimates with caution; a simpler model may help.",
            call. = FALSE)
  }

  if (!is.null(X)) {
    cd <- coefficient_diagnostics(fit)
    cd <- cd[match(colnames(X), cd$term), ]
    tab <- data.frame(
      type = attr(X, "outliers")$type, date = attr(X, "outliers")$date,
      index = attr(X, "outliers")$index,
      time = .time_values(y)[attr(X, "outliers")$index],
      estimate = cd$estimate, std.error = cd$std.error, statistic = cd$statistic,
      stringsAsFactors = FALSE
    )
    eff <- as.numeric(X %*% fit$coef[colnames(X)])
  } else {
    tab <- data.frame(type = character(0), date = character(0), index = integer(0),
                      time = .time_values(y)[0], estimate = numeric(0),
                      std.error = numeric(0), statistic = numeric(0))
    eff <- numeric(n)
  }
  eff <- stats::ts(eff, start = stats::start(y), frequency = stats::frequency(y))

  structure(
    list(outliers = tab, fit = fit, regressors = X, effects = eff,
         adjusted = y - eff, y = y, order = ord, cval = res$cval, delta = delta,
         auto = auto, converged = converged, warnings = warns),
    class = "arimatools_outliers"
  )
}

# A built-in drift is the time trend that forecast::Arima() creates itself
# with include.drift = TRUE. To accept only that, the stored column is compared
# with the one forecast builds for the same series (whatever its version
# uses), instead of merely checking that it looks linear: a shifted trend such
# as t + 10 is a different model when there is no mean. stats::arima has no
# built-in drift, so a 'drift' coefficient there always comes from a user
# regressor.
.is_builtin_drift <- function(fit) {
  X <- fit$xreg
  if (is.null(X) || !"drift" %in% colnames(X)) return(FALSE)
  d <- .orders(fit)$d
  if (d > 1) return(FALSE)
  y <- tryCatch(.response(fit), error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ref <- tryCatch(
    suppressWarnings(forecast::Arima(y, order = c(0, d, 0), include.drift = TRUE,
                                     include.mean = FALSE))$xreg,
    error = function(e) NULL)
  if (is.null(ref) || !"drift" %in% colnames(ref)) return(FALSE)
  a <- as.numeric(X[, "drift"]); b <- as.numeric(ref[, "drift"])
  length(a) == length(b) && isTRUE(all.equal(a, b, tolerance = 1e-10))
}

# Evaluate an expression, muffling but recording its warnings.
.collect_warnings <- function(expr) {
  warns <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(warns))
}

#' @export
print.arimatools_outliers <- function(x, digits = 4, ...) {
  cat("Outlier detection (Chen-Liu, via tsoutliers)\n")
  cn <- names(x$fit$coef)
  det <- if ("drift" %in% cn) " with drift" else if ("intercept" %in% cn) " with mean" else ""
  cat("Model:          ", paste0(.order_label(x$fit), det),
      if (x$auto) "(selected by auto.arima, BIC)" else "(order fixed)", "\n")
  cat("Critical value: ", format(x$cval), "\n")
  cat("Outliers found: ", nrow(x$outliers), "\n")
  if (isFALSE(x$converged)) cat("Converged:       NO (see warnings)\n")
  if (length(x$warnings)) cat("Warnings:       ", length(x$warnings), "(see $warnings)\n")
  if (nrow(x$outliers)) {
    cat("\n")
    tab <- x$outliers[, c("type", "date", "index", "estimate", "std.error", "statistic")]
    print(tab, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' Plot outlier effects
#'
#' Two panels in the style of `tsoutliers::plot.tsoutliers()`: the observed and
#' outlier-adjusted series, with outlier dates marked by type, and the combined
#' outlier effect over time.
#'
#' Works both for automatic detection and for outliers chosen by hand:
#' `x` can be the output of [detect_outliers()] or any model fitted with
#' regressors from [outlier_regressors()]. In the second case the effect is the
#' regressors multiplied by their estimated coefficients. Several outliers on
#' the same date share one marker (a diamond), whose tooltip lists them all.
#'
#' Models with a Box-Cox transformation (`lambda`) are rejected, because their
#' effects are on the transformed scale.
#'
#' @param x Output of [detect_outliers()], or a fitted `Arima` model whose
#'   regressors include columns named by [outlier_regressors()].
#' @param title Overall title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A patchwork object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_outliers <- function(x, title = "Outlier Effects", interactive = FALSE) {
  if (inherits(x, "arimatools_outliers")) {
    y <- x$y; X <- x$regressors; fit <- x$fit
  } else if (inherits(x, "Arima")) {
    fit <- x
    if (!is.null(fit$lambda)) {
      stop("The model uses a Box-Cox transformation (lambda): its outlier effects ",
           "are on the transformed scale and cannot be subtracted from the ",
           "original series. Transform the series yourself (e.g. log(y)) and fit ",
           "the model without `lambda`.", call. = FALSE)
    }
    y <- .response(fit)
    X <- fit$xreg
    keep <- grepl(.outlier_pattern, colnames(X))
    if (is.null(X) || !any(keep)) {
      stop("The model has no outlier regressors (columns named AO..., LS... or ",
           "TC..., as created by outlier_regressors()).", call. = FALSE)
    }
    X <- X[, keep, drop = FALSE]
  } else {
    stop("`x` must be the output of detect_outliers() or a fitted Arima model.",
         call. = FALSE)
  }

  n <- length(y)
  tt <- .time_values(y)
  if (!is.null(X) && ncol(X)) {
    eff <- as.numeric(as.matrix(X) %*% fit$coef[colnames(X)])
    first <- apply(as.matrix(X), 2, function(col) which(col != 0)[1])
    one <- sprintf("%s: %.4f", colnames(X), fit$coef[colnames(X)])
    # One marker per date; several outliers on the same date share a marker.
    pos <- sort(unique(first))
    grp <- lapply(pos, function(i) which(first == i))
    pts <- data.frame(
      t = tt[pos], y = as.numeric(y)[pos], eff = eff[pos],
      type = vapply(grp, function(g) {
        tp <- unique(substr(colnames(X)[g], 1, 2))
        if (length(tp) > 1) "Several" else tp
      }, character(1)),
      text = vapply(grp, function(g) paste(one[g], collapse = "<br>"), character(1)),
      stringsAsFactors = FALSE
    )
    subtitle <- NULL
  } else {
    eff <- numeric(n)
    pts <- NULL
    subtitle <- "No outliers detected"
  }

  series <- rbind(
    data.frame(t = tt, value = as.numeric(y), series = "Observed"),
    data.frame(t = tt, value = as.numeric(y) - eff, series = "Adjusted")
  )
  series$series <- factor(series$series, levels = c("Observed", "Adjusted"))
  series$text <- sprintf("%s<br>%s: %.4f", format(series$t), series$series, series$value)

  p1 <- ggplot2::ggplot(series, ggplot2::aes(x = .data$t, y = .data$value)) +
    suppressWarnings(ggplot2::geom_line(
      ggplot2::aes(colour = .data$series, group = .data$series, text = .data$text),
      linewidth = 0.6)) +
    ggplot2::scale_colour_manual(values = c(Observed = "grey55", Adjusted = "#1F5A96"),
                                 name = NULL) +
    ggplot2::labs(title = "Observed and Outlier-Adjusted Series", subtitle = subtitle,
                  x = NULL, y = .series_name(fit)) +
    .at_theme()

  edf <- data.frame(t = tt, eff = eff,
                    text = sprintf("%s<br>Outlier effect: %.4f", format(tt), eff))
  p2 <- ggplot2::ggplot(edf, ggplot2::aes(x = .data$t, y = .data$eff)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    suppressWarnings(ggplot2::geom_line(ggplot2::aes(text = .data$text, group = 1),
                                        colour = "#1F5A96", linewidth = 0.6)) +
    ggplot2::labs(title = "Combined Outlier Effect", x = NULL, y = "Effect") +
    .at_theme()

  if (!is.null(pts)) {
    shapes <- c(AO = 21, LS = 22, TC = 24, Several = 23)
    p1 <- p1 + suppressWarnings(ggplot2::geom_point(
      data = pts, ggplot2::aes(x = .data$t, y = .data$y, shape = .data$type,
                               text = .data$text),
      fill = "#B2182B", colour = "#B2182B", size = 2.6)) +
      ggplot2::scale_shape_manual(values = shapes, name = NULL)
    p2 <- p2 + suppressWarnings(ggplot2::geom_point(
      data = pts, ggplot2::aes(x = .data$t, y = .data$eff, shape = .data$type,
                               text = .data$text),
      fill = "#B2182B", colour = "#B2182B", size = 2.2)) +
      ggplot2::scale_shape_manual(values = shapes, guide = "none")
  }

  .stack_panels(list(p1, p2), heights = c(0.6, 0.4), title = title,
                interactive = interactive)
}
