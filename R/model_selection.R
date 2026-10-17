# Model comparison, order searches and automatic ARIMA selection -------------

#' Compare fitted ARIMA models by likelihood and information criteria
#'
#' With k = number of free coefficients + 1 (innovation variance) and n the
#' number of observations used (after differencing):
#' AIC = -2 logL + 2k, AICc = AIC + 2k(k+1)/(n-k-1), BIC = -2 logL + k log n,
#' HQIC = -2 logL + 2k log log n.
#'
#' Criteria are comparable only across models for the same dependent variable
#' and the same effective sample (in particular, the same d). Fits estimated by
#' conditional sum of squares (`method = "CSS"`) have no likelihood, so their
#' criteria are reported as `NA`.
#'
#' @param ... Two or more fitted ARIMA models (optionally named).
#' @param names Optional model labels; overrides argument names.
#' @return A data frame with columns `model`, `order`, `n`, `k`, `loglik`,
#'   `AIC`, `AICc`, `BIC` and `HQIC` (row names = model labels).
#' @export
compare_models <- function(..., names = NULL) {
  fits <- list(...)
  if (length(fits) < 2L) stop("At least two fitted models must be supplied.", call. = FALSE)
  lapply(fits, .check_fit)
  labels <- vapply(fits, .order_label, character(1))
  if (!is.null(base::names(fits))) {
    has <- nzchar(base::names(fits))
    labels[has] <- base::names(fits)[has]
  }
  if (!is.null(names)) {
    if (length(names) != length(fits)) {
      stop("The number of names must equal the number of models.", call. = FALSE)
    }
    labels <- as.character(names)
  }
  if (anyDuplicated(labels)) {
    stop("Model names must be unique (use `names =`).", call. = FALSE)
  }
  rows <- lapply(fits, function(f) {
    ic <- .ic_values(f)
    data.frame(order = .order_string(f), n = ic$n, k = ic$k, loglik = ic$loglik,
               AIC = ic$AIC, AICc = ic$AICc, BIC = ic$BIC, HQIC = ic$HQIC,
               stringsAsFactors = FALSE)
  })
  out <- cbind(model = labels, do.call(rbind, rows), stringsAsFactors = FALSE)
  rownames(out) <- labels
  out
}

#' Grid search over nonseasonal ARIMA(p,d,q) orders
#'
#' Estimates every ARIMA(p,d,q) with p = 0..p_max, q = 0..q_max and
#' d = d_min..d_max by [forecast::Arima()], and reports all four information
#' criteria. It estimates and reports; choosing a model is left to
#' [select_model()]. Rows are grouped by d and, within each d, sorted by
#' `criterion` with converged models first; models that did not converge or
#' failed stay visible.
#'
#' No unit-root test is performed: the range of d is the user's choice. A
#' range of d is allowed for inspection, but criteria are comparable only
#' within the same d, and [select_model()] requires `d =` for such tables.
#' With `method = "CSS"` there is no likelihood, so the criteria are `NA`.
#'
#' Related: `forecast::auto.arima(y, stepwise = FALSE, approximation = FALSE,
#' trace = TRUE)` prints a similar search but does not return the table.
#'
#' @param y A numeric vector or `ts` object.
#' @param p_max,q_max Maximum AR and MA orders.
#' @param d_min,d_max Range of integration orders.
#' @param criterion One of "AIC", "AICc", "BIC", "HQIC" (case-insensitive).
#' @param include.constant See [recursive_estimation()].
#' @param ... Further arguments for [forecast::Arima()] (e.g. `method = "ML"`).
#' @return A data frame with columns `model`, `p`, `d`, `q`, `n`, `k`,
#'   `loglik`, `AIC`, `AICc`, `BIC`, `HQIC`, `converged` and `status`.
#' @export
order_search <- function(y, p_max = 4, q_max = 4, d_min = 0, d_max = 0,
                         criterion = "BIC", include.constant = NULL, ...) {
  p_max <- .check_count(p_max, "p_max")
  q_max <- .check_count(q_max, "q_max")
  d_min <- .check_count(d_min, "d_min")
  d_max <- .check_count(d_max, "d_max")
  if (d_min > d_max) stop("`d_min` must be <= `d_max`.", call. = FALSE)
  crit <- .match_criterion(criterion)

  grid <- expand.grid(q = 0:q_max, p = 0:p_max, d = d_min:d_max)[, c("p", "d", "q")]
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    o <- c(grid$p[i], grid$d[i], grid$q[i])
    lab <- sprintf("ARIMA(%d,%d,%d)", o[1], o[2], o[3])
    res <- tryCatch({
      fit <- .fit_arima(y, o, include.constant, ...)
      ic <- .ic_values(fit)
      conv <- .converged(fit)
      status <- if (.is_css(fit)) "CSS fit: no likelihood, criteria not available"
        else if (!conv) "Not converged" else if (isTRUE(attr(fit, "ml_fallback")))
        "OK (ML after CSS-ML failed)" else "OK"
      data.frame(n = ic$n, k = ic$k, loglik = ic$loglik, AIC = ic$AIC,
                 AICc = ic$AICc, BIC = ic$BIC, HQIC = ic$HQIC,
                 converged = conv, status = status, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(n = NA_integer_, k = NA_integer_, loglik = NA_real_,
                 AIC = NA_real_, AICc = NA_real_, BIC = NA_real_, HQIC = NA_real_,
                 converged = FALSE, status = paste("Failed:", conditionMessage(e)),
                 stringsAsFactors = FALSE)
    })
    cbind(data.frame(model = lab, p = o[1], d = o[2], q = o[3],
                     stringsAsFactors = FALSE), res)
  })
  out <- do.call(rbind, rows)
  # Report, do not rank across d: group by d, then order by the criterion
  # within each d (converged fits first).
  out <- out[order(out$d, !out$converged, out[[crit]], na.last = TRUE), ]
  rownames(out) <- NULL
  attr(out, "criterion") <- crit
  if (d_max > d_min) {
    message("Note: information criteria are comparable only among models with ",
            "the same d (same effective sample). Choose d first; ",
            "select_model() then needs `d =`.")
  }
  out
}

#' Select the model that minimizes an information criterion
#'
#' This is the only function in the package that chooses a model;
#' [order_search()] and [compare_models()] only estimate and report. Before
#' choosing, it checks that the criteria are comparable. Non-converged models
#' (if a `converged` column exists) and non-finite values are excluded. The
#' minimum criterion says nothing about residual adequacy or parameter
#' stability.
#'
#' Information criteria are comparable only among models with the same order
#' of differencing d, because differencing changes the effective sample. The
#' order of integration should be chosen first (unit-root tests, graphs,
#' correlograms). If the table contains models with different d,
#' `select_model()` stops unless `d` says which d to select within. It also
#' stops if the remaining candidates were estimated on different numbers of
#' observations. It cannot check that the dependent variable is the same; that
#' remains the user's responsibility.
#'
#' Fits estimated by conditional sum of squares (`method = "CSS"`) have no
#' likelihood; their criteria are `NA` and they are never selected.
#'
#' @param table Output of [compare_models()] or [order_search()].
#' @param criterion One of "AIC", "AICc", "BIC", "HQIC".
#' @param d Order of differencing to select within. Needed only when the table
#'   contains several values of d.
#' @return The selected row, as a one-row data frame.
#' @export
select_model <- function(table, criterion = "BIC", d = NULL) {
  crit <- .match_criterion(criterion)
  if (!crit %in% names(table)) {
    stop("The table does not contain a '", crit, "' column.", call. = FALSE)
  }
  dd <- if ("d" %in% names(table)) {
    table$d
  } else if ("order" %in% names(table)) {
    as.integer(sub("^\\([0-9]+,([0-9]+),[0-9]+\\)$", "\\1", table$order))
  } else NULL
  keep <- rep(TRUE, nrow(table))
  if (!is.null(d)) {
    d <- .check_count(d, "d")
    if (is.null(dd)) stop("The table has no information on d.", call. = FALSE)
    keep <- !is.na(dd) & dd == d
    if (!any(keep)) stop("No models with d = ", d, " in the table.", call. = FALSE)
  } else if (!is.null(dd) && length(unique(stats::na.omit(dd))) > 1) {
    stop("The table contains models with different orders of differencing (d = ",
         paste(sort(unique(stats::na.omit(dd))), collapse = ", "), "). ",
         "Information criteria are not comparable across d, because the ",
         "effective sample differs. Choose d first and use `d =`, e.g. ",
         "select_model(table, \"", crit, "\", d = 1).", call. = FALSE)
  }
  ok <- keep & is.finite(table[[crit]])
  if ("converged" %in% names(table)) ok <- ok & table$converged
  if (!any(ok)) stop("No eligible models with finite ", crit, " values.", call. = FALSE)
  if ("n" %in% names(table) && length(unique(table$n[ok])) > 1) {
    stop("The candidate models were estimated on different numbers of ",
         "observations (n = ", paste(sort(unique(table$n[ok])), collapse = ", "),
         "). Information criteria are comparable only on the same sample.",
         call. = FALSE)
  }
  idx <- which(ok)[which.min(table[[crit]][ok])]
  table[idx, , drop = FALSE]
}

#' Automatic nonseasonal ARIMA selection (Hyndman-Khandakar)
#'
#' A thin wrapper around [forecast::auto.arima()] with the course defaults:
#' nonseasonal models and BIC. (StatsForecast's AutoARIMA, used by the Python
#' version, is a port of this same function.) Unlike in Python, the result is
#' an ordinary `Arima` object, so every other function in the package accepts
#' it directly.
#'
#' Intended as an automatic benchmark, not a replacement for specification
#' analysis and diagnostics.
#'
#' @param y A numeric vector or `ts` object.
#' @param criterion One of "AIC", "AICc", "BIC".
#' @param d Order of differencing; `NULL` lets auto.arima choose it by
#'   repeated KPSS tests.
#' @param max_p,max_q Maximum AR and MA orders.
#' @param stepwise Use the stepwise search? `FALSE` searches exhaustively.
#' @param allow_mean,allow_drift Allow a non-zero mean / a drift.
#' @param approximation Use the CSS approximation during the search?
#'   Default `FALSE`, as in StatsForecast (forecast's own default switches the
#'   approximation on for series longer than 150 observations, which can
#'   change the selected model).
#' @param ... Further arguments for [forecast::auto.arima()].
#' @return A fitted `forecast_ARIMA`/`Arima` object with attribute `ic`.
#' @export
auto_arima <- function(y, criterion = "BIC", d = NULL, max_p = 5, max_q = 5,
                       stepwise = TRUE, allow_mean = TRUE, allow_drift = TRUE,
                       approximation = FALSE, ...) {
  crit <- .match_criterion(criterion, c("AIC", "AICc", "BIC"))
  if (!is.null(d)) d <- .check_count(d, "d")
  max_p <- .check_count(max_p, "max_p")
  max_q <- .check_count(max_q, "max_q")
  if (any(!is.finite(y))) stop("`y` must contain only finite values.", call. = FALSE)
  fit <- forecast::auto.arima(
    y, d = if (is.null(d)) NA else d, max.p = max_p, max.q = max_q,
    seasonal = FALSE, ic = tolower(crit), stepwise = stepwise,
    allowmean = allow_mean, allowdrift = allow_drift,
    approximation = approximation, ...
  )
  attr(fit, "ic") <- crit
  fit
}

#' Compact summary of an automatically selected ARIMA model
#'
#' Reports the selected order and spells out the deterministic component: for
#' d = 0 a mean is the unconditional mean of the stationary process; for d = 1
#' a drift is a constant in the differenced model, i.e. a linear trend in the
#' level.
#'
#' @param model Output of [auto_arima()] (any fitted `Arima` object works; the
#'   selection criterion is then reported as `NA`).
#' @return A named list of class `arimatools_summary`.
#' @export
auto_arima_summary <- function(model) {
  .check_fit(model, "model")
  o <- .orders(model)
  cn <- names(model$coef)
  if ("drift" %in% cn) {
    det <- "Drift"; interp <- "Constant in differenced model; linear drift in level"
  } else if ("intercept" %in% cn) {
    det <- "Mean"; interp <- "Non-zero unconditional mean"
  } else if (o$d + o$D == 0) {
    det <- "None"; interp <- "Zero unconditional mean"
  } else {
    det <- "None"; interp <- "No deterministic component reported"
  }
  ic <- .ic_values(model)
  aicc <- if (!is.null(model$aicc)) model$aicc else ic$AICc
  bic <- if (!is.null(model$bic)) model$bic else ic$BIC
  crit <- attr(model, "ic")
  vals <- c(AIC = model$aic, AICc = aicc, BIC = bic)
  .at_summary(
    list(
      "Selected model" = .order_label(model),
      "Deterministic component" = det,
      "Interpretation" = interp,
      "Selection criterion" = if (is.null(crit)) NA_character_ else crit,
      "Criterion value" = if (is.null(crit)) NA_real_ else unname(vals[crit]),
      "AIC" = model$aic, "AICc" = aicc, "BIC" = bic,
      "Log likelihood" = model$loglik,
      "Observations" = as.integer(model$nobs),
      "Innovation variance" = model$sigma2
    ),
    "AutoARIMA summary"
  )
}
