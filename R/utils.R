# Internal helpers -------------------------------------------------------------
#
# None of these functions is exported. They centralize the conventions shared
# by the public API: validation of fitted models, ARIMA orders, polynomials,
# information criteria, robust estimation, time axes, and the plot theme.

# ---- Validation --------------------------------------------------------------

.check_fit <- function(fit, arg = "fit") {
  if (!inherits(fit, "Arima")) {
    stop(
      "`", arg, "` must be a fitted ARIMA model from forecast::Arima(), ",
      "forecast::auto.arima() or stats::arima().",
      call. = FALSE
    )
  }
  invisible(fit)
}

.check_level <- function(x, arg = "confidence") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must lie strictly between 0 and 1.", call. = FALSE)
  }
  invisible(x)
}

.is_count <- function(x, min = 0) {
  is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
}

.check_count <- function(x, arg, min = 0) {
  if (!.is_count(x, min)) {
    what <- if (min >= 1) "a positive integer" else "a non-negative integer"
    stop("`", arg, "` must be ", what, ".", call. = FALSE)
  }
  as.integer(x)
}

.match_criterion <- function(criterion, allowed = c("AIC", "AICc", "BIC", "HQIC")) {
  key <- match(tolower(criterion), tolower(allowed))
  if (length(criterion) != 1L || is.na(key)) {
    stop("`criterion` must be one of: ", paste(allowed, collapse = ", "), ".",
         call. = FALSE)
  }
  allowed[key]
}

# ---- Orders and polynomials --------------------------------------------------

# stats::arima stores arma = c(p, q, P, Q, s, d, D)
.orders <- function(fit) {
  a <- fit$arma
  list(p = a[1], q = a[2], P = a[3], Q = a[4], s = a[5], d = a[6], D = a[7])
}

.order_label <- function(fit) {
  o <- .orders(fit)
  lab <- sprintf("ARIMA(%d,%d,%d)", o$p, o$d, o$q)
  if (o$P + o$D + o$Q > 0) {
    lab <- paste0(lab, sprintf("(%d,%d,%d)[%d]", o$P, o$D, o$Q, o$s))
  }
  lab
}

.order_string <- function(fit) {
  o <- .orders(fit)
  sprintf("(%d,%d,%d)", o$p, o$d, o$q)
}

.trim_zeros <- function(x, tol = 0) {
  nz <- which(abs(x) > tol)
  if (length(nz) == 0L) numeric(0) else x[seq_len(max(nz))]
}

# Roots of the (expanded, seasonal-multiplied) AR and MA polynomials,
#   phi(z)   = 1 - phi_1 z - ... - phi_p z^p
#   theta(z) = 1 + theta_1 z + ... + theta_q z^q
.ar_roots <- function(fit) {
  phi <- .trim_zeros(fit$model$phi)
  if (length(phi) == 0L) return(complex(0))
  polyroot(c(1, -phi))
}

.ma_roots <- function(fit) {
  theta <- .trim_zeros(fit$model$theta)
  if (length(theta) == 0L) return(complex(0))
  polyroot(c(1, theta))
}

.conv <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    idx <- i:(i + length(b) - 1L)
    out[idx] <- out[idx] + a[i] * b
  }
  out
}

# psi_0, ..., psi_(n-1) of the full ARIMA transfer function, including ordinary
# and seasonal differencing. fit$model$Delta holds the differencing polynomial
# (1 - B)^d (1 - B^s)^D written as y_t - sum_i Delta_i y_(t-i).
.psi_weights <- function(fit, n) {
  ar_poly <- c(1, -fit$model$phi)
  delta <- fit$model$Delta
  if (length(delta)) ar_poly <- .conv(ar_poly, c(1, -delta))
  phi_full <- -ar_poly[-1]
  theta <- fit$model$theta
  if (n == 1L) return(1)
  phi_full <- .trim_zeros(phi_full)
  theta <- .trim_zeros(theta)
  if (length(phi_full) == 0L && length(theta) == 0L) {
    return(c(1, rep(0, n - 1L)))
  }
  c(1, stats::ARMAtoMA(ar = phi_full, ma = theta, lag.max = n - 1L))
}

# ---- Coefficients, residuals, information criteria ---------------------------

.free <- function(fit) {
  mask <- fit$mask
  if (is.null(mask)) mask <- rep(TRUE, length(fit$coef))
  as.logical(mask)
}

.free_coef <- function(fit) fit$coef[.free(fit)]

.resid <- function(fit) {
  r <- as.numeric(stats::residuals(fit))
  r[is.finite(r)]
}

# Residuals for diagnostics that rely on consecutive lags (ACF, PACF,
# Ljung-Box): gaps inside the sample would make non-adjacent residuals look
# adjacent, so they are rejected.
.resid_consecutive <- function(fit) {
  r <- as.numeric(stats::residuals(fit))
  ok <- is.finite(r)
  if (!any(ok)) stop("The model has no finite residuals.", call. = FALSE)
  span <- which(ok)
  if (any(!ok[min(span):max(span)])) {
    stop("The residuals contain missing values inside the sample (the series has ",
         "gaps). Autocorrelation diagnostics need a complete sequence.",
         call. = FALSE)
  }
  r[ok]
}

# Reject Box-Cox fits where a function's results would mix scales.
.reject_lambda <- function(fit, what) {
  if (!is.null(fit$lambda)) {
    stop(what, " does not support models with a Box-Cox transformation ",
         "(lambda): its results would mix the transformed and original scales. ",
         "Transform the series yourself (e.g. log(y)) and fit the model without ",
         "`lambda`.", call. = FALSE)
  }
  invisible(fit)
}

# Was the model estimated by conditional sum of squares only (no likelihood)?
.is_css <- function(fit) {
  is.null(fit$aic) || is.na(fit$aic)
}

# Names in `dots` resolved by partial matching against the formals of `fun`;
# unmatched names are returned unchanged.
.resolve_dot_names <- function(dots, fun) {
  nm <- names(dots)
  if (is.null(nm) || any(!nzchar(nm))) {
    if (length(dots)) stop("Arguments passed through `...` must be named.",
                           call. = FALSE)
    return(character(0))
  }
  fn <- names(formals(fun))
  hit <- pmatch(nm, fn, duplicates.ok = TRUE)
  ifelse(is.na(hit), nm, fn[hit])
}

.ic_values <- function(fit) {
  n <- as.numeric(fit$nobs)
  k <- sum(.free(fit)) + 1
  # CSS fits have no likelihood: criteria are not available.
  ll <- if (.is_css(fit)) NA_real_ else as.numeric(fit$loglik)
  aic <- -2 * ll + 2 * k
  aicc <- if (n - k - 1 > 0) aic + 2 * k * (k + 1) / (n - k - 1) else NA_real_
  list(
    n = as.integer(n), k = as.integer(k), loglik = ll,
    AIC = aic, AICc = aicc,
    BIC = -2 * ll + k * log(n),
    HQIC = -2 * ll + 2 * k * log(log(n))
  )
}

.converged <- function(fit) {
  code <- fit$code
  is.null(code) || isTRUE(code == 0L)
}

# Estimate with forecast::Arima, silencing routine warnings. If the default
# CSS-ML method fails (typically "non-stationary AR part from CSS"), retry with
# exact ML so that the grid or recursion is not interrupted.
.fit_arima <- function(y, order, include.constant = NULL, ...) {
  args <- list(y = y, order = order, ...)
  if (!is.null(include.constant)) args$include.constant <- include.constant
  fallback <- FALSE
  fit <- tryCatch(
    suppressWarnings(do.call(forecast::Arima, args)),
    error = function(e) {
      if (!is.null(args$method)) stop(e)
      fallback <<- TRUE
      args$method <- "ML"
      suppressWarnings(do.call(forecast::Arima, args))
    }
  )
  attr(fit, "ml_fallback") <- fallback
  fit
}

# ---- Time axes ---------------------------------------------------------------

# Calendar dates for quarterly, monthly and annual ts objects (period starts,
# as in the Python version); the index for zoo/xts; decimal time for other ts;
# and observation numbers otherwise.
.time_values <- function(x) {
  if (inherits(x, "zoo")) return(zoo::index(x))
  if (stats::is.ts(x)) {
    f <- stats::frequency(x)
    tt <- as.numeric(stats::time(x))
    if (f %in% c(4, 12)) return(zoo::as.Date(zoo::as.yearmon(tt)))
    if (f == 1 && all(abs(tt - round(tt)) < 1e-8) && min(tt) >= 1000 &&
        max(tt) <= 3000) {
      return(as.Date(sprintf("%d-01-01", as.integer(round(tt)))))
    }
    return(tt)
  }
  seq_along(x)
}

.response <- function(fit) {
  y <- forecast::getResponse(fit)
  if (is.null(y)) {
    stop("The observed series cannot be recovered from `fit`. ",
         "Estimate the model with forecast::Arima().", call. = FALSE)
  }
  y
}

.series_name <- function(fit) {
  nm <- fit$series
  if (is.null(nm) || !nzchar(nm)) "y" else nm
}

# ---- Plotting ----------------------------------------------------------------

.at_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey95")
    )
}

.need_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Install the 'plotly' package to use interactive = TRUE.",
         call. = FALSE)
  }
}

.as_interactive <- function(p, tooltip = "text") {
  .need_plotly()
  plotly::ggplotly(p, tooltip = tooltip)
}

# Stack single ggplot panels: patchwork (static) or plotly::subplot.
.stack_panels <- function(panels, heights, title, interactive) {
  if (interactive) {
    .need_plotly()
    ps <- lapply(panels, function(p) plotly::ggplotly(p, tooltip = "text"))
    out <- plotly::subplot(ps, nrows = length(ps), heights = heights,
                           shareX = FALSE, titleX = TRUE, titleY = TRUE,
                           margin = 0.05)
    return(plotly::layout(out, title = list(text = title), showlegend = FALSE))
  }
  patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15)
      )
    )
}

# ---- Printing of named summaries ---------------------------------------------

.at_summary <- function(x, title) {
  structure(x, class = "arimatools_summary", title = title)
}

#' @export
print.arimatools_summary <- function(x, digits = 6, ...) {
  title <- attr(x, "title")
  if (!is.null(title)) cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  vals <- vapply(unclass(x), function(v) {
    if (is.numeric(v) && !is.integer(v)) format(signif(v, digits)) else format(v)
  }, character(1))
  labs <- format(names(vals))
  cat(paste0(labs, "  ", vals), sep = "\n")
  invisible(x)
}

#' @export
as.data.frame.arimatools_summary <- function(x, ...) {
  data.frame(
    item = names(x),
    value = vapply(unclass(x), function(v) format(v), character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
