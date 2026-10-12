# Parameter inference and ARIMA root properties ------------------------------

#' AR and MA root diagnostics
#'
#' Roots of the AR and MA polynomials of the stationary ARMA component of a
#' fitted ARIMA model, with their moduli, inverse roots and inverse moduli.
#'
#' Seasonal models use the expanded (multiplied) polynomials. AR roots outside
#' the unit circle imply causality; MA roots outside the unit circle imply
#' invertibility. Equivalently, the inverse roots lie inside the unit circle.
#'
#' Related: `forecast::autoplot(fit)` plots the same inverse roots.
#'
#' @param fit A fitted model from [forecast::Arima()], [forecast::auto.arima()]
#'   or [stats::arima()].
#' @return A data frame with columns `type`, `root`, `value` (complex),
#'   `modulus`, `inverse_root` (complex) and `inverse_modulus`.
#' @examples
#' set.seed(12345)
#' y <- arima.sim(list(ar = 0.7, ma = 0.4), n = 240)
#' fit <- forecast::Arima(y, order = c(1, 0, 1), include.mean = FALSE)
#' root_diagnostics(fit)
#' @export
root_diagnostics <- function(fit) {
  .check_fit(fit)
  make <- function(type, roots) {
    if (length(roots) == 0L) return(NULL)
    data.frame(
      type = type,
      root = seq_along(roots),
      value = roots,
      modulus = Mod(roots),
      inverse_root = 1 / roots,
      inverse_modulus = Mod(1 / roots),
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(make("AR", .ar_roots(fit)), make("MA", .ma_roots(fit)))
  if (is.null(out)) {
    out <- data.frame(type = character(0), root = integer(0),
                      value = complex(0), modulus = numeric(0),
                      inverse_root = complex(0), inverse_modulus = numeric(0))
  }
  out
}

#' Is the AR component causal?
#'
#' `TRUE` when every AR root lies outside the unit circle. For ARIMA(p,d,q)
#' models this refers to the AR polynomial of the stationary ARMA component,
#' not to the stationarity of the original series when d > 0.
#'
#' @inheritParams root_diagnostics
#' @param tol Numerical tolerance when comparing root moduli with one.
#' @return A single logical value.
#' @export
is_causal <- function(fit, tol = 1e-10) {
  .check_fit(fit)
  r <- .ar_roots(fit)
  length(r) == 0L || all(Mod(r) > 1 + tol)
}

#' Is the MA component invertible?
#'
#' `TRUE` when every MA root lies outside the unit circle.
#'
#' @inheritParams is_causal
#' @return A single logical value.
#' @export
is_invertible <- function(fit, tol = 1e-10) {
  .check_fit(fit)
  r <- .ma_roots(fit)
  length(r) == 0L || all(Mod(r) > 1 + tol)
}

#' Common and near-common AR/MA roots
#'
#' Pairs of AR and MA inverse roots whose distance in the complex plane is at
#' most `tol`. Pairs within `exact_tol` are classified as "Common"; the rest as
#' "Near-common". Such pairs may indicate cancellation between the AR and MA
#' polynomials and possible overparameterization.
#'
#' @inheritParams root_diagnostics
#' @param tol Maximum distance between inverse roots to report a pair.
#' @param exact_tol Distance at or below which a pair is "Common".
#' @return A data frame with columns `ar_root`, `ma_root`, `ar_inverse_root`,
#'   `ma_inverse_root`, `distance` and `classification` (possibly zero rows).
#' @export
common_roots <- function(fit, tol = 0.05, exact_tol = 1e-10) {
  .check_fit(fit)
  ar_inv <- 1 / .ar_roots(fit)
  ma_inv <- 1 / .ma_roots(fit)
  empty <- data.frame(
    ar_root = integer(0), ma_root = integer(0),
    ar_inverse_root = complex(0), ma_inverse_root = complex(0),
    distance = numeric(0), classification = character(0),
    stringsAsFactors = FALSE
  )
  if (length(ar_inv) == 0L || length(ma_inv) == 0L) return(empty)
  grid <- expand.grid(ma_root = seq_along(ma_inv), ar_root = seq_along(ar_inv))
  grid <- grid[, c("ar_root", "ma_root")]
  grid$ar_inverse_root <- ar_inv[grid$ar_root]
  grid$ma_inverse_root <- ma_inv[grid$ma_root]
  grid$distance <- Mod(grid$ar_inverse_root - grid$ma_inverse_root)
  grid <- grid[grid$distance <= tol, , drop = FALSE]
  if (nrow(grid) == 0L) return(empty)
  grid$classification <- ifelse(grid$distance <= exact_tol, "Common", "Near-common")
  rownames(grid) <- NULL
  grid
}

#' Plot inverse AR and MA roots
#'
#' Inverse roots against the unit circle: AR as squares, MA as circles, and
#' common or near-common pairs marked with a cross at the midpoint of the pair.
#' The axis range expands automatically if a root lies outside the circle.
#'
#' Related: `forecast::autoplot(fit)` draws the inverse roots in separate AR and
#' MA panels, without the common-root markers.
#'
#' @inheritParams common_roots
#' @param title Plot title.
#' @param interactive If `TRUE`, return a plotly widget (requires 'plotly').
#' @return A ggplot object (or a plotly widget when `interactive = TRUE`).
#' @export
plot_inverse_roots <- function(fit, tol = 0.05, exact_tol = 1e-10,
                               title = "Inverse Roots", interactive = FALSE) {
  .check_fit(fit)
  ar_inv <- 1 / .ar_roots(fit)
  ma_inv <- 1 / .ma_roots(fit)

  theta <- seq(0, 2 * pi, length.out = 500)
  circle <- data.frame(x = cos(theta), y = sin(theta))

  pts <- rbind(
    if (length(ar_inv)) data.frame(x = Re(ar_inv), y = Im(ar_inv),
                                   kind = "AR inverse roots", mod = Mod(ar_inv)),
    if (length(ma_inv)) data.frame(x = Re(ma_inv), y = Im(ma_inv),
                                   kind = "MA inverse roots", mod = Mod(ma_inv))
  )
  pairs <- common_roots(fit, tol = tol, exact_tol = exact_tol)

  lim <- max(1.1, 1.05 * c(Mod(ar_inv), Mod(ma_inv)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_path(data = circle, ggplot2::aes(x = .data$x, y = .data$y),
                       colour = "grey40")

  if (!is.null(pts)) {
    pts$text <- sprintf("%s<br>Real: %.4f<br>Imaginary: %.4f<br>Modulus: %.4f",
                        sub("s$", "", pts$kind), pts$x, pts$y, pts$mod)
    p <- p + suppressWarnings(ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data$x, y = .data$y, shape = .data$kind,
                   colour = .data$kind, text = .data$text),
      size = 3.5
    ))
  }
  if (nrow(pairs)) {
    mid <- (pairs$ar_inverse_root + pairs$ma_inverse_root) / 2
    cr <- data.frame(
      x = Re(mid), y = Im(mid), kind = "Common / near-common",
      text = sprintf("%s roots<br>Real: %.4f<br>Imaginary: %.4f<br>Distance: %.4f",
                     pairs$classification, Re(mid), Im(mid), pairs$distance)
    )
    p <- p + suppressWarnings(ggplot2::geom_point(
      data = cr,
      ggplot2::aes(x = .data$x, y = .data$y, shape = .data$kind,
                   colour = .data$kind, text = .data$text),
      size = 5, stroke = 1.3
    ))
  }

  p <- p +
    ggplot2::scale_shape_manual(
      values = c("AR inverse roots" = 15, "MA inverse roots" = 16,
                 "Common / near-common" = 4), name = NULL) +
    ggplot2::scale_colour_manual(
      values = c("AR inverse roots" = "#1F5A96", "MA inverse roots" = "#D55E00",
                 "Common / near-common" = "black"), name = NULL) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(title = title, x = "Real", y = "Imaginary") +
    .at_theme()

  if (interactive) .as_interactive(p) else p
}

#' Summary of root properties
#'
#' Causality, invertibility, and counts of common and near-common AR/MA root
#' pairs, with the minimum AR-MA inverse-root distance among reported pairs.
#'
#' @inheritParams common_roots
#' @return A named list of class `arimatools_summary` (prints as a two-column
#'   display; use `as.data.frame()` for a table).
#' @export
root_summary <- function(fit, tol = 0.05, exact_tol = 1e-10) {
  .check_fit(fit)
  pairs <- common_roots(fit, tol = tol, exact_tol = exact_tol)
  .at_summary(
    list(
      "Causal AR component" = is_causal(fit),
      "Invertible MA component" = is_invertible(fit),
      "Common root pairs" = sum(pairs$classification == "Common"),
      "Near-common root pairs" = sum(pairs$classification == "Near-common"),
      "Minimum AR-MA distance" = if (nrow(pairs)) min(pairs$distance) else NA_real_
    ),
    "Root diagnostics"
  )
}

#' Coefficient table with asymptotic inference
#'
#' Estimates, standard errors, z-statistics, two-sided asymptotic p-values and
#' normal confidence intervals for the freely estimated coefficients.
#'
#' No significance stars or automatic classification are provided: in dynamic
#' models, coefficient relevance should not be judged mechanically by a fixed
#' p-value threshold.
#'
#' Related: `lmtest::coeftest(fit)` and `confint(fit)` give the same numbers.
#'
#' @inheritParams root_diagnostics
#' @param confidence Confidence level for the intervals.
#' @return A data frame with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `conf.low` and `conf.high` (broom conventions).
#' @export
coefficient_diagnostics <- function(fit, confidence = 0.95) {
  .check_fit(fit)
  .check_level(confidence)
  est <- .free_coef(fit)
  se <- sqrt(diag(fit$var.coef))
  z <- est / se
  crit <- stats::qnorm(1 - (1 - confidence) / 2)
  out <- data.frame(
    term = names(est),
    estimate = unname(est),
    std.error = unname(se),
    statistic = unname(z),
    p.value = 2 * stats::pnorm(-abs(unname(z))),
    conf.low = unname(est - crit * se),
    conf.high = unname(est + crit * se),
    stringsAsFactors = FALSE
  )
  attr(out, "confidence") <- confidence
  out
}

#' Correlation matrix of parameter estimates
#'
#' Large correlations may indicate that individual coefficients are
#' imprecisely identified or that the model is overparameterized.
#'
#' Unlike Statsmodels, R does not report a covariance for the innovation
#' variance, so `sigma2` is never included. Equivalent to
#' `cov2cor(vcov(fit))`.
#'
#' @inheritParams root_diagnostics
#' @return A numeric correlation matrix with coefficient names as dimnames.
#' @export
parameter_correlation <- function(fit) {
  .check_fit(fit)
  stats::cov2cor(fit$var.coef)
}

#' Wald test that a group of coefficients is jointly zero
#'
#' @inheritParams root_diagnostics
#' @param parameters Character vector of coefficient names (as in
#'   `names(coef(fit))`, e.g. `c("ar2", "ma2")`).
#' @return A named list of class `arimatools_summary` with the Wald statistic,
#'   degrees of freedom and asymptotic chi-squared p-value.
#' @export
joint_significance <- function(fit, parameters) {
  .check_fit(fit)
  parameters <- as.character(parameters)
  est <- .free_coef(fit)
  if (length(parameters) == 0L) {
    stop("At least one parameter must be specified.", call. = FALSE)
  }
  unknown <- setdiff(parameters, names(est))
  if (length(unknown)) {
    stop("Unknown parameter(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(est), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(parameters)) {
    stop("Parameter names must not be repeated.", call. = FALSE)
  }
  b <- est[parameters]
  V <- fit$var.coef[parameters, parameters, drop = FALSE]
  w <- as.numeric(crossprod(b, solve(V, b)))
  df <- length(parameters)
  .at_summary(
    list("Wald statistic" = w, "Degrees of freedom" = df,
         "P-value" = stats::pchisq(w, df, lower.tail = FALSE)),
    "Joint significance test"
  )
}
