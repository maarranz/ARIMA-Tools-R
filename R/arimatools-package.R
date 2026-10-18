#' arimatools: teaching-oriented tools for empirical ARIMA modelling
#'
#' R port of the Python package ARIMA-Tools (ATSE, UC3M). Every function takes
#' a model fitted with [forecast::Arima()], [forecast::auto.arima()] or
#' [stats::arima()], except the functions that estimate models themselves
#' ([order_search()], [recursive_estimation()], [auto_arima()]).
#'
#' Part of the ATSE software suite, with ARIMA-Tools (Python), ARMA-Theory,
#' ADL-Theory, Forecasting-Tools and ADL-Tools (forthcoming).
#'
#' Function names match the Python API. Tables follow R/broom conventions
#' (`term`, `estimate`, `std.error`, `p.value`, ...). Figures are ggplot2
#' objects; pass `interactive = TRUE` for a plotly widget.
#'
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"
