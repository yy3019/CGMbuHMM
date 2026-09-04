#' CGMmfpt: transition and HMM methods for CGM data
#'
#' CGMmfpt estimates participant-level empirical or hidden Markov transition
#' matrices and summarizes them using mean first-passage time (MFPT). The HMM
#' uses beta-uniform mixture emission templates with optional overlap around
#' adjacent state boundaries.
#'
#' @useDynLib CGMmfpt, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head tail
#' @keywords internal
"_PACKAGE"
