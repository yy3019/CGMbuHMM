#' Assign CGM measurements to three ordered states
#'
#' The first state contains values less than or equal to the lower cutoff, the
#' second contains values above the lower cutoff and less than or equal to the
#' upper cutoff, and the third contains values above the upper cutoff.
#'
#' @param glucose Numeric glucose measurements.
#' @param cutoffs Two increasing state boundaries. The clinical default is
#'   70 and 180 mg/dL.
#' @param labels If TRUE, return an ordered factor with levels lower, middle,
#'   and upper. Otherwise return integer state indices.
#'
#' @return An integer vector or ordered factor with missing values retained.
#' @export
#'
#' @examples
#' assign_cgm_states(c(55, 70, 120, 180, 240))
assign_cgm_states <- function(glucose, cutoffs = c(70, 180), labels = FALSE) {
  cutoffs <- .validate_cutoff_vector(cutoffs)
  if (!is.numeric(glucose)) {
    .stopf("glucose must be numeric.")
  }
  states <- rep.int(NA_integer_, length(glucose))
  finite <- is.finite(glucose)
  states[finite & glucose <= cutoffs[1L]] <- 1L
  states[finite & glucose > cutoffs[1L] & glucose <= cutoffs[2L]] <- 2L
  states[finite & glucose > cutoffs[2L]] <- 3L
  if (labels) {
    return(ordered(states, levels = 1:3, labels = .state_names))
  }
  states
}

#' Count observed transitions between CGM states
#'
#' Observations are ordered by time. A transition is counted only when the gap
#' between consecutive usable observations does not exceed max_gap_minutes.
#'
#' @param states Integer state labels from 1 through n_states.
#' @param time Observation times. Numeric values are interpreted as minutes.
#' @param max_gap_minutes Largest gap that is treated as continuous.
#' @param n_states Number of states.
#'
#' @return An n_states by n_states transition-count matrix.
#' @export
#'
#' @examples
#' transition_counts(c(1, 2, 2, 3), c(0, 5, 10, 15))
transition_counts <- function(states, time, max_gap_minutes = 6,
                              n_states = 3L) {
  if (length(states) != length(time)) {
    .stopf("states and time must have the same length.")
  }
  if (!is.numeric(max_gap_minutes) || length(max_gap_minutes) != 1L ||
      !is.finite(max_gap_minutes) || max_gap_minutes < 0) {
    .stopf("max_gap_minutes must be one non-negative finite number.")
  }
  n_states <- as.integer(n_states)
  if (length(n_states) != 1L || is.na(n_states) || n_states < 2L) {
    .stopf("n_states must be an integer of at least 2.")
  }
  time_minutes <- .time_to_minutes(time)
  keep <- !is.na(states) & is.finite(time_minutes)
  states <- as.integer(states[keep])
  time_minutes <- time_minutes[keep]
  if (any(states < 1L | states > n_states)) {
    .stopf("states must be integers between 1 and n_states.")
  }
  ord <- order(time_minutes, seq_along(time_minutes))
  states <- states[ord]
  time_minutes <- time_minutes[ord]
  counts <- matrix(
    0,
    nrow = n_states,
    ncol = n_states,
    dimnames = {
      labels <- if (n_states == 3L) {
        .state_names
      } else {
        paste0("state", seq_len(n_states))
      }
      list(labels, labels)
    }
  )
  if (length(states) < 2L) {
    return(counts)
  }
  gaps <- diff(time_minutes)
  valid <- is.finite(gaps) & gaps >= 0 & gaps <= max_gap_minutes
  if (any(valid)) {
    tab <- table(
      factor(head(states, -1L)[valid], levels = seq_len(n_states)),
      factor(tail(states, -1L)[valid], levels = seq_len(n_states))
    )
    counts[,] <- unclass(tab)
  }
  counts
}
