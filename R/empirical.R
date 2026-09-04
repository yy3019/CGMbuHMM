#' Fit an empirical three-state CGM transition model
#'
#' Glucose measurements are assigned to lower, middle, and upper states using
#' hard cutoffs. The transition matrix is estimated from consecutive
#' observations within continuous segments. Sparse transition counts can be
#' stabilized by adding a user-supplied reference transition matrix as
#' pseudocounts.
#'
#' @param glucose Numeric glucose measurements.
#' @param time Observation times. Numeric values are interpreted as minutes.
#' @param cutoffs Two increasing state boundaries. The clinical default is
#'   70 and 180 mg/dL.
#' @param domain Admissible glucose domain. Values at or outside its endpoints
#'   are excluded.
#' @param reference_transition Optional three by three row-stochastic reference
#'   transition matrix.
#' @param prior_strength Non-negative pseudocount total added to each row of
#'   the observed transition counts. Requires reference_transition when
#'   positive.
#' @param max_gap_minutes Largest gap retained within a continuous segment.
#' @param step_minutes Minutes represented by one Markov transition.
#' @param zero_row How to handle a row with no observed or prior transitions.
#' @param return_states Return the ordered observations and assigned states.
#'
#' @return A cgm_empirical_fit object containing counts, the transition matrix,
#'   and MFPT values.
#' @export
#'
#' @examples
#' fit_empirical_cgm(
#'   glucose = c(65, 72, 110, 190, 175, 140),
#'   time = seq(0, 25, by = 5)
#' )
fit_empirical_cgm <- function(
    glucose,
    time,
    cutoffs = c(70, 180),
    domain = c(40, 401),
    reference_transition = NULL,
    prior_strength = 0,
    max_gap_minutes = 6,
    step_minutes = 5,
    zero_row = c("self", "uniform", "error"),
    return_states = FALSE) {
  zero_row <- match.arg(zero_row)
  domain <- .validate_domain(domain)
  cutoffs <- .validate_cutoff_vector(cutoffs, domain)
  if (!is.numeric(glucose) || length(glucose) != length(time)) {
    .stopf("glucose must be numeric and have the same length as time.")
  }
  if (!is.numeric(prior_strength) || length(prior_strength) != 1L ||
      !is.finite(prior_strength) || prior_strength < 0) {
    .stopf("prior_strength must be one non-negative finite number.")
  }
  if (!is.numeric(max_gap_minutes) || length(max_gap_minutes) != 1L ||
      !is.finite(max_gap_minutes) || max_gap_minutes < 0) {
    .stopf("max_gap_minutes must be one non-negative finite number.")
  }
  if (prior_strength > 0 && is.null(reference_transition)) {
    .stopf("reference_transition is required when prior_strength is positive.")
  }
  if (!is.null(reference_transition)) {
    reference_transition <- .validate_probability_matrix(reference_transition)
    if (!identical(dim(reference_transition), c(3L, 3L))) {
      .stopf("reference_transition must be a three by three matrix.")
    }
  }

  time_minutes <- .time_to_minutes(time)
  keep <- is.finite(glucose) & is.finite(time_minutes) &
    glucose > domain[1L] & glucose < domain[2L]
  if (!any(keep)) {
    .stopf("No usable observations remain after filtering.")
  }
  original_rows <- which(keep)
  glucose_used <- glucose[keep]
  time_used <- time[keep]
  time_minutes <- time_minutes[keep]
  ord <- order(time_minutes, seq_along(time_minutes))
  original_rows <- original_rows[ord]
  glucose_used <- glucose_used[ord]
  time_used <- time_used[ord]
  time_minutes <- time_minutes[ord]
  states <- assign_cgm_states(glucose_used, cutoffs)

  counts <- transition_counts(
    states,
    time_minutes,
    max_gap_minutes = max_gap_minutes,
    n_states = 3L
  )
  smoothed_counts <- counts
  if (prior_strength > 0) {
    smoothed_counts <- smoothed_counts + prior_strength * reference_transition
  }
  transition_matrix <- .normalize_rows(smoothed_counts, zero_row = zero_row)
  dimnames(transition_matrix) <- list(.state_names, .state_names)
  dimnames(smoothed_counts) <- list(.state_names, .state_names)
  segments <- .segment_bounds(time_minutes, max_gap_minutes)

  state_data <- NULL
  if (isTRUE(return_states)) {
    state_data <- data.frame(
      original_row = original_rows,
      time = time_used,
      glucose = glucose_used,
      state = ordered(states, levels = 1:3, labels = .state_names),
      check.names = FALSE
    )
  }

  structure(
    list(
      transition_matrix = transition_matrix,
      transition_counts = counts,
      smoothed_counts = smoothed_counts,
      mfpt = mfpt(
        transition_matrix,
        target = "middle",
        step_minutes = step_minutes
      ),
      states = state_data,
      n_obs = length(glucose_used),
      n_segments = length(segments$starts),
      n_transitions = sum(counts),
      cutoffs = cutoffs,
      domain = domain,
      settings = list(
        prior_strength = prior_strength,
        max_gap_minutes = max_gap_minutes,
        step_minutes = step_minutes,
        zero_row = zero_row
      ),
      call = match.call()
    ),
    class = "cgm_empirical_fit"
  )
}
