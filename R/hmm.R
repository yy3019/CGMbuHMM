#' Fit the beta-uniform HMM to one CGM series
#'
#' The emission template is held fixed while the common segment-start
#' distribution and transition matrix are estimated by the Baum-Welch
#' expectation-maximization algorithm. Gaps longer than max_gap_minutes start
#' new segments, so transitions are not counted across those gaps.
#'
#' @param glucose Numeric glucose measurements.
#' @param time Observation times. Numeric values are interpreted as minutes.
#' @param templates Emission templates returned by fit_emission_templates, or
#'   one cgm_emission_template.
#' @param stratum Stratum identifying which template to use.
#' @param max_gap_minutes Largest gap retained within a continuous segment.
#' @param step_minutes Minutes represented by one Markov transition.
#' @param max_iter Maximum Baum-Welch iterations.
#' @param tolerance Absolute scaled log-likelihood convergence tolerance.
#' @param initial_probabilities Optional initial state distribution.
#' @param transition_matrix Optional initial transition matrix.
#' @param prevalence_weight Apply baseline prevalence weights in overlaps.
#' @param return_posterior Return observation-level posterior probabilities and
#'   decoded states.
#'
#' @return A cgm_hmm_fit object containing the transition matrix, initial
#'   probabilities, MFPT values, convergence information, and optional
#'   observation-level posterior probabilities.
#' @export
#'
#' @examples
#' data <- simulate_cgm_data(n_participants = 4, n_per_period = 96, seed = 2)
#' baseline <- data[data$period == "baseline", ]
#' templates <- fit_emission_templates(
#'   baseline,
#'   stratum = "age_group",
#'   cutoffs = list(adults = c(70, 180), children = c(70, 180)),
#'   min_state_n = 5
#' )
#' one <- data[data$participant_id == "P001" & data$period == "rct", ]
#' fit_hmm_cgm(
#'   one$glucose,
#'   one$timestamp,
#'   templates,
#'   stratum = one$age_group[1]
#' )
fit_hmm_cgm <- function(
    glucose,
    time,
    templates,
    stratum = NULL,
    max_gap_minutes = 6,
    step_minutes = 5,
    max_iter = 50L,
    tolerance = 1e-4,
    initial_probabilities = NULL,
    transition_matrix = NULL,
    prevalence_weight = TRUE,
    return_posterior = FALSE) {
  if (!is.numeric(glucose) || length(glucose) != length(time)) {
    .stopf("glucose must be numeric and have the same length as time.")
  }
  if (!is.numeric(max_gap_minutes) || length(max_gap_minutes) != 1L ||
      !is.finite(max_gap_minutes) || max_gap_minutes < 0) {
    .stopf("max_gap_minutes must be one non-negative finite number.")
  }
  max_iter <- as.integer(max_iter)
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    .stopf("max_iter must be a positive integer.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    .stopf("tolerance must be one positive finite number.")
  }
  template <- .resolve_template(templates, stratum)
  time_minutes <- .time_to_minutes(time)
  keep <- is.finite(glucose) & is.finite(time_minutes) &
    glucose > template$domain[1L] & glucose < template$domain[2L]
  if (!any(keep)) {
    .stopf("No usable observations remain after filtering.")
  }
  original_rows <- which(keep)
  glucose_used <- glucose[keep]
  time_used <- time[keep]
  time_minutes <- time_minutes[keep]
  ord <- order(time_minutes, seq_along(time_minutes))
  glucose_used <- glucose_used[ord]
  time_used <- time_used[ord]
  time_minutes <- time_minutes[ord]
  original_rows <- original_rows[ord]
  segments <- .segment_bounds(time_minutes, max_gap_minutes)

  log_scores <- emission_log_scores(
    glucose_used,
    template,
    prevalence_weight = prevalence_weight
  )
  row_max <- apply(log_scores, 1L, max)
  B <- exp(sweep(log_scores, 1L, row_max, FUN = "-"))
  B <- pmax(B, 1e-300)

  if (is.null(initial_probabilities)) {
    initial_probabilities <- rep(1 / 3, 3)
  }
  if (!is.numeric(initial_probabilities) ||
      length(initial_probabilities) != 3L ||
      any(!is.finite(initial_probabilities)) ||
      any(initial_probabilities < 0) ||
      abs(sum(initial_probabilities) - 1) > 1e-8) {
    .stopf("initial_probabilities must be three non-negative values summing to 1.")
  }
  if (is.null(transition_matrix)) {
    transition_matrix <- matrix(1 / 3, nrow = 3, ncol = 3)
  }
  transition_matrix <- .validate_probability_matrix(transition_matrix)

  fitted <- baum_welch_cpp(
    B,
    segments$starts,
    segments$ends,
    as.numeric(initial_probabilities),
    transition_matrix,
    max_iter,
    tolerance,
    isTRUE(return_posterior)
  )
  dimnames(fitted$transition_matrix) <- list(.state_names, .state_names)
  names(fitted$initial_probabilities) <- .state_names
  names(fitted$posterior_state_counts) <- .state_names
  passage <- mfpt(
    fitted$transition_matrix,
    target = "middle",
    step_minutes = step_minutes
  )

  posterior_data <- NULL
  if (isTRUE(return_posterior)) {
    colnames(fitted$posterior) <- paste0("posterior_", .state_names)
    posterior_data <- data.frame(
      original_row = original_rows,
      time = time_used,
      glucose = glucose_used,
      observed_state = assign_cgm_states(
        glucose_used,
        template$cutoffs,
        labels = TRUE
      ),
      decoded_state = ordered(
        max.col(fitted$posterior, ties.method = "first"),
        levels = 1:3,
        labels = .state_names
      ),
      fitted$posterior,
      check.names = FALSE
    )
  }

  structure(
    list(
      transition_matrix = fitted$transition_matrix,
      initial_probabilities = fitted$initial_probabilities,
      mfpt = passage,
      posterior_state_counts = fitted$posterior_state_counts,
      posterior = posterior_data,
      n_obs = length(glucose_used),
      n_segments = length(segments$starts),
      loglik_scaled = fitted$loglik_scaled,
      iterations = fitted$iterations,
      converged = fitted$converged,
      stratum = if (is.null(stratum)) ".all" else as.character(stratum),
      cutoffs = template$cutoffs,
      settings = list(
        max_gap_minutes = max_gap_minutes,
        step_minutes = step_minutes,
        max_iter = max_iter,
        tolerance = tolerance,
        prevalence_weight = prevalence_weight
      ),
      call = match.call()
    ),
    class = "cgm_hmm_fit"
  )
}
