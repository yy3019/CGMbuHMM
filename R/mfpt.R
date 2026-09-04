#' Mean first-passage time to a target state
#'
#' Calculates expected time to first reach a target state from every state in a
#' finite, discrete-time Markov chain. An infinite value is returned when the
#' target is not reached with probability one.
#'
#' @param transition_matrix Row-stochastic transition matrix.
#' @param target Index or name of the target state.
#' @param step_minutes Minutes represented by one transition.
#' @param tolerance Numerical tolerance used to assess hitting probability.
#'
#' @return A named numeric vector in minutes. The target state's value is zero.
#' @export
#'
#' @examples
#' A <- matrix(c(
#'   0.7, 0.3, 0,
#'   0.1, 0.8, 0.1,
#'   0, 0.4, 0.6
#' ), nrow = 3, byrow = TRUE)
#' dimnames(A) <- list(c("lower", "middle", "upper"),
#'                     c("lower", "middle", "upper"))
#' mfpt(A, target = "middle", step_minutes = 5)
mfpt <- function(transition_matrix, target = 2L, step_minutes = 5,
                 tolerance = 1e-8) {
  A <- .validate_probability_matrix(transition_matrix)
  if (!is.numeric(step_minutes) || length(step_minutes) != 1L ||
      !is.finite(step_minutes) || step_minutes <= 0) {
    .stopf("step_minutes must be one positive finite number.")
  }
  state_names <- rownames(A)
  if (is.character(target)) {
    if (is.null(state_names) || length(target) != 1L ||
        !target %in% state_names) {
      .stopf("Character target must match a transition-matrix row name.")
    }
    target <- match(target, state_names)
  }
  target <- as.integer(target)
  if (length(target) != 1L || is.na(target) ||
      target < 1L || target > nrow(A)) {
    .stopf("target must identify one state in transition_matrix.")
  }

  other <- setdiff(seq_len(nrow(A)), target)
  result <- rep(Inf, nrow(A))
  result[target] <- 0

  can_reach <- rep(FALSE, nrow(A))
  can_reach[target] <- TRUE
  repeat {
    previous <- can_reach
    can_reach <- can_reach | rowSums(A[, can_reach, drop = FALSE]) > tolerance
    if (identical(can_reach, previous)) {
      break
    }
  }
  reachable <- other[can_reach[other]]
  if (length(reachable)) {
    Q_reachable <- A[reachable, reachable, drop = FALSE]
    to_target <- A[reachable, target]
    hitting <- tryCatch(
      as.numeric(solve(diag(length(reachable)) - Q_reachable, to_target)),
      error = function(e) rep(NA_real_, length(reachable))
    )
    certain <- reachable[is.finite(hitting) & hitting >= 1 - tolerance]
    if (length(certain)) {
      passage <- tryCatch(
        as.numeric(solve(
          diag(length(certain)) - A[certain, certain, drop = FALSE],
          rep(1, length(certain))
        )),
        error = function(e) rep(NA_real_, length(certain))
      )
      finite <- is.finite(passage) & passage >= 0
      result[certain[finite]] <- passage[finite] * step_minutes
    }
  }
  if (!is.null(state_names)) {
    names(result) <- state_names
  } else {
    names(result) <- paste0("state", seq_len(nrow(A)))
  }
  result
}
