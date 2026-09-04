#' @export
print.cgm_emission_templates <- function(x, ...) {
  cat("Pooled beta-uniform CGM emission templates\n")
  cat("  strata: ", paste(names(x$templates), collapse = ", "), "\n", sep = "")
  cat("  overlap half-width: ", x$settings$overlap, " mg/dL\n", sep = "")
  invisible(x)
}

#' @export
print.cgm_hmm_fit <- function(x, ...) {
  cat("Three-state CGM hidden Markov model\n")
  cat("  observations: ", x$n_obs, " in ", x$n_segments,
      " segment(s)\n", sep = "")
  cat("  EM iterations: ", x$iterations,
      if (isTRUE(x$converged)) " (converged)\n" else " (not converged)\n",
      sep = "")
  cat("  MFPT to middle state (minutes):\n")
  print(x$mfpt[c("lower", "upper")])
  invisible(x)
}

#' @export
print.cgm_empirical_fit <- function(x, ...) {
  cat("Three-state empirical CGM transition model\n")
  cat("  observations: ", x$n_obs, " in ", x$n_segments,
      " segment(s)\n", sep = "")
  cat("  observed transitions: ", x$n_transitions, "\n", sep = "")
  cat("  MFPT to middle state (minutes):\n")
  print(x$mfpt[c("lower", "upper")])
  invisible(x)
}

#' @export
print.cgm_mfpt_results <- function(x, ...) {
  method <- attr(x, "method")
  cat("CGM MFPT estimates", if (!is.null(method)) paste0(" (", method, ")"),
      "\n", sep = "")
  NextMethod("print")
  invisible(x)
}
