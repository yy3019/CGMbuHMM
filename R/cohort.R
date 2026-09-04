.result_identity <- function(rows, id, group, stratum) {
  columns <- unique(c(id, group, stratum))
  values <- lapply(columns, function(column) .one_value(rows[[column]], column))
  names(values) <- columns
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

.handle_fit_error <- function(error, on_error, label) {
  message <- conditionMessage(error)
  if (on_error == "stop") {
    .stopf("%s: %s", label, message)
  }
  if (on_error == "warn") {
    warning(sprintf("%s: %s", label, message), call. = FALSE)
  }
  message
}

#' Estimate HMM MFPTs for participant-level CGM series
#'
#' Splits a data frame into participant analysis series, selects the appropriate
#' pooled emission template, and fits one HMM to each series.
#'
#' @param data Data frame containing CGM observations.
#' @param templates Emission templates returned by fit_emission_templates.
#' @param id Name of the participant identifier column.
#' @param time Name of the observation-time column.
#' @param glucose Name of the numeric glucose column.
#' @param group Optional columns that distinguish repeated analysis series,
#'   such as study period and treatment.
#' @param stratum Optional column used to select an emission template.
#' @param min_obs Minimum usable observations required in an analysis series.
#' @param keep_fits Include fitted model objects in a list column.
#' @param on_error Whether a failed series produces a warning, stops the
#'   analysis, or quietly returns a row with missing estimates.
#' @param ... Additional arguments passed to fit_hmm_cgm.
#'
#' @return A cgm_mfpt_results data frame with one row per analysis series.
#' @export
#'
#' @examples
#' data <- simulate_cgm_data(n_participants = 6, n_per_period = 96, seed = 4)
#' templates <- fit_emission_templates(
#'   data[data$period == "baseline", ],
#'   stratum = "age_group",
#'   min_state_n = 5
#' )
#' estimates <- estimate_hmm_mfpt(
#'   data[data$participant_id %in% c("P001", "P002"), ],
#'   templates,
#'   id = "participant_id",
#'   time = "timestamp",
#'   glucose = "glucose",
#'   group = "period",
#'   stratum = "age_group",
#'   min_obs = 40
#' )
#' estimates
estimate_hmm_mfpt <- function(
    data,
    templates,
    id,
    time,
    glucose,
    group = NULL,
    stratum = NULL,
    min_obs = 50L,
    keep_fits = FALSE,
    on_error = c("warn", "stop", "na"),
    ...) {
  on_error <- match.arg(on_error)
  .assert_scalar_character(id, "id")
  .assert_scalar_character(time, "time")
  .assert_scalar_character(glucose, "glucose")
  .assert_scalar_character(stratum, "stratum", allow_null = TRUE)
  if (!is.null(group) && (!is.character(group) || anyNA(group) ||
      any(!nzchar(group)))) {
    .stopf("group must be NULL or a character vector of column names.")
  }
  .assert_columns(data, unique(c(id, time, glucose, group, stratum)))
  min_obs <- as.integer(min_obs)
  if (length(min_obs) != 1L || is.na(min_obs) || min_obs < 2L) {
    .stopf("min_obs must be an integer of at least 2.")
  }

  split_columns <- unique(c(id, group))
  groups <- .split_indices(data, split_columns)
  output <- vector("list", length(groups))
  fits <- vector("list", length(groups))
  group_names <- names(groups)

  for (index in seq_along(groups)) {
    rows <- data[groups[[index]], , drop = FALSE]
    identity <- .result_identity(rows, id, group, stratum)
    template_key <- if (is.null(stratum)) NULL else .one_value(
      as.character(rows[[stratum]]), stratum
    )
    time_minutes <- .time_to_minutes(rows[[time]])
    usable <- is.finite(rows[[glucose]]) & is.finite(time_minutes)
    fit <- NULL
    error_message <- NA_character_
    if (sum(usable) < min_obs) {
      error_message <- sprintf(
        "Only %d usable observations; min_obs is %d.",
        sum(usable), min_obs
      )
      if (on_error == "stop") {
        .stopf("Series '%s': %s", group_names[index], error_message)
      }
      if (on_error == "warn") {
        warning(
          sprintf("Series '%s': %s", group_names[index], error_message),
          call. = FALSE
        )
      }
    } else {
      fit <- tryCatch(
        fit_hmm_cgm(
          rows[[glucose]],
          rows[[time]],
          templates,
          stratum = template_key,
          ...
        ),
        error = function(error) {
          error_message <<- .handle_fit_error(
            error,
            on_error,
            sprintf("Series '%s'", group_names[index])
          )
          NULL
        }
      )
    }
    fits[[index]] <- fit
    metrics <- if (is.null(fit)) {
      data.frame(
        n_obs = sum(usable),
        n_segments = NA_integer_,
        mfpt_lower_to_middle = NA_real_,
        mfpt_upper_to_middle = NA_real_,
        converged = NA,
        iterations = NA_integer_,
        status = error_message,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        n_obs = fit$n_obs,
        n_segments = fit$n_segments,
        mfpt_lower_to_middle = unname(fit$mfpt["lower"]),
        mfpt_upper_to_middle = unname(fit$mfpt["upper"]),
        converged = fit$converged,
        iterations = fit$iterations,
        status = if (isTRUE(fit$converged)) "ok" else "maximum iterations reached",
        stringsAsFactors = FALSE
      )
    }
    output[[index]] <- cbind(identity, metrics)
  }

  result <- do.call(rbind, output)
  rownames(result) <- NULL
  if (isTRUE(keep_fits)) {
    result$fit <- I(fits)
  }
  class(result) <- c("cgm_mfpt_results", class(result))
  attr(result, "method") <- "HMM"
  result
}

.pooled_reference_transitions <- function(data, series_groups,
                                          reference_by, stratum, time, glucose,
                                          cutoffs, domain, max_gap_minutes) {
  counts_by_series <- vector("list", length(series_groups))
  reference_keys <- character(length(series_groups))
  for (index in seq_along(series_groups)) {
    rows <- data[series_groups[[index]], , drop = FALSE]
    stratum_value <- if (is.list(cutoffs)) {
      if (is.null(stratum)) {
        .stopf("stratum is required when cutoffs is a named list.")
      }
      .one_value(as.character(rows[[stratum]]), stratum)
    } else {
      NULL
    }
    series_cutoffs <- .resolve_cutoffs(cutoffs, stratum_value, domain)
    keep <- is.finite(rows[[glucose]]) &
      is.finite(.time_to_minutes(rows[[time]])) &
      rows[[glucose]] > domain[1L] & rows[[glucose]] < domain[2L]
    states <- assign_cgm_states(rows[[glucose]][keep], series_cutoffs)
    counts_by_series[[index]] <- transition_counts(
      states,
      rows[[time]][keep],
      max_gap_minutes = max_gap_minutes
    )
    reference_keys[index] <- .group_key(rows[1L, , drop = FALSE], reference_by)
  }

  references <- lapply(
    split(seq_along(counts_by_series), reference_keys),
    function(indices) {
      pooled <- Reduce(`+`, counts_by_series[indices])
      .normalize_rows(pooled, zero_row = "uniform")
    }
  )
  list(
    counts_by_series = counts_by_series,
    reference_keys = reference_keys,
    references = references
  )
}

#' Estimate empirical MFPTs for participant-level CGM series
#'
#' Each glucose series is converted to hard states and its empirical transition
#' matrix is estimated. If prior_strength is positive, each transition row is
#' smoothed toward a pooled transition matrix calculated within reference_by
#' strata.
#'
#' @inheritParams estimate_hmm_mfpt
#' @param cutoffs Two increasing cutoffs or a named list of cutoff pairs indexed
#'   by stratum.
#' @param domain Admissible glucose domain. Values at or outside its endpoints
#'   are excluded.
#' @param reference_by Columns defining the pooled smoothing reference. The
#'   default uses stratum and all group columns.
#' @param prior_strength Non-negative pseudocount total added to each transition
#'   row. Set to zero for unsmoothed empirical estimates.
#' @param max_gap_minutes Largest gap retained within a continuous segment.
#' @param step_minutes Minutes represented by one Markov transition.
#' @param zero_row How to handle an unsmoothed row with no transitions.
#'
#' @return A cgm_mfpt_results data frame with one row per analysis series.
#' @export
#'
#' @examples
#' data <- simulate_cgm_data(n_participants = 6, n_per_period = 96, seed = 5)
#' estimate_empirical_mfpt(
#'   data,
#'   id = "participant_id",
#'   time = "timestamp",
#'   glucose = "glucose",
#'   group = c("period", "treatment"),
#'   stratum = "age_group",
#'   reference_by = c("age_group", "period"),
#'   prior_strength = 20,
#'   min_obs = 40
#' )
estimate_empirical_mfpt <- function(
    data,
    id,
    time,
    glucose,
    group = NULL,
    stratum = NULL,
    cutoffs = c(70, 180),
    domain = c(40, 401),
    reference_by = unique(c(stratum, group)),
    prior_strength = 0,
    min_obs = 50L,
    max_gap_minutes = 6,
    step_minutes = 5,
    zero_row = c("self", "uniform", "error"),
    keep_fits = FALSE,
    on_error = c("warn", "stop", "na")) {
  on_error <- match.arg(on_error)
  zero_row <- match.arg(zero_row)
  domain <- .validate_domain(domain)
  .assert_scalar_character(id, "id")
  .assert_scalar_character(time, "time")
  .assert_scalar_character(glucose, "glucose")
  .assert_scalar_character(stratum, "stratum", allow_null = TRUE)
  if (!is.null(group) && (!is.character(group) || anyNA(group) ||
      any(!nzchar(group)))) {
    .stopf("group must be NULL or a character vector of column names.")
  }
  if (is.null(reference_by)) {
    reference_by <- character()
  }
  if (!is.character(reference_by) || anyNA(reference_by) ||
      any(!nzchar(reference_by))) {
    .stopf("reference_by must be a character vector of column names.")
  }
  .assert_columns(
    data,
    unique(c(id, time, glucose, group, stratum, reference_by))
  )
  if (!is.numeric(prior_strength) || length(prior_strength) != 1L ||
      !is.finite(prior_strength) || prior_strength < 0) {
    .stopf("prior_strength must be one non-negative finite number.")
  }
  min_obs <- as.integer(min_obs)
  if (length(min_obs) != 1L || is.na(min_obs) || min_obs < 2L) {
    .stopf("min_obs must be an integer of at least 2.")
  }

  split_columns <- unique(c(id, group))
  groups <- .split_indices(data, split_columns)
  pooled <- NULL
  if (prior_strength > 0) {
    pooled <- .pooled_reference_transitions(
      data,
      groups,
      reference_by,
      stratum,
      time,
      glucose,
      cutoffs,
      domain,
      max_gap_minutes
    )
  }

  output <- vector("list", length(groups))
  fits <- vector("list", length(groups))
  group_names <- names(groups)
  for (index in seq_along(groups)) {
    rows <- data[groups[[index]], , drop = FALSE]
    identity <- .result_identity(rows, id, group, stratum)
    stratum_value <- if (is.null(stratum)) NULL else .one_value(
      as.character(rows[[stratum]]), stratum
    )
    series_cutoffs <- .resolve_cutoffs(cutoffs, stratum_value, domain)
    usable <- is.finite(rows[[glucose]]) &
      is.finite(.time_to_minutes(rows[[time]])) &
      rows[[glucose]] > domain[1L] & rows[[glucose]] < domain[2L]
    fit <- NULL
    error_message <- NA_character_
    if (sum(usable) < min_obs) {
      error_message <- sprintf(
        "Only %d usable observations; min_obs is %d.",
        sum(usable), min_obs
      )
      if (on_error == "stop") {
        .stopf("Series '%s': %s", group_names[index], error_message)
      }
      if (on_error == "warn") {
        warning(
          sprintf("Series '%s': %s", group_names[index], error_message),
          call. = FALSE
        )
      }
    } else {
      reference <- if (prior_strength > 0) {
        pooled$references[[pooled$reference_keys[index]]]
      } else {
        NULL
      }
      fit <- tryCatch(
        fit_empirical_cgm(
          rows[[glucose]],
          rows[[time]],
          cutoffs = series_cutoffs,
          domain = domain,
          reference_transition = reference,
          prior_strength = prior_strength,
          max_gap_minutes = max_gap_minutes,
          step_minutes = step_minutes,
          zero_row = zero_row
        ),
        error = function(error) {
          error_message <<- .handle_fit_error(
            error,
            on_error,
            sprintf("Series '%s'", group_names[index])
          )
          NULL
        }
      )
    }
    fits[[index]] <- fit
    metrics <- if (is.null(fit)) {
      data.frame(
        n_obs = sum(usable),
        n_segments = NA_integer_,
        n_transitions = NA_integer_,
        mfpt_lower_to_middle = NA_real_,
        mfpt_upper_to_middle = NA_real_,
        status = error_message,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        n_obs = fit$n_obs,
        n_segments = fit$n_segments,
        n_transitions = fit$n_transitions,
        mfpt_lower_to_middle = unname(fit$mfpt["lower"]),
        mfpt_upper_to_middle = unname(fit$mfpt["upper"]),
        status = "ok",
        stringsAsFactors = FALSE
      )
    }
    output[[index]] <- cbind(identity, metrics)
  }

  result <- do.call(rbind, output)
  rownames(result) <- NULL
  if (isTRUE(keep_fits)) {
    result$fit <- I(fits)
  }
  class(result) <- c("cgm_mfpt_results", class(result))
  attr(result, "method") <- "empirical"
  attr(result, "prior_strength") <- prior_strength
  result
}

#' Estimate HMM and empirical CGM MFPTs
#'
#' Convenience wrapper that applies both methods to the same participant-level
#' analysis series. Empirical cutoffs default to those stored in templates.
#'
#' @inheritParams estimate_hmm_mfpt
#' @param empirical_cutoffs Cutoffs used by the empirical method. If NULL, use
#'   the cutoffs stored in templates.
#' @param empirical_domain Admissible glucose domain for the empirical method.
#'   If NULL, use the domain stored in templates.
#' @param empirical_reference_by Columns defining empirical smoothing strata.
#' @param empirical_prior_strength Empirical transition pseudocount strength.
#' @param hmm_args Named list of additional arguments for estimate_hmm_mfpt.
#' @param empirical_args Named list of additional arguments for
#'   estimate_empirical_mfpt.
#'
#' @return A list with hmm and empirical cgm_mfpt_results data frames.
#' @export
estimate_cgm_mfpt <- function(
    data,
    templates,
    id,
    time,
    glucose,
    group = NULL,
    stratum = NULL,
    min_obs = 50L,
    empirical_cutoffs = NULL,
    empirical_domain = NULL,
    empirical_reference_by = unique(c(stratum, group)),
    empirical_prior_strength = 0,
    hmm_args = list(),
    empirical_args = list()) {
  if (!is.list(hmm_args) || !is.list(empirical_args)) {
    .stopf("hmm_args and empirical_args must be lists.")
  }
  if (is.null(empirical_cutoffs)) {
    empirical_cutoffs <- .cutoffs_from_templates(templates)
  }
  if (is.null(empirical_domain)) {
    if (!inherits(templates, "cgm_emission_templates")) {
      .stopf("templates must be created by fit_emission_templates().")
    }
    empirical_domain <- templates$settings$domain
  }
  common <- list(
    data = data,
    id = id,
    time = time,
    glucose = glucose,
    group = group,
    stratum = stratum,
    min_obs = min_obs
  )
  list(
    hmm = do.call(
      estimate_hmm_mfpt,
      c(common, list(templates = templates), hmm_args)
    ),
    empirical = do.call(
      estimate_empirical_mfpt,
      c(
        common,
        list(
          cutoffs = empirical_cutoffs,
          domain = empirical_domain,
          reference_by = empirical_reference_by,
          prior_strength = empirical_prior_strength
        ),
        empirical_args
      )
    )
  )
}
