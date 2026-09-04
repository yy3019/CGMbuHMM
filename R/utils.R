.state_names <- c("lower", "middle", "upper")

.stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

.assert_scalar_character <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(invisible(TRUE))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .stopf("%s must be a single non-empty column name.", name)
  }
  invisible(TRUE)
}

.assert_columns <- function(data, columns) {
  if (!is.data.frame(data)) {
    .stopf("data must be a data frame.")
  }
  columns <- unique(columns[!is.na(columns) & nzchar(columns)])
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .stopf("Missing required column(s): %s.", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

.validate_cutoff_vector <- function(cutoffs, domain = NULL) {
  if (!is.numeric(cutoffs) || length(cutoffs) != 2L ||
      any(!is.finite(cutoffs)) || cutoffs[1L] >= cutoffs[2L]) {
    .stopf("cutoffs must contain two increasing finite numbers.")
  }
  if (!is.null(domain) &&
      (cutoffs[1L] <= domain[1L] || cutoffs[2L] >= domain[2L])) {
    .stopf("Both cutoffs must lie strictly inside domain.")
  }
  as.numeric(cutoffs)
}

.validate_domain <- function(domain) {
  if (!is.numeric(domain) || length(domain) != 2L ||
      any(!is.finite(domain)) || domain[1L] >= domain[2L]) {
    .stopf("domain must contain two increasing finite numbers.")
  }
  as.numeric(domain)
}

.resolve_cutoffs <- function(cutoffs, stratum = NULL, domain = NULL) {
  if (is.numeric(cutoffs)) {
    return(.validate_cutoff_vector(cutoffs, domain))
  }
  if (!is.list(cutoffs) || is.null(names(cutoffs)) || is.null(stratum)) {
    .stopf("cutoffs must be a numeric pair or a named list indexed by stratum.")
  }
  key <- as.character(stratum)
  if (!key %in% names(cutoffs)) {
    .stopf("No cutoffs were supplied for stratum '%s'.", key)
  }
  .validate_cutoff_vector(cutoffs[[key]], domain)
}

.time_to_minutes <- function(time) {
  if (inherits(time, "POSIXt")) {
    return(as.numeric(time) / 60)
  }
  if (inherits(time, "Date")) {
    return(as.numeric(time) * 24 * 60)
  }
  if (inherits(time, "difftime")) {
    return(as.numeric(time, units = "mins"))
  }
  if (is.numeric(time)) {
    return(as.numeric(time))
  }
  .stopf("time must be POSIXct, POSIXlt, Date, difftime, or numeric minutes.")
}

.segment_bounds <- function(time_minutes, max_gap_minutes) {
  n <- length(time_minutes)
  if (!n) {
    return(list(starts = integer(), ends = integer()))
  }
  if (n == 1L) {
    return(list(starts = 1L, ends = 1L))
  }
  gaps <- diff(time_minutes)
  break_after <- which(!is.finite(gaps) | gaps < 0 | gaps > max_gap_minutes)
  list(
    starts = as.integer(c(1L, break_after + 1L)),
    ends = as.integer(c(break_after, n))
  )
}

.normalize_rows <- function(counts, zero_row = c("self", "uniform", "error")) {
  zero_row <- match.arg(zero_row)
  counts <- as.matrix(counts)
  out <- matrix(0, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  totals <- rowSums(counts)
  for (j in seq_len(nrow(counts))) {
    if (is.finite(totals[j]) && totals[j] > 0) {
      out[j, ] <- counts[j, ] / totals[j]
    } else if (zero_row == "self") {
      out[j, j] <- 1
    } else if (zero_row == "uniform") {
      out[j, ] <- 1 / ncol(counts)
    } else {
      .stopf("Transition row %d has no observed or prior transitions.", j)
    }
  }
  out
}

.group_key <- function(data, columns) {
  if (!length(columns)) {
    return(rep(".all", nrow(data)))
  }
  if (anyNA(data[columns])) {
    .stopf("Grouping columns cannot contain missing values.")
  }
  do.call(
    paste,
    c(lapply(data[columns], function(x) as.character(x)), sep = "\034")
  )
}

.split_indices <- function(data, columns) {
  split(seq_len(nrow(data)), .group_key(data, columns), drop = TRUE)
}

.one_value <- function(x, name) {
  values <- unique(x[!is.na(x)])
  if (length(values) != 1L) {
    .stopf("Each analysis group must contain exactly one value of %s.", name)
  }
  values[[1L]]
}

.resolve_template <- function(templates, stratum = NULL) {
  if (inherits(templates, "cgm_emission_template")) {
    return(templates)
  }
  if (!inherits(templates, "cgm_emission_templates")) {
    .stopf("templates must be created by fit_emission_templates().")
  }
  key <- if (is.null(stratum)) ".all" else as.character(stratum)
  if (!key %in% names(templates$templates)) {
    .stopf("No emission template was fitted for stratum '%s'.", key)
  }
  templates$templates[[key]]
}

.cutoffs_from_templates <- function(templates) {
  if (!inherits(templates, "cgm_emission_templates")) {
    .stopf("templates must be created by fit_emission_templates().")
  }
  values <- lapply(templates$templates, function(x) x$cutoffs)
  if (identical(names(values), ".all")) {
    values[[1L]]
  } else {
    values
  }
}

.validate_probability_matrix <- function(A, tolerance = 1e-8) {
  A <- as.matrix(A)
  if (!is.numeric(A) || nrow(A) != ncol(A) || nrow(A) < 2L ||
      any(!is.finite(A)) || any(A < -tolerance)) {
    .stopf("transition_matrix must be a finite, non-negative square matrix.")
  }
  A[A < 0] <- 0
  totals <- rowSums(A)
  if (any(abs(totals - 1) > tolerance)) {
    .stopf("Each row of transition_matrix must sum to 1.")
  }
  A / totals
}
