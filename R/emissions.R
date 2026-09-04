.fit_beta_uniform <- function(x, support, shape_min, epsilon,
                              boundary_nudge, max_fit_n, seed) {
  lower <- support[1L]
  upper <- support[2L]
  x <- x[is.finite(x) & x >= lower & x <= upper]
  if (length(x) < 5L) {
    .stopf("At least five observations are required to fit each state.")
  }
  n_total <- length(x)
  if (length(x) > max_fit_n) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    x <- sample(x, max_fit_n)
  }

  width <- upper - lower
  nudge <- min(boundary_nudge, width / 4)
  adjusted <- x
  adjusted[adjusted <= lower] <- lower + nudge
  adjusted[adjusted >= upper] <- upper - nudge
  u <- (adjusted - lower) / width
  u <- pmin(pmax(u, epsilon), 1 - epsilon)

  objective <- function(par) {
    uniform_weight <- stats::plogis(par[1L])
    alpha <- shape_min[1L] + exp(par[2L])
    beta <- shape_min[2L] + exp(par[3L])
    density <- uniform_weight +
      (1 - uniform_weight) * stats::dbeta(u, alpha, beta)
    -sum(log(pmax(density, 1e-12)))
  }

  fitted <- stats::optim(
    c(stats::qlogis(0.3), log(1.2), log(1.2)),
    objective,
    method = "Nelder-Mead"
  )
  list(
    uniform_weight = stats::plogis(fitted$par[1L]),
    alpha = shape_min[1L] + exp(fitted$par[2L]),
    beta = shape_min[2L] + exp(fitted$par[3L]),
    n = n_total,
    n_fitted = length(x),
    convergence = fitted$convergence,
    objective = fitted$value
  )
}

.beta_uniform_density <- function(x, support, fit, epsilon,
                                  boundary_nudge) {
  lower <- support[1L]
  upper <- support[2L]
  density <- numeric(length(x))
  inside <- is.finite(x) & x >= lower & x <= upper
  if (!any(inside)) {
    return(density)
  }
  width <- upper - lower
  nudge <- min(boundary_nudge, width / 4)
  adjusted <- x[inside]
  adjusted[adjusted <= lower] <- lower + nudge
  adjusted[adjusted >= upper] <- upper - nudge
  u <- (adjusted - lower) / width
  u <- pmin(pmax(u, epsilon), 1 - epsilon)
  density[inside] <- (
    fit$uniform_weight +
      (1 - fit$uniform_weight) * stats::dbeta(u, fit$alpha, fit$beta)
  ) / width
  density
}

#' Fit pooled beta-uniform emission templates
#'
#' Fits one beta-uniform mixture density for each of three ordered CGM states.
#' Observations are assigned to states using hard cutoffs for template fitting,
#' while the fitted supports extend into adjacent states by overlap mg/dL.
#' State prevalences in the reference data are retained for overlap weighting.
#'
#' @param data Reference data frame, usually pooled baseline CGM observations.
#' @param glucose Name of the numeric glucose column.
#' @param stratum Optional name of a column defining separate templates, such
#'   as an adult/child indicator.
#' @param cutoffs Two increasing cutoffs or a named list of cutoff pairs
#'   indexed by stratum.
#' @param domain Admissible glucose domain. Values at or outside its endpoints
#'   are excluded.
#' @param overlap Half-width in mg/dL added to each side of an adjacent-state
#'   boundary.
#' @param epsilon Numerical boundary adjustment on the unit interval.
#' @param boundary_nudge Amount in mg/dL used to move exact support endpoints
#'   into the support before beta-density evaluation.
#' @param shape_min Three by two matrix of lower bounds for alpha and beta.
#' @param max_fit_n Maximum observations used to optimize each state density.
#' @param min_state_n Minimum reference observations required in each state.
#' @param seed Seed used only when a state is subsampled for fitting.
#'
#' @return A cgm_emission_templates object containing one template per stratum.
#' @export
#'
#' @examples
#' data <- simulate_cgm_data(n_participants = 6, n_per_period = 96, seed = 1)
#' baseline <- data[data$period == "baseline", ]
#' templates <- fit_emission_templates(
#'   baseline,
#'   glucose = "glucose",
#'   stratum = "age_group",
#'   cutoffs = list(adults = c(70, 180), children = c(70, 180)),
#'   max_fit_n = 5000,
#'   min_state_n = 5
#' )
#' templates
fit_emission_templates <- function(
    data,
    glucose = "glucose",
    stratum = NULL,
    cutoffs = c(70, 180),
    domain = c(40, 401),
    overlap = 10,
    epsilon = 1e-6,
    boundary_nudge = 0.5,
    shape_min = matrix(
      c(0.2, 1.0, 0.8, 0.8, 0.8, 0.8),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(.state_names, c("alpha", "beta"))
    ),
    max_fit_n = 200000L,
    min_state_n = 20L,
    seed = 1L) {
  .assert_scalar_character(glucose, "glucose")
  .assert_scalar_character(stratum, "stratum", allow_null = TRUE)
  .assert_columns(data, c(glucose, stratum))
  domain <- .validate_domain(domain)
  if (!is.numeric(overlap) || length(overlap) != 1L ||
      !is.finite(overlap) || overlap < 0) {
    .stopf("overlap must be one non-negative finite number.")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      epsilon <= 0 || epsilon >= 0.5) {
    .stopf("epsilon must lie between 0 and 0.5.")
  }
  if (!is.numeric(boundary_nudge) || length(boundary_nudge) != 1L ||
      !is.finite(boundary_nudge) || boundary_nudge <= 0) {
    .stopf("boundary_nudge must be one positive finite number.")
  }
  shape_min <- as.matrix(shape_min)
  if (!is.numeric(shape_min) || !identical(dim(shape_min), c(3L, 2L)) ||
      any(!is.finite(shape_min)) || any(shape_min < 0)) {
    .stopf("shape_min must be a finite non-negative 3 by 2 matrix.")
  }
  max_fit_n <- as.integer(max_fit_n)
  min_state_n <- as.integer(min_state_n)
  seed <- as.integer(seed)
  if (length(max_fit_n) != 1L || is.na(max_fit_n) || max_fit_n < 5L ||
      length(min_state_n) != 1L || is.na(min_state_n) || min_state_n < 5L) {
    .stopf("max_fit_n and min_state_n must be integers of at least 5.")
  }
  if (length(seed) != 1L || is.na(seed)) {
    .stopf("seed must be one integer.")
  }

  x_all <- data[[glucose]]
  if (!is.numeric(x_all)) {
    .stopf("The glucose column must be numeric.")
  }
  if (is.null(stratum)) {
    strata <- ".all"
  } else {
    if (anyNA(data[[stratum]])) {
      .stopf("The stratum column cannot contain missing values.")
    }
    strata <- unique(as.character(data[[stratum]]))
  }

  templates <- vector("list", length(strata))
  names(templates) <- strata
  for (index in seq_along(strata)) {
    key <- strata[index]
    rows <- if (is.null(stratum)) {
      rep(TRUE, nrow(data))
    } else {
      as.character(data[[stratum]]) == key
    }
    x <- x_all[rows]
    x <- x[is.finite(x) & x > domain[1L] & x < domain[2L]]
    state_cutoffs <- .resolve_cutoffs(
      cutoffs,
      if (is.null(stratum)) NULL else key,
      domain
    )
    states <- assign_cgm_states(x, state_cutoffs)
    state_n <- tabulate(states, nbins = 3L)
    if (any(state_n < min_state_n)) {
      .stopf(
        "Stratum '%s' has fewer than %d observations in at least one state.",
        key,
        min_state_n
      )
    }
    supports <- rbind(
      lower = c(domain[1L], min(domain[2L], state_cutoffs[1L] + overlap)),
      middle = c(
        max(domain[1L], state_cutoffs[1L] - overlap),
        min(domain[2L], state_cutoffs[2L] + overlap)
      ),
      upper = c(max(domain[1L], state_cutoffs[2L] - overlap), domain[2L])
    )
    colnames(supports) <- c("lower", "upper")

    fits <- vector("list", 3L)
    names(fits) <- .state_names
    for (state in 1:3) {
      fits[[state]] <- .fit_beta_uniform(
        x[states == state],
        supports[state, ],
        shape_min[state, ],
        epsilon,
        boundary_nudge,
        max_fit_n,
        seed + index * 10L + state
      )
    }

    templates[[key]] <- structure(
      list(
        stratum = key,
        cutoffs = state_cutoffs,
        domain = domain,
        overlap = overlap,
        supports = supports,
        fits = fits,
        prevalence = stats::setNames(state_n / sum(state_n), .state_names),
        n_reference = length(x),
        epsilon = epsilon,
        boundary_nudge = boundary_nudge
      ),
      class = "cgm_emission_template"
    )
  }

  structure(
    list(
      templates = templates,
      glucose = glucose,
      stratum = stratum,
      settings = list(
        cutoffs = cutoffs,
        domain = domain,
        overlap = overlap,
        epsilon = epsilon,
        boundary_nudge = boundary_nudge,
        shape_min = shape_min,
        max_fit_n = max_fit_n,
        min_state_n = min_state_n,
        seed = seed
      ),
      call = match.call()
    ),
    class = "cgm_emission_templates"
  )
}

#' Evaluate HMM emission log-scores
#'
#' Evaluates the fitted beta-uniform densities. Within each overlap, the scores
#' for the two adjacent states are multiplied by their normalized reference
#' prevalences.
#'
#' @param glucose Numeric glucose measurements.
#' @param templates Object returned by fit_emission_templates, or one individual
#'   cgm_emission_template.
#' @param stratum Stratum identifying the template to use.
#' @param prevalence_weight If TRUE, apply normalized prevalence weights in
#'   overlap regions.
#' @param min_density Positive lower bound used before taking logarithms.
#'
#' @return A matrix of log-scores with lower, middle, and upper columns.
#' @export
emission_log_scores <- function(glucose, templates, stratum = NULL,
                                prevalence_weight = TRUE,
                                min_density = 10^-300) {
  if (!is.numeric(glucose)) {
    .stopf("glucose must be numeric.")
  }
  if (!is.numeric(min_density) || length(min_density) != 1L ||
      !is.finite(min_density) || min_density <= 0) {
    .stopf("min_density must be one positive finite number.")
  }
  template <- .resolve_template(templates, stratum)
  raw <- vapply(
    seq_len(3L),
    function(state) {
      .beta_uniform_density(
        glucose,
        template$supports[state, ],
        template$fits[[state]],
        template$epsilon,
        template$boundary_nudge
      )
    },
    numeric(length(glucose))
  )
  colnames(raw) <- .state_names
  log_scores <- log(pmax(raw, min_density))

  if (isTRUE(prevalence_weight)) {
    p <- pmax(template$prevalence, min_density)
    weights_12 <- p[1:2] / sum(p[1:2])
    weights_23 <- p[2:3] / sum(p[2:3])
    overlap_12 <- glucose >= max(template$supports[1L, 1L],
                                 template$supports[2L, 1L]) &
      glucose <= min(template$supports[1L, 2L],
                     template$supports[2L, 2L])
    overlap_23 <- glucose >= max(template$supports[2L, 1L],
                                 template$supports[3L, 1L]) &
      glucose <= min(template$supports[2L, 2L],
                     template$supports[3L, 2L])
    log_scores[overlap_12, 1L] <- log_scores[overlap_12, 1L] +
      log(weights_12[1L])
    log_scores[overlap_12, 2L] <- log_scores[overlap_12, 2L] +
      log(weights_12[2L])
    log_scores[overlap_23, 2L] <- log_scores[overlap_23, 2L] +
      log(weights_23[1L])
    log_scores[overlap_23, 3L] <- log_scores[overlap_23, 3L] +
      log(weights_23[2L])
  }
  log_scores
}

#' Summarize fitted emission templates
#'
#' @param object A cgm_emission_templates object.
#' @param ... Unused.
#'
#' @return A data frame with supports, mixture parameters, and prevalences.
#' @export
summary.cgm_emission_templates <- function(object, ...) {
  rows <- list()
  counter <- 0L
  for (key in names(object$templates)) {
    template <- object$templates[[key]]
    for (state in seq_len(3L)) {
      counter <- counter + 1L
      fit <- template$fits[[state]]
      rows[[counter]] <- data.frame(
        stratum = key,
        state = .state_names[state],
        support_lower = template$supports[state, 1L],
        support_upper = template$supports[state, 2L],
        cutoff_lower = template$cutoffs[1L],
        cutoff_upper = template$cutoffs[2L],
        prevalence = unname(template$prevalence[state]),
        uniform_weight = fit$uniform_weight,
        alpha = fit$alpha,
        beta = fit$beta,
        n = fit$n,
        convergence = fit$convergence,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
