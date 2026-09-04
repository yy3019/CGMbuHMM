#' Simulate example three-state CGM data
#'
#' Generates synthetic participant-period CGM series for examples and testing.
#' The data do not represent trial participants or reproduce the study results.
#'
#' @param n_participants Number of synthetic participants.
#' @param n_per_period Number of observations in each baseline and RCT period.
#' @param step_minutes Minutes between scheduled observations.
#' @param seed Random-number seed.
#'
#' @return A data frame with participant, age-group, treatment, period, time,
#'   glucose, and latent-state columns.
#' @export
#'
#' @examples
#' synthetic_cgm <- simulate_cgm_data(n_participants = 4, n_per_period = 48)
#' head(synthetic_cgm)
simulate_cgm_data <- function(n_participants = 12L, n_per_period = 144L,
                              step_minutes = 5, seed = 1L) {
  n_participants <- as.integer(n_participants)
  n_per_period <- as.integer(n_per_period)
  if (length(n_participants) != 1L || is.na(n_participants) ||
      n_participants < 2L) {
    .stopf("n_participants must be an integer of at least 2.")
  }
  if (length(n_per_period) != 1L || is.na(n_per_period) ||
      n_per_period < 12L) {
    .stopf("n_per_period must be an integer of at least 12.")
  }
  if (!is.numeric(step_minutes) || length(step_minutes) != 1L ||
      !is.finite(step_minutes) || step_minutes <= 0) {
    .stopf("step_minutes must be one positive finite number.")
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    .stopf("seed must be one integer.")
  }

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

  ids <- sprintf("P%03d", seq_len(n_participants))
  age_groups <- rep(
    c("adults", "adults", "children", "children"),
    length.out = n_participants
  )
  treatments <- rep(
    c("Control", "BP", "Control", "BP"),
    length.out = n_participants
  )
  output <- vector("list", n_participants * 2L)
  counter <- 0L

  for (participant in seq_len(n_participants)) {
    for (period_index in seq_along(c("baseline", "rct"))) {
      period <- c("baseline", "rct")[period_index]
      age_group <- age_groups[participant]
      treatment <- treatments[participant]
      A <- matrix(
        c(
          0.80, 0.19, 0.01,
          0.03, 0.93, 0.04,
          0.01, 0.17, 0.82
        ),
        nrow = 3,
        byrow = TRUE
      )
      if (age_group == "children") {
        A[2L, ] <- c(0.04, 0.90, 0.06)
        A[3L, ] <- c(0.01, 0.14, 0.85)
      }
      if (period == "rct" && treatment == "BP") {
        A[1L, ] <- c(0.73, 0.26, 0.01)
        A[3L, ] <- c(0.01, 0.24, 0.75)
      }

      states <- integer(n_per_period)
      states[1L] <- (participant + period_index - 2L) %% 3L + 1L
      for (index in 2:n_per_period) {
        states[index] <- sample.int(3L, 1L, prob = A[states[index - 1L], ])
      }
      states[seq_len(9L)] <- rep(1:3, each = 3L)
      glucose <- numeric(n_per_period)
      for (state in 1:3) {
        where <- which(states == state)
        if (!length(where)) {
          next
        }
        if (state == 1L) {
          values <- stats::rnorm(length(where), 61, 7)
          glucose[where] <- pmin(pmax(values, 41), 70)
        } else if (state == 2L) {
          values <- stats::rnorm(
            length(where),
            if (age_group == "adults") 126 else 132,
            25
          )
          glucose[where] <- pmin(pmax(values, 70.1), 180)
        } else {
          values <- stats::rnorm(
            length(where),
            if (age_group == "adults") 220 else 232,
            30
          )
          glucose[where] <- pmin(pmax(values, 180.1), 400)
        }
      }
      glucose <- round(glucose + stats::rnorm(n_per_period, 0, 1.5), 1)
      glucose <- pmin(pmax(glucose, 40.1), 400.9)
      counter <- counter + 1L
      start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC") +
        (counter - 1L) * 24 * 60 * 60
      output[[counter]] <- data.frame(
        participant_id = ids[participant],
        age_group = age_group,
        treatment = treatment,
        period = period,
        timestamp = start + seq(0, by = step_minutes * 60,
                                length.out = n_per_period),
        glucose = glucose,
        latent_state = ordered(states, levels = 1:3, labels = .state_names),
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}
