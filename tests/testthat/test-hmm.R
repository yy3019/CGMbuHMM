test_that("emission templates and one HMM fit work together", {
  data <- simulate_cgm_data(
    n_participants = 4,
    n_per_period = 60,
    seed = 19
  )
  baseline <- data[data$period == "baseline", ]
  templates <- fit_emission_templates(
    baseline,
    stratum = "age_group",
    min_state_n = 5,
    max_fit_n = 1000
  )
  expect_s3_class(templates, "cgm_emission_templates")
  expect_equal(nrow(summary(templates)), 6)

  one <- data[data$participant_id == "P001" & data$period == "rct", ]
  fit <- fit_hmm_cgm(
    one$glucose,
    one$timestamp,
    templates,
    stratum = "adults",
    max_iter = 30,
    return_posterior = TRUE
  )
  expect_s3_class(fit, "cgm_hmm_fit")
  expect_equal(
    unname(rowSums(fit$transition_matrix)),
    rep(1, 3),
    tolerance = 1e-8
  )
  expect_equal(nrow(fit$posterior), nrow(one))
  expect_named(fit$mfpt, c("lower", "middle", "upper"))
})

test_that("combined wrapper returns both methods", {
  data <- simulate_cgm_data(
    n_participants = 4,
    n_per_period = 48,
    seed = 23
  )
  templates <- fit_emission_templates(
    data[data$period == "baseline", ],
    stratum = "age_group",
    min_state_n = 5,
    max_fit_n = 1000
  )
  one <- data[data$participant_id == "P001", ]
  result <- estimate_cgm_mfpt(
    one,
    templates,
    id = "participant_id",
    time = "timestamp",
    glucose = "glucose",
    group = "period",
    stratum = "age_group",
    min_obs = 24,
    hmm_args = list(max_iter = 20),
    empirical_prior_strength = 0
  )
  expect_named(result, c("hmm", "empirical"))
  expect_equal(nrow(result$hmm), 2)
  expect_equal(nrow(result$empirical), 2)
})
