test_that("empirical fitting returns expected transitions", {
  fit <- fit_empirical_cgm(
    glucose = c(60, 65, 100, 120, 200, 170),
    time = seq(0, 25, by = 5),
    zero_row = "self",
    return_states = TRUE
  )
  expect_s3_class(fit, "cgm_empirical_fit")
  expect_equal(sum(fit$transition_counts), 5)
  expect_equal(unname(fit$transition_counts[1, 1]), 1)
  expect_equal(unname(fit$transition_counts[1, 2]), 1)
  expect_equal(unname(fit$transition_counts[3, 2]), 1)
  expect_equal(nrow(fit$states), 6)
})

test_that("empirical smoothing adds reference pseudocounts", {
  reference <- matrix(1 / 3, 3, 3)
  fit <- fit_empirical_cgm(
    glucose = c(60, 65, 100, 120, 200, 170),
    time = seq(0, 25, by = 5),
    reference_transition = reference,
    prior_strength = 6
  )
  expect_equal(
    unname(fit$smoothed_counts),
    unname(fit$transition_counts + 2)
  )
  expect_equal(unname(rowSums(fit$transition_matrix)), rep(1, 3))
})

test_that("cohort empirical estimates do not cross series", {
  data <- data.frame(
    id = rep(c("A", "B"), each = 3),
    period = "baseline",
    age = "adults",
    time = rep(c(0, 5, 10), 2),
    glucose = c(60, 100, 120, 200, 170, 150)
  )
  result <- estimate_empirical_mfpt(
    data,
    id = "id",
    time = "time",
    glucose = "glucose",
    group = "period",
    stratum = "age",
    min_obs = 2,
    on_error = "stop"
  )
  expect_equal(nrow(result), 2)
  expect_equal(result$n_transitions, c(2, 2))
})
