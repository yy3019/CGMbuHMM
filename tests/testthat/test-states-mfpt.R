test_that("state assignment follows the two cutoffs", {
  expect_equal(
    assign_cgm_states(c(55, 70, 70.1, 180, 181, NA_real_)),
    c(1L, 1L, 2L, 2L, 3L, NA_integer_)
  )
})

test_that("transition counts respect segment gaps", {
  counts <- transition_counts(
    states = c(1, 2, 3, 2),
    time = c(0, 5, 20, 25),
    max_gap_minutes = 6
  )
  expect_equal(unname(counts[1, 2]), 1)
  expect_equal(unname(counts[3, 2]), 1)
  expect_equal(sum(counts), 2)
})

test_that("MFPT solves first-step equations", {
  A <- matrix(
    c(
      0.7, 0.3, 0,
      0.1, 0.8, 0.1,
      0, 0.4, 0.6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("lower", "middle", "upper"),
      c("lower", "middle", "upper")
    )
  )
  expect_equal(unname(mfpt(A, "middle")), c(50 / 3, 0, 12.5))
})

test_that("MFPT is infinite when the target cannot be reached", {
  A <- diag(3)
  dimnames(A) <- list(
    c("lower", "middle", "upper"),
    c("lower", "middle", "upper")
  )
  result <- mfpt(A, "middle")
  expect_true(is.infinite(result["lower"]))
  expect_true(is.infinite(result["upper"]))
  expect_equal(unname(result["middle"]), 0)
})
