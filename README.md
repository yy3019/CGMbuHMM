# CGMbuHMM

`CGMbuHMM` implements the three-state continuous glucose monitoring (CGM)
methods used in our manuscript. It provides:

- an empirical transition-matrix estimator based on hard glucose cutoffs;
- a hidden Markov model (HMM) with pooled beta-uniform mixture emissions and
  overlapping state supports; and
- mean first-passage times (MFPTs) from the lower and upper states to the
  middle state.

The package accepts clinical or user-specified cutoffs. Emission templates can
be estimated separately for analysis strata, such as adults and children, from
a pooled reference period. The package contains synthetic example data only.

## Installation

After this folder is uploaded as a GitHub repository, install it with:

```r
install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/CGMbuHMM")
```

For local development, run:

```r
install.packages(".", repos = NULL, type = "source")
```

The package compiles a small C++ routine through `Rcpp`, so a working R build
toolchain is required.

## Data format

Long-format data should include one row per CGM observation and columns for:

- participant ID;
- observation time, as `POSIXct`, `Date`, `difftime`, or numeric minutes;
- glucose in mg/dL;
- any columns that identify an analysis series, such as period or treatment;
  and
- an optional template stratum, such as age group.

Observations separated by more than `max_gap_minutes` are placed in different
continuous segments. No transition is estimated across a segment boundary.
By default, values at or outside 40 and 401 mg/dL are excluded from both
methods, matching the analysis domain used in the manuscript.

## Quick start

```r
library(CGMbuHMM)

cgm <- simulate_cgm_data(
  n_participants = 8,
  n_per_period = 144,
  seed = 2026
)

# Fit emissions using pooled baseline observations within age group.
baseline <- cgm[cgm$period == "baseline", ]
templates <- fit_emission_templates(
  baseline,
  glucose = "glucose",
  stratum = "age_group",
  cutoffs = list(
    adults = c(70, 180),
    children = c(70, 180)
  ),
  overlap = 10,
  min_state_n = 5
)

# Apply the fixed templates to each participant-period series.
hmm_results <- estimate_hmm_mfpt(
  cgm,
  templates,
  id = "participant_id",
  time = "timestamp",
  glucose = "glucose",
  group = c("period", "treatment"),
  stratum = "age_group",
  min_obs = 50
)

# Estimate hard-state empirical transitions. A strength of 20 adds one pooled
# reference transition row, with total weight 20, to each participant row.
empirical_results <- estimate_empirical_mfpt(
  cgm,
  id = "participant_id",
  time = "timestamp",
  glucose = "glucose",
  group = c("period", "treatment"),
  stratum = "age_group",
  cutoffs = list(
    adults = c(70, 180),
    children = c(70, 180)
  ),
  reference_by = c("age_group", "period"),
  prior_strength = 20,
  min_obs = 50
)

head(hmm_results)
head(empirical_results)
```

Both result tables report MFPT in minutes as
`mfpt_lower_to_middle` and `mfpt_upper_to_middle`.

## Fit one series

The lower-level functions return transition matrices and diagnostic details:

```r
one <- cgm[cgm$participant_id == "P001" & cgm$period == "rct", ]

hmm_fit <- fit_hmm_cgm(
  one$glucose,
  one$timestamp,
  templates,
  stratum = one$age_group[1],
  return_posterior = TRUE
)
hmm_fit$transition_matrix
hmm_fit$mfpt
head(hmm_fit$posterior)

empirical_fit <- fit_empirical_cgm(
  one$glucose,
  one$timestamp,
  cutoffs = c(70, 180),
  return_states = TRUE
)
empirical_fit$transition_matrix
empirical_fit$mfpt
```

## Custom or data-driven cutoffs

Any pair of increasing cutoffs can replace 70 and 180 mg/dL. When cutoffs vary
by stratum, pass a named list to both template fitting and empirical estimation:

```r
study_cutoffs <- list(
  adults = c(82, 166),
  children = c(85, 172)
)

custom_templates <- fit_emission_templates(
  baseline,
  stratum = "age_group",
  cutoffs = study_cutoffs,
  overlap = 10,
  min_state_n = 5
)
```

`CGMbuHMM` applies supplied cutoffs but does not optimize them. This keeps cutoff
selection separate from transition and MFPT estimation and helps prevent using
the same participant-period outcomes both to choose and evaluate thresholds.

## Method details

For the empirical method, `assign_cgm_states()` applies hard cutoffs and
`transition_counts()` counts only consecutive observations in the same
continuous segment. With smoothing strength `lambda`, a participant's count
row is augmented by `lambda` times the corresponding row of the pooled
reference transition matrix before normalization.

For the HMM, `fit_emission_templates()` fits a beta-uniform mixture within each
state support using reference data. Adjacent supports overlap around each
cutoff. Within an overlap, the two adjacent emission scores are multiplied by
their normalized reference-state prevalences. `fit_hmm_cgm()` then holds these
emissions fixed and estimates the segment-start distribution and transition
matrix by the scaled Baum-Welch expectation-maximization algorithm.

`mfpt()` solves the finite-state first-step equations and reports minutes. It
returns `Inf` when the target state is not reached with probability one.

## Reproducibility

Package tests can be run with:

```r
testthat::test_local()
```

A manuscript can cite the implementation with wording such as:

> All analyses were implemented in R using the CGMbuHMM package and
> author-written analysis scripts. The package source code is available at
> [GitHub repository URL].

Replace the placeholder installation path and repository URL after the GitHub
repository is created.

## License

MIT
