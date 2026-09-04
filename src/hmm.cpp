#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

struct EStepResult {
  double loglik;
  NumericMatrix xi_sum;
  NumericVector pi_sum;
  NumericVector post_n;
  NumericMatrix posterior;
};

EStepResult run_estep(
    const NumericMatrix& B,
    const IntegerVector& starts,
    const IntegerVector& ends,
    const NumericVector& pi,
    const NumericMatrix& A,
    bool store_posterior) {
  const int T = B.nrow();
  const int K = B.ncol();
  const int S = starts.size();
  NumericMatrix xi_sum(K, K);
  NumericVector pi_sum(K);
  NumericVector post_n(K);
  NumericMatrix posterior(store_posterior ? T : 0, store_posterior ? K : 0);
  double loglik = 0.0;

  for (int segment = 0; segment < S; ++segment) {
    const int first = starts[segment] - 1;
    const int last = ends[segment] - 1;
    const int length = last - first + 1;
    if (length <= 0) {
      continue;
    }

    NumericMatrix alpha(length, K);
    NumericMatrix beta(length, K);
    NumericVector scale(length);

    double total = 0.0;
    for (int k = 0; k < K; ++k) {
      alpha(0, k) = pi[k] * B(first, k);
      total += alpha(0, k);
    }
    if (!(total > 0.0) || !R_finite(total)) {
      stop("The forward probability was not finite at a segment start.");
    }
    scale[0] = total;
    for (int k = 0; k < K; ++k) {
      alpha(0, k) /= total;
    }

    for (int t = 1; t < length; ++t) {
      total = 0.0;
      for (int k = 0; k < K; ++k) {
        double prediction = 0.0;
        for (int j = 0; j < K; ++j) {
          prediction += alpha(t - 1, j) * A(j, k);
        }
        alpha(t, k) = prediction * B(first + t, k);
        total += alpha(t, k);
      }
      if (!(total > 0.0) || !R_finite(total)) {
        stop("The forward probability was not finite.");
      }
      scale[t] = total;
      for (int k = 0; k < K; ++k) {
        alpha(t, k) /= total;
      }
    }

    for (int t = 0; t < length; ++t) {
      loglik += std::log(scale[t]);
    }
    for (int k = 0; k < K; ++k) {
      beta(length - 1, k) = 1.0;
    }
    for (int t = length - 2; t >= 0; --t) {
      for (int j = 0; j < K; ++j) {
        double value = 0.0;
        for (int k = 0; k < K; ++k) {
          value += A(j, k) * B(first + t + 1, k) * beta(t + 1, k);
        }
        beta(t, j) = value / scale[t + 1];
      }
    }

    for (int t = 0; t < length; ++t) {
      double gamma_total = 0.0;
      for (int k = 0; k < K; ++k) {
        gamma_total += alpha(t, k) * beta(t, k);
      }
      if (!(gamma_total > 0.0) || !R_finite(gamma_total)) {
        stop("The posterior state probability was not finite.");
      }
      for (int k = 0; k < K; ++k) {
        const double gamma = alpha(t, k) * beta(t, k) / gamma_total;
        post_n[k] += gamma;
        if (t == 0) {
          pi_sum[k] += gamma;
        }
        if (store_posterior) {
          posterior(first + t, k) = gamma;
        }
      }
    }

    for (int t = 0; t < length - 1; ++t) {
      double xi_total = 0.0;
      for (int j = 0; j < K; ++j) {
        for (int k = 0; k < K; ++k) {
          xi_total += alpha(t, j) * A(j, k) *
            B(first + t + 1, k) * beta(t + 1, k);
        }
      }
      if (!(xi_total > 0.0) || !R_finite(xi_total)) {
        stop("The posterior transition probability was not finite.");
      }
      for (int j = 0; j < K; ++j) {
        for (int k = 0; k < K; ++k) {
          xi_sum(j, k) += alpha(t, j) * A(j, k) *
            B(first + t + 1, k) * beta(t + 1, k) / xi_total;
        }
      }
    }
  }

  return EStepResult{
    loglik,
    xi_sum,
    pi_sum,
    post_n,
    posterior
  };
}

// [[Rcpp::export]]
List baum_welch_cpp(
    NumericMatrix B,
    IntegerVector starts,
    IntegerVector ends,
    NumericVector initial,
    NumericMatrix transition,
    int max_iter = 50,
    double tolerance = 1e-4,
    bool return_posterior = false) {
  const int T = B.nrow();
  const int K = B.ncol();
  const int S = starts.size();
  if (T < 1 || K < 2) {
    stop("B must have at least one row and two columns.");
  }
  if (starts.size() != ends.size() || S < 1) {
    stop("starts and ends must define at least one segment.");
  }
  if (initial.size() != K ||
      transition.nrow() != K || transition.ncol() != K) {
    stop("Initial and transition dimensions do not match B.");
  }
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < K; ++k) {
      if (!(B(t, k) > 0.0) || !R_finite(B(t, k))) {
        stop("All emission scores must be positive and finite.");
      }
    }
  }
  for (int s = 0; s < S; ++s) {
    if (starts[s] < 1 || ends[s] > T || starts[s] > ends[s]) {
      stop("Invalid segment bounds.");
    }
    if (s > 0 && starts[s] <= ends[s - 1]) {
      stop("Segments must be ordered and non-overlapping.");
    }
  }

  NumericVector pi = clone(initial);
  NumericMatrix A = clone(transition);
  double previous_loglik = R_NegInf;
  bool converged = false;
  int iterations = 0;

  for (int iteration = 0; iteration < max_iter; ++iteration) {
    EStepResult step = run_estep(B, starts, ends, pi, A, false);
    iterations = iteration + 1;
    double pi_total = 0.0;
    for (int k = 0; k < K; ++k) {
      pi_total += step.pi_sum[k];
    }
    if (pi_total > 0.0) {
      for (int k = 0; k < K; ++k) {
        pi[k] = step.pi_sum[k] / pi_total;
      }
    }

    for (int j = 0; j < K; ++j) {
      double row_total = 0.0;
      for (int k = 0; k < K; ++k) {
        row_total += step.xi_sum(j, k);
      }
      if (row_total > 0.0) {
        for (int k = 0; k < K; ++k) {
          A(j, k) = step.xi_sum(j, k) / row_total;
        }
      }
    }
    if (R_finite(previous_loglik) &&
        std::fabs(step.loglik - previous_loglik) < tolerance) {
      converged = true;
      break;
    }
    previous_loglik = step.loglik;
  }

  EStepResult final_step = run_estep(
    B,
    starts,
    ends,
    pi,
    A,
    return_posterior
  );
  return List::create(
    Named("transition_matrix") = A,
    Named("initial_probabilities") = pi,
    Named("posterior_state_counts") = final_step.post_n,
    Named("posterior") = final_step.posterior,
    Named("loglik_scaled") = final_step.loglik,
    Named("iterations") = iterations,
    Named("converged") = converged
  );
}
