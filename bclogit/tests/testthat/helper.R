# ============================================================================
# Helper: generate simulated matched-pairs data
# ============================================================================
generate_test_data <- function(n_pairs = 100, p = 2, seed = 42) {
  set.seed(seed)
  n <- 2 * n_pairs
  strata <- rep(seq_len(n_pairs), each = 2)
  treatment <- rep(c(0, 1), n_pairs)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))

  # Generate outcomes: roughly half concordant, half discordant
  beta_true <- c(0.5, rep(0.3, p))
  eta <- cbind(treatment, X) %*% beta_true
  prob <- 1 / (1 + exp(-eta))
  y <- rbinom(n, 1, prob)

  list(
    y = y, X = X, treatment = treatment, strata = strata,
    n = n, p = p, n_pairs = n_pairs
  )
}

# Build a data.frame version from the list
make_df <- function(dat) {
  data.frame(y = dat$y, dat$X, treatment = dat$treatment, strata = dat$strata)
}

# Helper to call bclogit.default (unexported S3 method)
bclogit_default <- function(...) bclogit:::bclogit.default(...)

# Skip slow Stan tests unless BCLOGIT_RUN_STAN_TESTS env var is set
skip_if_no_stan <- function() {
  if (identical(Sys.getenv("BCLOGIT_RUN_STAN_TESTS"), "true")) return(invisible())
  skip("Stan tests skipped. Set BCLOGIT_RUN_STAN_TESTS=true to run.")
}
