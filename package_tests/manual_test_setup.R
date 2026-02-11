# Manual Test Setup Script
# This script sets up the data environment similar to package_tests/parallel_pacakage_tests_from_sratch.R
# It creates the variables: y, X, w, strat, and runs a basic bclogit model.

options(error = recover)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, data.table, MASS, bclogit)

# --- Parameters (Default values from parallel test script) ---
n <- 500 # Total sample size (number of individuals)
p <- 6 # Number of covariates
beta_T <- 0.5 # Treatment effect size
X_style <- "non-correlated" # Options: "correlated", "non-correlated"
true_function <- "linear" # Options: "linear", "non-linear"
regress_on_X <- "some" # Options: "all", "one", "some", "none"

# --- Data Generation Logic ---

# 1. Generate X matrix
if (X_style == "correlated") {
    Sigma <- 1 * (matrix(0.5, nrow = p, ncol = p) + diag(1 - 0.5, p))
    X_raw <- MASS::mvrnorm(n / 2, rep(0, p), Sigma)
    X_raw <- pnorm(X_raw)
    X_raw <- matrix(2 * X_raw - 1, ncol = p)
} else {
    # non-correlated
    X_raw <- matrix(runif((n / 2) * p, min = -1, max = 1), ncol = p)
    X_plus_eps <- X_raw + matrix(rnorm((n / 2) * p, 0, 0.05), ncol = p)
    combined <- rbind(X_raw, X_plus_eps)

    # Shuffle pairs to mix them up but keep pairs together logically?
    # Wait, the original code orders by c(1:(n/2), 1:(n/2)) which interleaves them effectively if we consider strat logic later.
    # Let's match original exactly:
    ids <- order(c(1:(n / 2), 1:(n / 2)))
    X_raw <- combined[ids, ]

    # Special modification for column 1
    X_raw[, 1] <- runif(n, min = -1, max = 1)
}

# 2. Generate w (treatment) and strat (strata)
# w is randomized within pairs (0,1) or (1,0)
w <- c(rbind(replicate(n / 2, sample(c(0, 1)), simplify = TRUE)))
strat <- rep(1:(n / 2), each = 2)

# 3. Calculate Probabilities and Outcome (y)
if (true_function == "linear") {
    beta_X_value <- if (p == 6) {
        1.25
    } else {
        0.75
    }
    beta_X <- rep(beta_X_value, p)
    beta_0 <- -0.5
    probs <- 1 / (1 + exp(-(beta_0 + (as.matrix(X_raw) %*% beta_X) + beta_T * w)))
} else {
    # non-linear
    f_x <- sin(pi * X_raw[, 1] * X_raw[, 2]) + X_raw[, 3]^3 + X_raw[, 4]^2 + X_raw[, 5]^2
    probs <- 1 / (1 + exp(-(f_x + beta_T * w)))
}

y <- rbinom(n, 1, probs)

# 4. Select Covariates for Regression (X)
if (regress_on_X == "one") {
    X <- X_raw[, 1, drop = FALSE]
} else if (regress_on_X == "some") {
    if (p == 20) {
        X <- X_raw[, 1:5, drop = FALSE]
    } else {
        X <- X_raw[, c(1, 2), drop = FALSE]
    }
} else if (regress_on_X == "all") {
    X <- X_raw
} else {
    # none
    X <- NULL
}


 fit <- bclogit(y ~ X, treatment = w, strat = strat)
 print(summary(fit))
