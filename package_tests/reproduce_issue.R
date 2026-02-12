# Reproduction script for checking Median in summary
options(error = recover)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, data.table, MASS, devtools)

# Load the local package
devtools::load_all("bclogit")

# --- Parameters (Default values from parallel test script) ---
n <- 200 # Smaller sample size for speed
p <- 2 # Fewer covariates for speed
beta_T <- 0.5
X_style <- "non-correlated"
true_function <- "linear"
regress_on_X <- "some"

# --- Data Generation Logic ---

# 1. Generate X matrix
X_raw <- matrix(runif((n / 2) * p, min = -1, max = 1), ncol = p)
X_plus_eps <- X_raw + matrix(rnorm((n / 2) * p, 0, 0.05), ncol = p)
combined <- rbind(X_raw, X_plus_eps)
ids <- order(c(1:(n / 2), 1:(n / 2)))
X_raw <- combined[ids, ]
X_raw[, 1] <- runif(n, min = -1, max = 1)

# 2. Generate w (treatment) and strat (strata)
w <- c(rbind(replicate(n / 2, sample(c(0, 1)), simplify = TRUE)))
strat <- rep(1:(n / 2), each = 2)

# 3. Calculate Probabilities and Outcome (y)
beta_X <- rep(0.75, p)
beta_0 <- -0.5
probs <- 1 / (1 + exp(-(beta_0 + (as.matrix(X_raw) %*% beta_X) + beta_T * w)))
y <- rbinom(n, 1, probs)

# 4. Select Covariates for Regression (X)
X <- X_raw

fit <- bclogit(y ~ X, treatment = w, strat = strat)
summ <- summary(fit)
print(summ)

# Check if Median column exists
if ("Median" %in% colnames(summ$coefficients)) {
    print("SUCCESS: Median column found in summary.")
} else {
    print("FAILURE: Median column NOT found in summary.")
}
