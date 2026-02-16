#!/usr/bin/env Rscript
# Benchmark bclogit::clogit vs survival::clogit

library(bclogit)
library(survival)

set.seed(2024)

generate_matched_pair_data <- function(n_pairs, n_covariates = 3) {
  n <- 2 * n_pairs
  strata <- rep(1:n_pairs, each = 2)
  treatment <- rep(c(0, 1), times = n_pairs)

  X <- matrix(rnorm(n * n_covariates), nrow = n, ncol = n_covariates)
  colnames(X) <- paste0("x", 1:n_covariates)

  beta_trt <- 0.5
  beta_x <- seq(0.3, by = 0.2, length.out = n_covariates)

  linpred <- beta_trt * treatment + X %*% beta_x
  prob <- plogis(linpred)
  y <- rbinom(n, 1, prob)

  data.frame(y = y, treatment = treatment, X, strata = strata)
}

pair_sizes <- c(500, 1000, 5000, 10000, 50000)
n_reps <- 5

results <- data.frame(
  n_pairs = integer(),
  bclogit_all_time = numeric(),
  bclogit_trt_time = numeric(),
  bclogit_none_time = numeric(),
  survival_time = numeric(),
  speedup_all = numeric(),
  speedup_trt = numeric(),
  speedup_none = numeric(),
  max_coef_diff = numeric(),
  max_se_diff = numeric(),
  stringsAsFactors = FALSE
)

cat("=============================================================\n")
cat("Benchmark: bclogit::clogit vs survival::clogit\n")
cat("=============================================================\n\n")

for (np in pair_sizes) {
  cat(sprintf("--- %d matched pairs ---\n", np))

  dat <- generate_matched_pair_data(np)
  covariate_cols <- grep("^x", names(dat), value = TRUE)
  fmla_str <- paste("y ~", paste(covariate_cols, collapse = " + "))

  # bclogit::clogit with do_inference_on_var = "all"
  bc_all_times <- numeric(n_reps)
  for (r in 1:n_reps) {
    bc_all_times[r] <- system.time({
      fit_bc <- bclogit::clogit(
        as.formula(fmla_str),
        data = dat,
        treatment = treatment,
        strata = strata,
        do_inference_on_var = "all"
      )
    })["elapsed"]
  }

  # bclogit::clogit with do_inference_on_var = 1 (treatment only)
  bc_trt_times <- numeric(n_reps)
  for (r in 1:n_reps) {
    bc_trt_times[r] <- system.time({
      fit_bc_trt <- bclogit::clogit(
        as.formula(fmla_str),
        data = dat,
        treatment = treatment,
        strata = strata,
        do_inference_on_var = 1
      )
    })["elapsed"]
  }

  # bclogit::clogit with do_inference_on_var = "none"
  bc_none_times <- numeric(n_reps)
  for (r in 1:n_reps) {
    bc_none_times[r] <- system.time({
      fit_bc_none <- bclogit::clogit(
        as.formula(fmla_str),
        data = dat,
        treatment = treatment,
        strata = strata,
        do_inference_on_var = "none"
      )
    })["elapsed"]
  }

  # survival::clogit
  surv_fmla <- as.formula(paste("y ~", paste(covariate_cols, collapse = " + "), "+ treatment + strata(strata)"))
  sv_times <- numeric(n_reps)
  for (r in 1:n_reps) {
    sv_times[r] <- system.time({
      fit_sv <- survival::clogit(surv_fmla, data = dat)
    })["elapsed"]
  }

  # Compare coefficients — align by common names
  bc_coefs <- coef(fit_bc)
  sv_coefs <- coef(fit_sv)
  common_names <- intersect(names(sv_coefs), names(bc_coefs))
  coef_diff <- max(abs(bc_coefs[common_names] - sv_coefs[common_names]))

  # Compare SEs — align by common names
  bc_se <- fit_bc$se
  names(bc_se) <- names(bc_coefs)
  sv_se <- sqrt(diag(vcov(fit_sv)))
  se_diff <- max(abs(bc_se[common_names] - sv_se[common_names]))

  bc_all_med <- median(bc_all_times)
  bc_trt_med <- median(bc_trt_times)
  bc_none_med <- median(bc_none_times)
  sv_med <- median(sv_times)

  cat(sprintf("  bclogit (all SEs)  median: %.4f s  (%.1fx vs survival)\n", bc_all_med, sv_med / bc_all_med))
  cat(sprintf("  bclogit (trt SE)   median: %.4f s  (%.1fx vs survival)\n", bc_trt_med, sv_med / bc_trt_med))
  cat(sprintf("  bclogit (no SEs)   median: %.4f s  (%.1fx vs survival)\n", bc_none_med, sv_med / bc_none_med))
  cat(sprintf("  survival::clogit   median: %.4f s\n", sv_med))
  cat(sprintf("  Max |coef diff|: %.2e\n", coef_diff))
  cat(sprintf("  Max |SE diff|:   %.2e\n\n", se_diff))

  results <- rbind(results, data.frame(
    n_pairs = np,
    bclogit_all_time = bc_all_med,
    bclogit_trt_time = bc_trt_med,
    bclogit_none_time = bc_none_med,
    survival_time = sv_med,
    speedup_all = sv_med / bc_all_med,
    speedup_trt = sv_med / bc_trt_med,
    speedup_none = sv_med / bc_none_med,
    max_coef_diff = coef_diff,
    max_se_diff = se_diff
  ))
}

cat("\n=============================================================\n")
cat("Summary Table\n")
cat("=============================================================\n")
print(results, row.names = FALSE)

# Spot-check: print coefficients side by side for the last run
cat("\n--- Coefficient comparison (last run, n_pairs =", tail(pair_sizes, 1), ") ---\n")
cat("\nbclogit::clogit coefficients:\n")
print(coef(fit_bc))
cat("\nsurvival::clogit coefficients:\n")
print(coef(fit_sv))
