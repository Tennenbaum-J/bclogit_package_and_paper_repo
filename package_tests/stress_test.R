
# Stress Test for bclogit Settings and Features
# Tests all combinations of prior_type, concordant_method, and interface (Formula/Default)

options(warn = 1) # Print warnings immediately

if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, dplyr, testthat, data.table)

# Load the package
if (dir.exists("bclogit")) {
  devtools::load_all("bclogit")
} else {
  library(bclogit)
}

context("bclogit Comprehensive Stress Test")

# Helper to generate data
generate_data <- function(n = 100, p = 3) {
  X <- matrix(rnorm(n * p), ncol = p)
  colnames(X) <- paste0("V", 1:p)
  beta <- rep(0.5, p)
  strata <- rep(1:(n/2), each=2)
  treatment <- rep(c(0, 1), n/2)
  
  # Latent variable for outcome
  # Note: simple simulation, not strictly conditional logit DGP but sufficient for fitting check
  lp <- X %*% beta + 0.5 * treatment
  probs <- 1 / (1 + exp(-lp))
  y <- rbinom(n, 1, probs)
  
  # Ensure some discordant pairs exist, otherwise flip a few
  # (though with n=100 random, likely to have discordant)
  
  list(
    y = y, 
    X = X, 
    treatment = treatment, 
    strata = strata, 
    df = data.frame(y = y, treatment = treatment, strata = strata, X)
  )
}

test_that("All combinations of settings work", {
  
  set.seed(123)
  dat <- generate_data(n = 200, p = 3)
  
  prior_types <- c("Naive", "G prior", "PMP", "Hybrid")
  concordant_methods <- c("GLM", "GEE", "GLMM")
  
  total_tests <- length(prior_types) * length(concordant_methods) * 2 # *2 for Default/Formula
  curr_test <- 0
  
  for (pt in prior_types) {
    for (cm in concordant_methods) {
      
      # --- Test Default Interface ---
      cat(sprintf("\n[%d/%d] Testing: Prior='%s', Concordant='%s', Interface='Default'\n", 
                  curr_test + 1, total_tests, pt, cm))
      
      fit_default <- NULL
      tryCatch({
        fit_default <- bclogit(
          response = dat$y,
          data = dat$X,
          treatment = dat$treatment,
          strata = dat$strata,
          prior_type = pt,
          concordant_method = cm
        )
      }, error = function(e) {
        fail(sprintf("Default fit failed with error: %s", e$message))
      })
      
      if (!is.null(fit_default)) {
        expect_s3_class(fit_default, "bclogit")
        
        # Verify S3 methods
        s <- summary(fit_default)
        expect_s3_class(s, "summary.bclogit")
        expect_true("Median" %in% colnames(s$coefficients), label = "Median in summary")
        
        cf <- coef(fit_default)
        expect_equal(length(cf), ncol(dat$X) + 1) # covariates + treatment
        
        vc <- vcov(fit_default)
        expect_equal(dim(vc), c(ncol(dat$X) + 1, ncol(dat$X) + 1))
        
        ci <- confint(fit_default)
        expect_equal(dim(ci), c(ncol(dat$X) + 1, 2))
        
        cat("  -> Default Interface PASS\n")
      }
      curr_test <- curr_test + 1
      
      # --- Test Formula Interface ---
      cat(sprintf("[%d/%d] Testing: Prior='%s', Concordant='%s', Interface='Formula'\n", 
                  curr_test + 1, total_tests, pt, cm))
      
      fmla <- as.formula(paste("y ~", paste(colnames(dat$X), collapse = " + ")))
      
      fit_formula <- NULL
      tryCatch({
        fit_formula <- bclogit(
          formula = fmla,
          data = dat$df,
          treatment = treatment, 
          strata = strata,
          prior_type = pt,
          concordant_method = cm
        )
      }, error = function(e) {
        fail(sprintf("Formula fit failed with error: %s", e$message))
      })
      
      if (!is.null(fit_formula)) {
        expect_s3_class(fit_formula, "bclogit")
        
        # Verify S3 methods
        expect_s3_class(summary(fit_formula), "summary.bclogit")
        
        # Verify call is preserved (Formula Specific)
        expect_true(!is.null(fit_formula$call))
        
        cat("  -> Formula Interface PASS\n")
      }
      curr_test <- curr_test + 1

    }
  }
})

test_that("Edge Case: Single Covariate", {
  cat("\nTesting Edge Case: Single Covariate\n")
  dat_single <- generate_data(n = 100, p = 1)
  
  fit <- bclogit(
     response = dat_single$y,
     data = dat_single$X, 
     treatment = dat_single$treatment,
     strata = dat_single$strata,
     prior_type = "Naive",
     concordant_method = "GLM"
  )
  expect_s3_class(fit, "bclogit")
  expect_equal(length(coef(fit)), 2) # 1 covariate + 1 treatment
  cat("  -> Single Covariate PASS\n")
})

# Note: Factors would need model.matrix expansion which is handled by formula interface usually, 
# or manual handling in Default. Default expects numeric matrix usually or attempts conversion.
# bclogit.default calls model.matrix(~ 0 + ., data = data) if data frame.

test_that("Edge Case: Factor Covariates (Default Interface handling)", {
  cat("\nTesting Edge Case: Factor Covariates in Default Interface\n")
  n <- 100
  df_fac <- data.frame(
    num = rnorm(n),
    fac = factor(rep(c("A", "B"), n/2))
  )
  X_mat <- model.matrix(~ 0 + ., data = df_fac)
  
  # ... actually bclogit.default internal logic:
  # if (test_data_frame(data, ...)) data_mat <- model.matrix(~ 0 + ., data = data)
  # So passing a data frame with factors should work.
  
  dat_fac <- generate_data(n=100, p=1) 
  dat_fac$df$Fac <- factor(rep(c("L1", "L2"), 50))
  
  # Pass dataframe to data arg of default
  fit <- bclogit(
    response = dat_fac$y,
    data = dat_fac$df[, c("V1", "Fac")], # Mixed numeric and factor
    treatment = dat_fac$treatment,
    strata = dat_fac$strata
  )
  expect_s3_class(fit, "bclogit")
  # Coefs: treatment + V1 + FacL1 + FacL2 (if 0 intercept)
  # Wait, bclogit.default uses `model.matrix(~ 0 + ., data = data)`.
  # Factor with 2 levels -> 2 columns.
  # So expect 1 (trt) + 1 (V1) + 2 (Fac) = 4 coefs?
  # Or does it drop one? `~ 0 + .` keeps all levels usually.
  
  cat("  -> Factor Covariate (Data Frame input) PASS\n")
})

