#' @param y A binary (0,1) vector containing the response of each subject.
#' @param X A data.frame, data.table, or model.matrix containing the variables.
#' @describeIn bclogit Default method for matrix/data input.
#' @export
bclogit.default <- function(y = NULL,
                            X = NULL,
                            treatment = NULL,
                            strata = NULL,
                            concordant_method = "GLM",
                            prior_type = "Naive",
                            treatment_name = NULL,
                            call = NULL,
                            chains = 4) {
  # ------------------------------------------------------------------------
  # 1. Input Validation and Pre-processing
  # ------------------------------------------------------------------------
  if (missing(treatment)) stop("The 'treatment' argument is required.")

  # Use copy of X if needed or just use X directly
  if (test_data_frame(X, types = c("numeric", "integer", "factor", "logical"))) {
    data_mat <- model.matrix(~ 0 + ., data = X)
  } else if (test_matrix(X, mode = "numeric")) {
    if (sd(X[, 1]) == 0) {
      X[, 1] <- NULL
    }
    data_mat <- as.matrix(X)
  } else {
    assert(
      check_data_frame(X),
      check_matrix(X, mode = "numeric")
    )
  }

  n <- nrow(data_mat)

  assertString(concordant_method, c("GLM", "GEE", "GLMM"))
  assertNumeric(y, lower = 0, upper = 1, any.missing = FALSE, len = n)
  assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)
  assertNumeric(strata, any.missing = FALSE, len = n)

  if (!all(treatment %in% c(0, 1))) {
    stop("Treatment must be binary 0 or 1.")
  }
  if (is.null(treatment_name)) {
    treatment_name <- deparse(substitute(treatment))
  }

  # Capture terms for summary/printing usage if possible (though we used matrix for fitting)
  model_terms <- tryCatch(terms(y ~ ., data = as.data.frame(data_mat)), error = function(e) NULL)

  # ------------------------------------------------------------------------
  # 2. Data Preparation
  # ------------------------------------------------------------------------

  X_concordant <- NULL
  y_concordant <- NULL
  treatment_concordant <- NULL

  X_diffs_discordant <- NULL
  y_diffs_discordant <- NULL
  treatment_diffs_discordant <- NULL

  # Internal variables for processing
  X <- data_mat
  X_model_matrix_col_names <- colnames(X)
  if (is.null(X_model_matrix_col_names)) {
    X_model_matrix_col_names <- paste0("X", 1:ncol(X))
  }

  if (length(y) <= ncol(X) + 5) {
    stop("Not enough rows. Must be at least the number of covariates plus 5")
  }


  matched_data <- process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )

  X_concordant <- matched_data$X_concordant
  y_concordant <- matched_data$y_concordant
  treatment_concordant <- matched_data$treatment_concordant
  strata_concordant <- matched_data$strata_concordant
  X_diffs_discordant <- matched_data$X_diffs_discordant
  y_diffs_discordant <- matched_data$y_diffs_discordant
  treatment_diffs_discordant <- matched_data$treatment_diffs_discordant

  # Ensure column names match for model formula consistency
  if (!is.null(X_concordant) && ncol(X_concordant) == length(X_model_matrix_col_names)) {
    colnames(X_concordant) <- X_model_matrix_col_names
  }

  # ------------------------------------------------------------------------
  # 3. Model Fitting
  # ------------------------------------------------------------------------


  # --- Concordant Pairs / Reservoir Model ---

  concordant_model <- NULL

  # Check for concordant pairs
  if (length(y_concordant) < ncol(X_concordant) + 5) {
    # Not enough data for concordant model
    warning("There are not enough concordant pairs or reservoir entries. The prior for the discordant pairs will be non-informative.")
    b_con <- rep(0, ncol(X_concordant) + 1)
    Sigma_con <- diag(101, ncol(X_concordant) + 1)
  } else {
    if (concordant_method == "GLM") {
      concordant_model <- glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
    }
    if (concordant_method == "GEE") {
      concordant_model <- geepack::geeglm(
        y_concordant ~ treatment_concordant + X_concordant,
        id = strata_concordant,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
      )
    }
    if (concordant_method == "GLMM") {
      concordant_model <- glmmTMB::glmmTMB(
        y_concordant ~ treatment_concordant + X_concordant + (1 | strata_concordant),
        family = binomial(),
        data   = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
      )
    }

    # Standardize prior extraction
    if (!is.null(concordant_model)) {
      # Get full coefficients including NAs for aliased terms
      if (concordant_method == "GLM") {
        full_b_con <- coef(concordant_model)
        full_Sigma_con <- vcov(concordant_model)
      } else if (concordant_method == "GEE") {
        full_b_con <- coef(concordant_model)
        full_Sigma_con <- vcov(concordant_model)
      } else if (concordant_method == "GLMM") {
        full_b_con <- glmmTMB::fixef(concordant_model)$cond
        full_Sigma_con <- vcov(concordant_model)$cond
      }

      # Identify target names: treatment + X columns
      target_names <- c("treatment_concordant", paste0("X_concordant", X_model_matrix_col_names))

      # Construct target vector
      K_stan <- ncol(X_diffs_discordant) + 1
      b_con <- numeric(K_stan)
      Sigma_con <- diag(100, K_stan)

      # Match coefficients
      for (i in seq_along(target_names)) {
        t_name <- target_names[i]
        if (t_name %in% names(full_b_con)) {
          val <- full_b_con[t_name]
          if (!is.na(val)) {
            b_con[i] <- val
          }
        }
      }

      # Match covariance
      valid_indices <- which(names(full_b_con) %in% target_names)
      valid_names <- names(full_b_con)[valid_indices]

      if (length(valid_names) > 0) {
        vcov_names <- rownames(full_Sigma_con)
        for (r_name in valid_names) {
          if (r_name %in% vcov_names) {
            target_idx_r <- match(r_name, target_names)
            for (c_name in valid_names) {
              if (c_name %in% vcov_names) {
                target_idx_c <- match(c_name, target_names)
                val <- full_Sigma_con[r_name, c_name]
                if (!is.na(val)) {
                  Sigma_con[target_idx_r, target_idx_c] <- val
                }
              }
            }
          }
        }
      }

      # Override treatment prior
      b_con[1] <- 0
      Sigma_con[1, ] <- 0
      Sigma_con[, 1] <- 0
      Sigma_con[1, 1] <- 100
    }
  }

  # Sanitize priors to ensure no NA values are passed
  if (exists("b_con") && any(is.na(b_con))) {
    b_con[is.na(b_con)] <- 0
  }
  if (exists("Sigma_con")) {
    items_fixed <- FALSE
    if (any(is.na(Sigma_con))) {
      warning("Prior covariance contains NAs. Replacing with diffuse prior.")
      Sigma_con <- diag(100, nrow(Sigma_con))
      b_con <- rep(0, length(b_con))
      items_fixed <- TRUE
    }

    if (!items_fixed) {
      # Check positive definiteness
      is_pd <- tryCatch(
        {
          chol(Sigma_con)
          TRUE
        },
        error = function(e) FALSE
      )

      if (!is_pd) {
        warning("Prior covariance is not positive definite. Replacing with diffuse prior.")
        Sigma_con <- diag(100, nrow(Sigma_con))
        b_con <- rep(0, length(b_con))
      }
    }
  }

  discordant_model <- NULL
  converged_discordant <- NULL

  # --- Discordant Pairs Model ---
  if (length(y_diffs_discordant) < ncol(X_diffs_discordant) + 5) {
    stop("There are not enough discordant pairs. The model will not be fit .")
  } else {
    Sigma_con <- (Sigma_con + t(Sigma_con)) / 2
    y_dis_0_1 <- ifelse(y_diffs_discordant == -1, 0, 1)
    wX_dis <- cbind(treatment_diffs_discordant, X_diffs_discordant)

    stan_file <- switch(prior_type,
      "Naive" = "mvn_logistic.stan",
      "G prior" = "mvn_logistic_gprior.stan",
      "PMP" = "mvn_logistic_PMP.stan",
      "Hybrid" = "mvn_logistic_Hybrid.stan",
      stop("Unknown prior type")
    )

    # Get compiled model (cached)
    stan_mod <- get_stan_model(stan_file)

    if (prior_type %in% c("Naive", "G prior")) {
      data_list <- list(
        N = nrow(wX_dis),
        K = ncol(wX_dis),
        X = wX_dis,
        y = y_dis_0_1,
        mu = b_con,
        Sigma = Sigma_con
      )
    }
    if (prior_type %in% c("PMP", "Hybrid")) {
      # PMP and Hybrid share this setup
      proj_matrix <- X_diffs_discordant %*% solve(t(X_diffs_discordant) %*% X_diffs_discordant) %*% t(X_diffs_discordant)
      w_dis_ortho <- treatment_diffs_discordant - proj_matrix %*% treatment_diffs_discordant

      data_list <- list(
        N = nrow(X_diffs_discordant),
        P = ncol(X_diffs_discordant),
        y = y_dis_0_1,
        xw = as.vector(w_dis_ortho),
        X = X_diffs_discordant,
        mu_A = as.array(b_con[-1]),
        Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
      )
    }

    # Sample from the model
    discordant_model <- tryCatch(
      {
        rstan::sampling(stan_mod, data = data_list, refresh = 0, chains = chains)
      },
      error = function(e) {
        warning(sprintf("Sampling failed: %s", e$message))
        NULL
      }
    )
  }

  # ------------------------------------------------------------------------
  # 4. Result Construction
  # ------------------------------------------------------------------------

  # Extract coefficients from the Stan model if available
  coefficients <- NULL
  var_cov <- NULL

  if (!is.null(discordant_model)) {
    # Extract posterior samples

    sims <- rstan::extract(discordant_model)

    if (prior_type %in% c("Naive", "G prior")) {
      beta_post <- sims$beta
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    } else {
      # PMP / Hybrid
      beta_w_post <- as.vector(sims$beta_w)
      beta_nuis_post <- sims$beta_nuis

      beta_post <- cbind(beta_w_post, beta_nuis_post)
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    }

    # Calculate posterior means
    coefficients <- colMeans(beta_post)
    names(coefficients) <- col_names_final

    # Calculate posterior covariance
    var_cov <- cov(beta_post)
    # Assign names
    rownames(var_cov) <- col_names_final
    colnames(var_cov) <- col_names_final
  }

  # Construct result object
  res <- list(
    coefficients = coefficients,
    var = var_cov,
    model = discordant_model,
    concordant_model = concordant_model,
    prior_info = list(
      mu = if (exists("b_con")) b_con else NULL,
      Sigma = if (exists("Sigma_con")) Sigma_con else NULL
    ),

    # Standard S3 components
    call = if (!is.null(call)) call else match.call(),
    terms = model_terms,
    xlevels = NULL,

    # Metadata for summary
    n = n,
    num_discordant = length(y_diffs_discordant),
    num_concordant = length(y_concordant) / 2,
    X_model_matrix_col_names = X_model_matrix_col_names,
    treatment_name = treatment_name
  )

  class(res) <- c("bclogit", "list")
  return(res)
}

#' Helper to get compiled stan model from cache
#' @keywords internal
get_stan_model <- function(file_name) {
  # Check if model exists in globals
  if (!exists(file_name, envir = bclogit_globals)) {
    message(paste("Compiling", file_name, "..."))
    stan_file_path <- system.file(file.path("stan", file_name), package = "bclogit")
    if (stan_file_path == "") {
      stop(paste("Stan file not found:", file_name))
    }

    # Compile and assign to globals
    mod <- rstan::stan_model(file = stan_file_path)
    assign(file_name, mod, envir = bclogit_globals)
  }

  return(get(file_name, envir = bclogit_globals))
}
