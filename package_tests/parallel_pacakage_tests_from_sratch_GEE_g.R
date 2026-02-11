pacman::p_load(dplyr, tidyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2, geepack, glmmTMB, rstan, binaryMM, rstanarm) # doParallel

if (!require("bclogit", character.only = TRUE)) {
  remotes::install_local("bclogit", dependencies = FALSE, force = TRUE, upgrade = "never")
  library(bclogit)
}

num_cores <- availableCores() - 6
ps <- c(6)
Nsim <- 100
external_nsim <- 100000
ns <- c(100, 250, 500)
beta_Ts <- c(0, 0.5)
X_styles <- c("non-correlated") # "correlated",
true_funtions <- c("linear", "non-linear")
regress_on_Xs <- c("one", "some") # "all", "one", "none"

# Bayesian_prior_for_betaTs = c(TRUE, FALSE)
#sm <- rstan::stan_model("mvn_logistic.stan")
sm_g <- rstan::stan_model("mvn_logistic_gprior.stan")
#sm_PMP <- rstan::stan_model("mvn_logistic_PMP.stan")
sm_hybrid <- rstan::stan_model("mvn_logistic_Hybrid.stan")

params <- expand.grid(
  nsim = 1:Nsim,
  p = ps,
  beta_T = beta_Ts,
  n = ns,
  X_style = X_styles
)

params <- params %>%
  arrange(nsim, p, beta_T, n, X_style)

############### CODE ###############

Bayesian_Clogit <- function(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit) {
  if (concordant_fit == "GLM") {
    fit_con <- glm(y_con ~ w_con + X_con, family = "binomial")
    b_con <- summary(fit_con)$coefficients[, 1]
    Sigma_con <- pmin(vcov(fit_con), 1e3)
    eps <- 1e-6 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

    b_con <- c(0, b_con[-c(1,2)])
    Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 1e3
  }
  if (concordant_fit == "GEE") {
    strat_con <- rep(1:(nrow(X_con) / 2), each = 2)
    
    warning_con = NA
    
    fit_con <- withCallingHandlers({
      geeglm(
        y_con ~ w_con + X_con,
        id = strat_con,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data = data.frame(y_con, w_con, X_con, strat_con)
      )
    }, warning = function(w) {
      warning_con <<- w$message
      invokeRestart("muffleWarning") 
    })

    b_con <- summary(fit_con)$coefficients[, 1]
    
    b_con_max_abs = max(abs(b_con[-1]))
    Sigma_con_max_abs = max(abs(vcov(fit_con)[-1,-1]))
    Sigma_con_min_abs = min(abs(vcov(fit_con)[-1,-1]))
    
    Sigma_con <- pmin(vcov(fit_con), 1e3)
    eps <- 1e-10 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

    b_con <- c(0, b_con[-c(1,2)])
    Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 1e3
    
    mtv <- as.numeric(summary(fit_con)$corr["Estimate"])
    pval_mtv <- 2 * pnorm(as.numeric(abs(mtv / summary(fit_con)$corr["Std.err"])), lower.tail = FALSE)
  }
  if (concordant_fit == "GLMM") {
    strat_con <- rep(1:(nrow(X_con) / 2), each = 2)

    fit_con <- glmmTMB(
      y_con ~  w_con + X_con + (1 | strat_con),
      family = binomial(),
      data   = data.frame(y_con, w_con, X_con, strat_con)
    )

    b_con <- summary(fit_con)$coefficients$cond[, 1]
    Sigma_con <- pmin(vcov(fit_con)$cond, 1e3)
    eps <- 1e-6 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

    b_con <- c(0, b_con[-c(1,2)])
    Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 1e3
  }

  Sigma_con = (Sigma_con + t(Sigma_con)) / 2
  y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
  wX_dis <- cbind(w_dis, X_dis)
  g = NA
  
  if (prior_type == "normal") {
    if (all(diag(Sigma_con) == 1e3)) {
      ret <- list()
      ret$betaT <- NA
      ret$ssq_beta_T <- NA
      ret$reject <- NA
      ret$pval <- NA
      return(ret) # model blew up
    }

    data_list <- list(
      N = nrow(wX_dis),
      K = ncol(wX_dis),
      X = wX_dis,
      y = y_dis_0_1,
      mu = b_con, # example prior mean
      Sigma = Sigma_con # example covariance (wide prior)
    )

    discordant_model <- tryCatch(
      {
        # for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
        summary(rstan::sampling(sm, data = data_list, refresh = 0, chains = 1))$summary["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "G prior") {
    would_have_canceled = all(diag(Sigma_con) == 1e3)
    # if (all(diag(Sigma_con) == 1e3)) {
    #   ret <- list()
    #   ret$betaT <- NA
    #   ret$ssq_beta_T <- NA
    #   ret$reject <- NA
    #   ret$pval <- NA
    #   return(ret) # model blew up
    # }

    data_list <- list(
      N = nrow(wX_dis),
      K = ncol(wX_dis),
      X = wX_dis,
      y = y_dis_0_1,
      mu = b_con,
      Sigma = Sigma_con
    )

    failior_dis = NA
    discordant_model <- tryCatch(
      {
        # for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
        mod = summary(rstan::sampling(sm_g, data = data_list, refresh = 0, chains = 1))$summary
        g = mod["g", "mean"]
        mod["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        failior_dis <<- e$message
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "PMP") {
    proj_matrix <- X_dis %*% solve(t(X_dis) %*% X_dis) %*% t(X_dis)
    w_dis_ortho <- w_dis - proj_matrix %*% w_dis

    # Prepare the data list for Stan (matching mvn_logistic_PMP.stan)
    data_list <- list(
      N = nrow(X_dis),
      P = ncol(X_dis),
      y = y_dis_0_1,
      xw = as.vector(w_dis_ortho),
      X = X_dis,
      mu_A = as.array(b_con[-1]),
      Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
    )

    discordant_model <- tryCatch(
      {
        # for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
        summary(rstan::sampling(sm_PMP, data = data_list, refresh = 0, chains = 1))$summary["beta_w", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "Hybrid") {
    would_have_canceled = all(diag(Sigma_con) == 1e3)
    proj_matrix <- X_dis %*% solve(t(X_dis) %*% X_dis) %*% t(X_dis)
    w_dis_ortho <- w_dis - proj_matrix %*% w_dis


    # Prepare the data list for Stan
    data_list <- list(
      N = nrow(X_dis),
      P = ncol(X_dis),
      y = y_dis_0_1,
      xw = as.vector(w_dis_ortho),
      X = X_dis,
      mu_A = as.array(b_con[-1]),
      Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
    )

    failior_dis = NA
    discordant_model <- tryCatch(
      {
        mod = summary(rstan::sampling(sm_hybrid, data = data_list, refresh = 0, chains = 1))$summary
        g = mod["g", "mean"]
        mod["beta_w", c("mean", "sd", "2.5%", "97.5%")]
        
        #summary(rstan::sampling(sm_hybrid, data = data_list, refresh = 0, chains = 1))$summary["beta_w", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        failior_dis <<- e$message
        warning(sprintf("Hybrid PMP failed: %s", e$message))
        NULL
      }
    )
  }
  
  ret <- list(
    betaT = NA,
    ssq_beta_T = NA,
    reject = NA,
    pval = NA,
    g = NA,
    mtv = NA,
    pval_mtv = NA,
    b_con_max_abs = NA,
    Sigma_con_max_abs = NA,
    Sigma_con_min_abs = NA,
    would_have_canceled = NA,
    warning_con = NA,
    failior_dis = NA
  )

  if (!is.null(discordant_model)) {
    ret$betaT <- discordant_model[1]
    ret$ssq_beta_T <- discordant_model[2]
    ret$reject <- !(discordant_model[3] < 0 & discordant_model[4] > 0)
    ret$pval <- if (ret$reject) 0 else 1
    ret$g = g
    ret$mtv = mtv
    ret$pval_mtv = pval_mtv
    ret$b_con_max_abs = b_con_max_abs
    ret$Sigma_con_max_abs = Sigma_con_max_abs
    ret$Sigma_con_min_abs = Sigma_con_min_abs
    ret$would_have_canceled <- would_have_canceled
    ret$warning_con <- warning_con
    ret$failior_dis <- failior_dis
  }
  

  return(ret)
}

Do_Inference <- function(y, X, w, strat, p, beta_T, n, X_style, true_funtion, regress_on_X) {
  res <- data.frame(
    n = numeric(),
    p = numeric(),
    beta_T = numeric(),
    X_style = character(),
    true_funtion = character(),
    regress_on_X = character(),
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric(),
    mtv = numeric(),
    pval_mtv = numeric(),
    g = numeric(),
    mean_y = numeric(),
    mean_y_con = numeric(),
    n_con = numeric(),
    n_dis = numeric(),
    percent_dis = numeric(),
    b_con_max_abs = numeric(),
    Sigma_con_max_abs = numeric(),
    Sigma_con_min_abs = numeric(),
    would_have_canceled = logical(),
    warning_con = character(),
    failior_dis = character()
  )
  
  
  matched_data <-
    bclogit:::process_matched_pairs_cpp(
      strata = strat,
      y = y,
      X = data.matrix(X),
      treatment = w
    )

  X_con <- matched_data$X_concordant
  y_con <- matched_data$y_concordant
  w_con <- matched_data$treatment_concordant
  X_dis <- matched_data$X_diffs_discordant
  y_dis <- matched_data$y_diffs_discordant
  w_dis <- matched_data$treatment_diffs_discordant
  dis_idx <- matched_data$discordant_idx + 1

  if (regress_on_X %in% c("one", "some")) { # "all",
    discordant_viabele <- if (length(y_dis) > ncol(X) + 7) {
      TRUE
    } else {
      FALSE
    }
    concordant_viabele <- if (length(y_con) > ncol(X) + 7) {
      TRUE
    } else {
      FALSE
    }

    ########################### CLOGIT  ###########################
    # beta_hat_T <- NA
    # ssq_beta_hat_T <- NA
    # pval <- NA
    # modle_tracking_value <- NA
    # pval_modle_tracking_value <- NA
    # if (discordant_viabele) {
    #   y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
    #   model <- summary(glm(y_dis_0_1 ~ 0 + w_dis + X_dis, family = "binomial"))$coefficients[1, c(1, 2)]
    #   beta_hat_T <- model[1]
    #   ssq_beta_hat_T <- model[2]
    #   pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    # }
    # res <- rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "clogit",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval,
    #   modle_tracking_value = modle_tracking_value,
    #   pval_modle_tracking_value = pval_modle_tracking_value,
    #   g = NA
    # ))

    ########################### LOGIT  ###########################
    # beta_hat_T <- NA
    # ssq_beta_hat_T <- NA
    # pval <- NA
    # modle_tracking_value <- NA
    # pval_modle_tracking_value <- NA
    # if (TRUE) {
    #   model <- summary(glm(y ~ w + X, family = "binomial"))$coefficients[2, c(1, 2)]
    #   beta_hat_T <- model[1]
    #   ssq_beta_hat_T <- model[2]
    #   pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    # }
    # res <- rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "GLM",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval,
    #   modle_tracking_value = modle_tracking_value,
    #   pval_modle_tracking_value = pval_modle_tracking_value,
    #   g = NA
    # ))
    
    ########################### BAYESIAN ###########################
    for (prior_type in c("G prior", "Hybrid")) { #"normal", "G prior", "PMP", "Hybrid"
      for (concordant_fit in c("GEE")) { #"GLM", "GEE", "GLMM"
        beta_hat_T <- NA ; ssq_beta_hat_T <- NA ; pval <- NA ; beta_hat_T <- NA
        ssq_beta_hat_T <- NA ; pval <- NA ; mtv <- NA ; pval_mtv <- NA ; g <- NA
        mean_y <- mean(y)
        mean_y_con <- mean(y_con)
        n_con <- length(y_con)
        n_dis <- n - n_con
        percent_dis <- n_dis / n
        
        if (discordant_viabele & concordant_viabele) {
          model <- Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit)
          beta_hat_T <- model$betaT
          ssq_beta_hat_T <- model$ssq_beta_T
          pval <- model$pval
          g <- model$g
          mtv <- model$mtv
          pval_mtv <- model$pval_mtv
          b_con_max_abs <- model$b_con_max_abs
          Sigma_con_max_abs <- model$Sigma_con_max_abs
          Sigma_con_min_abs <- model$Sigma_con_min_abs
          would_have_canceled <- model$would_have_canceled
          warning_con <- model$warning_con
          failior_dis <- model$failior_dis
        }
        res <- rbind(res, data.frame(
          n = n,
          p = p,
          beta_T = beta_T,
          X_style = X_style,
          true_funtion = true_funtion,
          regress_on_X = regress_on_X,
          inference = paste0("GEE + ", prior_type),
          beta_hat_T = beta_hat_T,
          pval = pval,
          mtv = mtv,
          pval_mtv = pval_mtv,
          g = g,
          mean_y = mean_y,
          mean_y_con = mean_y_con,
          n_con = n_con,
          n_dis = n_dis,
          percent_dis = percent_dis,
          b_con_max_abs = b_con_max_abs,
          Sigma_con_max_abs = Sigma_con_max_abs,
          Sigma_con_min_abs = Sigma_con_min_abs,
          would_have_canceled = would_have_canceled,
          warning_con = warning_con,
          failior_dis = failior_dis
        ))
      }
    }

    ########################### glmmTMB  ###########################
    # beta_hat_T <- NA
    # ssq_beta_hat_T <- NA
    # pval <- NA
    # modle_tracking_value <- NA
    # pval_modle_tracking_value <- NA
    # tryCatch(
    #   {
    #     fit_tmb <- glmmTMB(
    #       y ~ X + w + (1 | strat),
    #       family = binomial(),
    #       data   = data.frame(y, X, w, strat)
    #     )
    #     model <- summary(fit_tmb)$coefficients$cond["w", c("Estimate", "Std. Error")]
    #     beta_hat_T <- model[1]
    #     ssq_beta_hat_T <- model[2]
    #     pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    #     modle_tracking_value = VarCorr(fit_tmb)$cond$strat[1,1]
    #     
    #     m_null <- glmmTMB(
    #       y ~ X + w,
    #       family = binomial(),
    #       data   = data.frame(y, X, w, strat)
    #     )
    #     pval_modle_tracking_value <- anova(m_null, fit_tmb)[2, "Pr(>Chisq)"]
    #     
    #   },
    #   error = function(e) {
    #     beta_hat_T <<- NA
    #     ssq_beta_hat_T <<- NA
    #     pval <<- NA
    #   }
    # )
    # res <- rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "GLMM",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval,
    #   modle_tracking_value = modle_tracking_value,
    #   pval_modle_tracking_value = pval_modle_tracking_value,
    #   g = NA
    # ))

    ########################### GEE  ###########################
    # beta_hat_T <- NA
    # ssq_beta_hat_T <- NA
    # pval <- NA
    # modle_tracking_value <- NA
    # pval_modle_tracking_value <- NA
    # tryCatch(
    #   {
    #     fit_gee <- geeglm(
    #       y ~ X + w,
    #       id = strat,
    #       family = binomial(link = "logit"),
    #       corstr = "exchangeable",
    #       data = data.frame(y, X, w, strat)
    #     )
    #     model <- summary(fit_gee)$coefficients["w", c("Estimate", "Std.err")]
    #     beta_hat_T <- as.numeric(model[1])
    #     ssq_beta_hat_T <- as.numeric(model[2])
    #     pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    #     modle_tracking_value <- as.numeric(summary(fit_gee)$corr["Estimate"])
    #     pval_modle_tracking_value <- 2 * pnorm(as.numeric(abs(modle_tracking_value / summary(fit_gee)$corr["Std.err"])), lower.tail = FALSE)
    #     
    #   },
    #   error = function(e) {
    #     beta_hat_T <<- NA
    #     ssq_beta_hat_T <<- NA
    #     pval <<- NA
    #   }
    # )
    # res <- rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "GEE",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval,
    #   modle_tracking_value = modle_tracking_value,
    #   pval_modle_tracking_value = pval_modle_tracking_value,
    #   g = NA
    # ))

    ########################### STAN_GLMER ###########################
    # beta_hat_T <- NA
    # ssq_beta_hat_T <- NA
    # pval <- NA
    # tryCatch(
    #   {
    #     fit_stan <- stan_glmer(
    #       y ~ X + w + (1 | strat),
    #       family = binomial(link = "logit"),
    #       data = data.frame(y, X, w, strat),
    #       prior = normal(0, 2.5),
    #       prior_intercept = normal(0, 5),
    #       chains = 2,
    #       iter = 2000,
    #       refresh = 0
    #     )
    # 
    #     post <- as.matrix(fit_stan)
    #     w_draws <- post[, "w"]
    # 
    #     beta_hat_T <- mean(w_draws)
    #     ssq_beta_hat_T <- sd(w_draws)
    # 
    #     reject <- !(quantile(w_draws, 0.025) < 0 & quantile(w_draws, 0.975) > 0)
    #     pval <- if (reject) 0 else 1
    #   },
    #   error = function(e) {
    #     beta_hat_T <<- NA
    #     ssq_beta_hat_T <<- NA
    #     pval <<- NA
    #   }
    # )
    # res <- rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "stan_glmer",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval
    # ))

    ########################### binaryMM ###########################
    # beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    # tryCatch({
    #   # binaryMM requires a data frame and explicit column names
    #   df_mm = data.frame(y = y, w = w, strat = strat)
    #   df_mm = cbind(df_mm, X) # Add covariates
    #
    #   # Construct the formula dynamically based on X columns
    #   if (is.null(colnames(X))) {
    #     colnames(X) = paste0("V", 1:ncol(X))
    #   }
    #   x_names = colnames(X)
    #   mm_formula = as.formula(paste("y ~ w +", paste(x_names, collapse = " + ")))
    #
    #   fit_mm = mm(
    #     formula = mm_formula,
    #     id = strat,
    #     data = df_mm,
    #     # 'average' provides the population-averaged (marginal) effects
    #     # similar to GEE, but via a different estimation method
    #     type = "marginal"
    #   )
    #
    #   # Extracting coefficients for the treatment effect 'w'
    #   # binaryMM stores results in a 'coefficients' matrix
    #   mm_sum = summary(fit_mm)
    #   beta_hat_T = mm_sum$coefficients["w", "Estimate"]
    #   ssq_beta_hat_T = mm_sum$coefficients["w", "Std.Error"]
    #
    #   # Standard Wald-type p-value calculation
    #   pval = 2 * pnorm(-abs(beta_hat_T / ssq_beta_hat_T))
    #
    # }, error = function(e) {
    #   beta_hat_T <<- NA; ssq_beta_hat_T <<- NA; pval <<- NA
    # })
    #
    # res = rbind(res, data.frame(
    #   n = n,
    #   p = p,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "binaryMM",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval
    # ))
  } else { # if no x then remove the x parameter
    discordant_viabele <- if (length(y_dis) > 5) {
      TRUE
    } else {
      FALSE
    }

    ########################### CLOGIT  ###########################
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    if (discordant_viabele) {
      y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
      model <- summary(glm(y_dis_0_1 ~ 0 + w_dis, family = "binomial"))$coefficients[1, c(1, 2)]
      beta_hat_T <- model[1]
      ssq_beta_hat_T <- model[2]
      pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res <- rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "clogit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))

    ########################### LOGIT  ###########################
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    if (discordant_viabele) {
      model <- summary(glm(y ~ w, family = "binomial"))$coefficients[2, c(1, 2)]
      beta_hat_T <- model[1]
      ssq_beta_hat_T <- model[2]
      pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res <- rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "logit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))

    ########################### glmmTMB  ###########################
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    tryCatch(
      {
        fit_tmb <- glmmTMB(
          y ~ w + (1 | strat),
          family = binomial(),
          data   = data.frame(y, X, w, strat)
        )
        model <- summary(fit_tmb)$coefficients$cond["w", c("Estimate", "Std. Error")]
        beta_hat_T <- model[1]
        ssq_beta_hat_T <- model[2]
        pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
      },
      error = function(e) {
        beta_hat_T <<- NA
        ssq_beta_hat_T <<- NA
        pval <<- NA
      }
    )
    res <- rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "glmmTMB",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))

    ########################### GEE  ###########################
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    tryCatch(
      {
        fit_gee <- geeglm(
          y ~ w,
          id = strat,
          family = binomial(link = "logit"),
          corstr = "exchangeable",
          data = data.frame(y, X, w, strat)
        )
        model <- summary(fit_gee)$coefficients["w", c("Estimate", "Std.err")]
        beta_hat_T <- as.numeric(model[1])
        ssq_beta_hat_T <- as.numeric(model[2])
        pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
      },
      error = function(e) {
        beta_hat_T <<- NA
        ssq_beta_hat_T <<- NA
        pval <<- NA
      }
    )
    res <- rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "GEE",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
  }

  rownames(res) <- NULL
  return(res)
}

Run_sim <- function(p, beta_T, n, X_style) {
  BIG_res <- data.frame()
  y <- array(NA, n)
  probs <- array(NA, n)

  if (X_style == "correlated") { ############## correlated is the old style
    Sigma <- 1 * (matrix(0.5, nrow = p, ncol = p) + diag(1 - 0.5, p))
    X <- MASS::mvrnorm(n / 2, rep(0, p), Sigma)
    X <- pnorm(X)
    X <- matrix(2 * X - 1, ncol = p)
  } else {
    X <- matrix(runif((n / 2) * p, min = -1, max = 1), ncol = p)
    X_plus_eps <- X + matrix(rnorm((n / 2) * p, 0, 0.05), ncol = p)
    combined <- rbind(X, X_plus_eps)
    ids <- order(c(1:(n / 2), 1:(n / 2)))
    X <- combined[ids, ]
    X[, 1] <- runif(n, min = -1, max = 1)
    rm(X_plus_eps, combined, ids)
  }
  w <- c(rbind(replicate(n / 2, sample(c(0, 1)), simplify = TRUE)))
  strat <- rep(1:(n / 2), each = 2)


  for (true_funtion in true_funtions) {
    if (true_funtion == "linear") {
      beta_X_value <- if (p == 6) {
        1.25
      } else {
        0.75
      }
      beta_X <- rep(beta_X_value, p)
      beta_0 <- -0.5
      probs <- 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_X) + beta_T * w)))
    } else {
      f_x <- sin(pi * X[, 1] * X[, 2]) + X[, 3]^3 + X[, 4]^2 + X[, 5]^2
      probs <- 1 / (1 + exp(-(f_x + beta_T * w)))
    }
    y <- rbinom(n, 1, probs)

    # df = data.frame(y = y, probs = probs)
    # ggplot(df, aes(x = probs, fill = factor(y))) +
    #   geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    #   labs(x = "Predicted probability", fill = "Outcome") +
    #   scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
    #   theme_minimal()

    for (regress_on_X in regress_on_Xs) {
      if (regress_on_X == "one") {
        X_run <- X[, 1, drop = FALSE]
      } else if (regress_on_X == "some") {
        if (p == 20) {
          X_run <- X[, 1:5, drop = FALSE]
        } else {
          X_run <- X[, c(1, 2), drop = FALSE]
        }
      } else {
        X_run <- X
      }
      one_res <- Do_Inference(y, X_run, w, strat, p, beta_T, n, X_style, true_funtion, regress_on_X)
      BIG_res <- rbind(BIG_res, one_res)
    }
  }
  return(BIG_res)
}

results_2 = data.frame(
  n = numeric(),
  p = numeric(),
  beta_T = numeric(),
  X_style = character(),
  true_funtion = character(),
  regress_on_X = character(),
  inference = character(),
  beta_hat_T = numeric(),
  pval = numeric(),
  mtv = numeric(),
  pval_mtv = numeric(),
  g = numeric(),
  mean_y = numeric(),
  mean_y_con = numeric(),
  n_con = numeric(),
  n_dis = numeric(),
  percent_dis = numeric(),
  b_con_max_abs = numeric(),
  Sigma_con_max_abs = numeric(),
  Sigma_con_min_abs = numeric(),
  would_have_canceled = logical(),
  warning_con = character(),
  failior_dis = character()
)
for (j in 1:120) {
  cat("################", j, "################\n")
  p = params[j,]$p
  beta_T = params[j,]$beta_T
  n = params[j,]$n
  X_style = params[j,]$X_style
  results_2 = rbind(results_2, Run_sim(p = p, beta_T = beta_T, n = n, X_style = X_style))
}

############### SIM SET UP ###############

handlers(global = TRUE)
handlers("txtprogressbar")

registerDoFuture()
plan(multisession, workers = num_cores)

############### SIM ###############

for (e_nsim in 1:external_nsim) {
  with_progress({
    prog <- progressor(along = 1:nrow(params))

    results <- foreach(
      row = iter(params, by = "row"),
      .combine = rbind,
      .packages = c(
        "bclogit", "nbpMatching", "data.table",
        "dplyr", "MASS", "Rcpp", "rstanarm"
      )
    ) %dorng% {
      # extract parameters
      nsim <- row$nsim
      p <- row$p
      beta_T <- row$beta_T
      n <- row$n
      X_style <- row$X_style
      res <- tryCatch(
        {
          out <- Run_sim(p, beta_T, n, X_style)
          # cat("Successfully ran simulation")
          prog()
          out
        },
        error = function(e) {
          cat(glue::glue("Error in nsim={nsim}: {e$message}"), "\n")
          prog() # still update progress bar even if it fails
          NULL # return NULL if failed, will be dropped in rbind
        }
      )
    }
  })
  write.csv(results, file = paste0("C:/temp/clogitR_kap_test_from_scratch/", Nsim, "_", e_nsim, "_new.csv"), row.names = FALSE)
  rm(results)
  gc()
}

plan(sequential)


############### COMPILE RESULTS ###############

results <- read.csv("C:/temp/clogitR_kap_test_from_scratch/100_1_new.csv")

sum <- 1
for (i in 2:27) {
  file_path <- paste0("C:/temp/clogitR_kap_test_from_scratch/100_", i, "_new.csv")
  if (file.exists(file_path)) {
    sum <- sum + 1
    message("Reading file ", i)
    temp <- read.csv(file_path)
    results <- rbind(results, temp)
  } else {
    message("Skipping missing file ", i)
  }
}

# for (i in 2:475){
#   print(i)
#   results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test_from_scratch/1000_", i, ".csv")))
# }

results$X <- NULL

res_mod <- results %>%
  mutate(
    lower_ci = beta_hat_T - (1.96 * ssq_beta_hat_T),
    upper_ci = beta_hat_T + (1.96 * ssq_beta_hat_T),
    covered = (lower_ci <= beta_T) & (upper_ci >= beta_T),
    sq_err = (beta_hat_T - beta_T)^2,
    rej = pval < 0.05,
    rej_mtv = pval_modle_tracking_value < 0.05
  ) %>%
  group_by(p, beta_T, true_funtion, regress_on_X, n, inference, X_style) %>%
  summarize(
    num_na = sum(is.na(pval)),
    num_real = sum(!is.na(pval)),
    mse = mean(sq_err, na.rm = TRUE, trim = 0.001),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    coverage = mean(covered, na.rm = TRUE),
    mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_beta_hat_T = mean(ssq_beta_hat_T, trim = 0.001, na.rm = TRUE),
    significant = binom.test(sum(rej, na.rm = TRUE), n = (n() - num_na),  p = 0.05)$p.value,
    mean_mtv = mean(modle_tracking_value, trim = 0.001, na.rm = TRUE),
    num_na_mtv = sum(is.na(pval_modle_tracking_value)),
    percent_reject_mtv = sum(rej_mtv, na.rm = TRUE) / (n() - num_na_mtv),
    significant_mtv = ifelse(num_na_mtv < 0.9 * n(), binom.test(sum(rej_mtv, na.rm = TRUE), n = (n() - num_na_mtv),  p = 0.05)$p.value, NaN),
    mean_g = mean(g, na.rm = TRUE),
    .groups = "drop"
  )


write.csv(res_mod, file = "C:/temp/clogitR_kap_test_from_scratch/combined_2700_3.csv", row.names = FALSE)



table(res_bay$inference)

res_bay = res_mod[startsWith(res_mod$inference, "bayes") & (res_mod$beta_T != 0), ]

res_bay = res_bay %>% mutate(
  prior = sub("^[^_]*_([^_]*)_.*$", "\\1", inference),
  method = sub(".*_", "", inference),
  log_pow = log(percent_reject),
  prior = factor(prior, levels = c("normal", "G prior", "PMP", "Hybrid")),
  method = factor(method, levels = c("GLM", "GEE", "GLMM")),
  n = factor(n)
) %>%
  arrange(prior, method)

summary(lm(log_pow ~ n*method + regress_on_X + true_funtion + prior ,data = res_bay))


ggplot(data = data.frame(d = log(results[results$inference == "bayesian_G prior_GEE",]$beta_hat_T^2)), aes(x = d)) +
  geom_histogram()
