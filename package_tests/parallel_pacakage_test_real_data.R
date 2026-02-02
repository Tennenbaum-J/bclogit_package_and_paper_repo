pacman::p_load(dplyr, tidyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2, geepack, glmmTMB, rstan, binaryMM, rstanarm) # doParallel

if (!require("bclogit", character.only = TRUE)) {
  remotes::install_local("bclogit", dependencies = FALSE, force = TRUE, upgrade = "never")
  library(bclogit)
}

data(diabetic, package="survival")



num_cores <- availableCores() - 6
Nsim <- 100
external_nsim <- 100000
ns <- c(100, 250, 394)

# Bayesian_prior_for_betaTs = c(TRUE, FALSE)
sm <- rstan::stan_model("mvn_logistic.stan")
sm_g <- rstan::stan_model("mvn_logistic_gprior.stan")
sm_PMP <- rstan::stan_model("mvn_logistic_PMP.stan")
sm_hybrid <- rstan::stan_model("mvn_logistic_Hybrid.stan")

params <- expand.grid(
  nsim = 1:Nsim,
  n = ns
)

params <- params %>%
  arrange(nsim, n)

############### CODE ###############

Bayesian_Clogit <- function(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit) {
  if (concordant_fit == "GLM") {
    fit_con <- glm(y_con ~ X_con, family = "binomial")
    b_con <- summary(fit_con)$coefficients[, 1]
    Sigma_con <- pmin(vcov(fit_con), 20)
    eps <- 1e-6 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))
    
    b_con <- c(0, b_con[-c(1:3)])
    #Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 20
  } else if (concordant_fit == "GEE") {
    strat_con <- rep(1:(nrow(X_con) / 2), each = 2)
    fit_con <- geeglm(
      y_con ~  X_con,
      id = strat_con,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      data = data.frame(y_con, w_con, X_con, strat_con)
    )
    
    b_con <- summary(fit_con)$coefficients[, 1]
    Sigma_con <- pmin(vcov(fit_con), 20)
    eps <- 1e-6 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))
    
    b_con <- c(0, b_con[-c(1:3)])
    #Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 20
  } else if (concordant_fit == "GLMM") {
    strat_con <- rep(1:(nrow(X_con) / 2), each = 2)
    
    fit_con <- glmmTMB(
      y_con ~  X_con + (1 | strat_con),
      family = binomial(),
      data   = data.frame(y_con, w_con, X_con, strat_con)
    )
    
    b_con <- summary(fit_con)$coefficients$cond[, 1]
    Sigma_con <- pmin(vcov(fit_con)$cond, 20)
    eps <- 1e-6 # added for stabilibty
    Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))
    
    b_con <- c(0, b_con[-c(1:3)])
    #Sigma_con <- Sigma_con[-1, -1]
    Sigma_con[1, ] <- 0
    Sigma_con[, 1] <- 0
    Sigma_con[1, 1] <- 20
  }
  
  y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
  wX_dis <- cbind(w_dis, X_dis)
  
  if (prior_type == "normal") {
    if (all(diag(Sigma_con) == 20)) {
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
  } else if (prior_type == "G prior") {
    if (all(diag(Sigma_con) == 20)) {
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
      mu = b_con,
      Sigma = Sigma_con
    )
    
    discordant_model <- tryCatch(
      {
        # for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
        summary(rstan::sampling(sm_g, data = data_list, refresh = 0, chains = 1))$summary["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  } else if (prior_type == "PMP") {
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
  } else if (prior_type == "Hybrid") {
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
    
    # # 2. Build the data list to match the Stan model variables
    # data_list = list(
    #   N = nrow(wX_dis),
    #   K = ncol(wX_dis),
    #   X = wX_dis,
    #   y = y_dis_0_1,
    #   mu_T = 0,         # Center treatment prior at 0
    #   V_T = 20,        # High variance for treatment (non-informative)
    #   mu_X = mu_X_prior,
    #   Sigma_X = Sigma_X_prior
    # )
    
    discordant_model <- tryCatch(
      {
        summary(rstan::sampling(sm_hybrid, data = data_list, refresh = 0, chains = 1))$summary["beta_w", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("Hybrid PMP failed: %s", e$message))
        NULL
      }
    )
  }
  
  ret <- list(
    betaT = NA,
    ssq_beta_T = NA,
    reject = NA,
    pval = NA
  )
  
  if (!is.null(discordant_model)) {
    ret$betaT <- discordant_model[1]
    ret$ssq_beta_T <- discordant_model[2]
    ret$reject <- !(discordant_model[3] < 0 & discordant_model[4] > 0)
    ret$pval <- if (ret$reject) 0 else 1
  }
  
  return(ret)
}

Do_Inference <- function(y, X, w, strat, n) {
  res <- data.frame(
    n = numeric(),
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric()
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
  beta_hat_T <- NA
  ssq_beta_hat_T <- NA
  pval <- NA
  if (discordant_viabele) {
    y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
    model <- summary(glm(y_dis_0_1 ~ 0 + w_dis + X_dis, family = "binomial"))$coefficients[1, c(1, 2)]
    beta_hat_T <- model[1]
    ssq_beta_hat_T <- model[2]
    pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
  }
  res <- rbind(res, data.frame(
    n = n,
    inference = "clogit",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval
  ))
  
  ########################### LOGIT  ###########################
  beta_hat_T <- NA
  ssq_beta_hat_T <- NA
  pval <- NA
  if (TRUE) {
    model <- summary(glm(y ~ w + X, family = "binomial"))$coefficients[2, c(1, 2)]
    beta_hat_T <- model[1]
    ssq_beta_hat_T <- model[2]
    pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))
  }
  res <- rbind(res, data.frame(
    n = n,
    inference = "logit",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval
  ))
  
  ########################### BAYESIAN ###########################
  for (prior_type in c("normal", "G prior", "PMP", "Hybrid")) {
    for (concordant_fit in c("GLM", "GEE", "GLMM")) {
      beta_hat_T <- NA
      ssq_beta_hat_T <- NA
      pval <- NA
      pval_freq <- NA
      if (discordant_viabele & concordant_viabele) {
        model <- Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit)
        beta_hat_T <- model$betaT
        ssq_beta_hat_T <- model$ssq_beta_T
        pval <- model$pval
        # pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
      }
      res <- rbind(res, data.frame(
        n = n,
        inference = paste0("bayesian_", prior_type, "_", concordant_fit),
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval
      ))
    }
  }
  
  ########################### glmmTMB  ###########################
  beta_hat_T <- NA
  ssq_beta_hat_T <- NA
  pval <- NA
  tryCatch(
    {
      fit_tmb <- glmmTMB(
        y ~ X + w + (1 | strat),
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
    inference = "GLMM",
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
        y ~ X + w,
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
    inference = "GEE",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval
  ))
  
}

Run_sim <- function(n) {
  BIG_res <- data.frame()
  y <- array(NA, n)
  probs <- array(NA, n)
  
  pairs = sample(unique(diabetic$id), size = n/2)
  current_data = diabetic[diabetic$id %in% pairs, , drop = TRUE]
  strat = current_data$id
  w = current_data$trt
  y = current_data$status
  X = current_data[, c("eye", "risk", "time")]
  X$eye = factor(X$eye)
  X = model.matrix(~ 0 + ., X)
  
  one_res <- Do_Inference(y, X, w, strat, n)
  BIG_res <- rbind(BIG_res, one_res)

  return(BIG_res)
}


for (j in 1:120) {
  cat("################", j, "################\n")
  n = params[j,]$n
  print(Run_sim(n = n)); cat('\n')
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
      n <- row$n
      res <- tryCatch(
        {
          out <- Run_sim(n)
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
  write.csv(results, file = paste0("C:/temp/clogitR_kap_test_from_scratch/", Nsim, "_", e_nsim, "real.csv"), row.names = FALSE)
  rm(results)
  gc()
}

plan(sequential)


############### COMPILE RESULTS ###############

results <- read.csv("C:/temp/clogitR_kap_test_from_scratch/100_1.csv")

sum <- 1
for (i in 2:475) {
  file_path <- paste0("C:/temp/clogitR_kap_test_from_scratch/100_", i, ".csv")
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
    rej = pval < 0.05
  ) %>%
  group_by(p, beta_T, true_funtion, regress_on_X, n, inference, X_style) %>%
  arrange(desc(sq_err)) %>%
  slice(-(1:200)) %>%
  summarize(
    num_na = sum(is.na(pval)),
    num_real = sum(!is.na(pval)),
    mse = mean(sq_err, na.rm = TRUE, trim = 0.1),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    coverage = mean(covered, na.rm = TRUE),
    mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_beta_hat_T = mean(ssq_beta_hat_T, trim = 0.001, na.rm = TRUE),
    significant = binom.test(sum(rej, na.rm = TRUE), n = (n() - num_na),  p = 0.05)$p.value,
    .groups = "drop"
  )


write.csv(res_mod, file = "C:/temp/clogitR_kap_test_from_scratch/combined_2700.csv", row.names = FALSE)



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
