pacman::p_load(dplyr, tidyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2, geepack, glmmTMB, rstan, binaryMM, rstanarm) # doParallel
rm(list = ls())

if (!require("bclogit", character.only = TRUE)) {
  remotes::install_local("bclogit", dependencies = FALSE, force = TRUE, upgrade = "never")
  library(bclogit)
}

############### VARIABLES ###############

num_cores <- availableCores() - 4
ps <- c(6)                                                  #number of covareates
Nsim <- 100                                                 #number of simulations per block of simulations
external_nsim <- 100                                        #number of blocks of simulations
ns <- c(100, 250, 500)                                      #n values
beta_Ts <- c(0, 0.5)                                        #treatment effects
true_funtions <- c("linear", "non-linear")                  #true underlying function either linear or nonlinear
regress_on_Xs <- c("one", "two")                            #number of covareats the model gets to see

# Priors  || Have to set proper directory
sm <- rstan::stan_model("mvn_logistic.stan")                #naive
sm_g <- rstan::stan_model("mvn_logistic_gprior.stan")       #g prior
sm_PMP <- rstan::stan_model("mvn_logistic_PMP.stan")        #PMP
sm_hybrid <- rstan::stan_model("mvn_logistic_Hybrid.stan")  #Hybrid


params <- expand.grid(
  nsim = 1:Nsim,
  p = ps,
  beta_T = beta_Ts,
  n = ns
)
params <- params %>%
  arrange(nsim, p, beta_T, n)

############### CODE ###############

# Runs bclogit without the package. This is how it works at the time of simulation for the paper
# called from `Do_Inference`
Bayesian_Clogit <- function(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit) {
  
  # NOTE ON THE CONCORDANT MODLES:
  # Since we include `w` in the concordant regression we have to get rid of a column from the estimate to match the dimension of the discordant model.
  # Additionaly we assume a flat prior for the treatment effect, hence setting it to mean 0 and high var.
  
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
    fit_con <- geeglm(
      y_con ~ w_con + X_con,
      id = strat_con,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      data = data.frame(y_con, w_con, X_con, strat_con)
    )

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

  # conditional logistic regression is the same as logistic regression as long as the response is set to be 0 or 1.
  y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
  wX_dis <- cbind(w_dis, X_dis)
  g = NA
  
  if (prior_type == "naive") {
    if (all(diag(Sigma_con) == 1e3 + eps)) {
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
        summary(rstan::sampling(sm, data = data_list, refresh = 0, chains = 1))$summary["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "G prior") {
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
      mu = b_con,
      Sigma = Sigma_con
    )

    discordant_model <- tryCatch(
      {
        mod = summary(rstan::sampling(sm_g, data = data_list, refresh = 0, chains = 1))$summary
        g = mod["g", "mean"]
        mod["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "PMP") {
    proj_matrix <- X_dis %*% solve(t(X_dis) %*% X_dis) %*% t(X_dis)
    w_dis_ortho <- w_dis - proj_matrix %*% w_dis

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
        summary(rstan::sampling(sm_PMP, data = data_list, refresh = 0, chains = 1))$summary["beta_w", c("mean", "sd", "2.5%", "97.5%")]
      },
      error = function(e) {
        warning(sprintf("stan_glm failed: %s", e$message))
        NULL
      }
    )
  }
  if (prior_type == "Hybrid") {
    proj_matrix <- X_dis %*% solve(t(X_dis) %*% X_dis) %*% t(X_dis)
    w_dis_ortho <- w_dis - proj_matrix %*% w_dis

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
        mod = summary(rstan::sampling(sm_hybrid, data = data_list, refresh = 0, chains = 1))$summary
        g = mod["g", "mean"]
        mod["beta_w", c("mean", "sd", "2.5%", "97.5%")]
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
    pval = NA,
    g = NA
  )

  if (!is.null(discordant_model)) {
    ret$betaT <- discordant_model[1]
    ret$ssq_beta_T <- discordant_model[2]
    ret$reject <- !(discordant_model[3] < 0 & discordant_model[4] > 0)
    ret$pval <- if (ret$reject) 0 else 1
  }
  
  if (!is.null(g)) { ret$g = g }

  return(ret)
}

# Runs all considered types of inference on the simulated data
# called from `Run_sim`
Do_Inference <- function(y, X, w, strat, p, beta_T, n, true_funtion, regress_on_X) {
  # result data.frame
  res <- data.frame(
    n = numeric(),
    p = numeric(),
    beta_T = numeric(),
    true_funtion = character(),
    regress_on_X = character(),
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric(),
    g = numeric()
  )
  
  # prepare the matched data
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

  # All types of inference in the paper
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
    p = p,
    beta_T = beta_T,
    true_funtion = true_funtion,
    regress_on_X = regress_on_X,
    inference = "clogit",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval,
    g = NA
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
    p = p,
    beta_T = beta_T,
    true_funtion = true_funtion,
    regress_on_X = regress_on_X,
    inference = "GLM",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval,
    g = NA
  ))
  
  ########################### BAYESIAN ###########################
  for (prior_type in c("naive", "G prior", "PMP", "Hybrid")) {
    for (concordant_fit in c("GLM", "GEE", "GLMM")) {
      beta_hat_T <- NA
      ssq_beta_hat_T <- NA
      pval <- NA
      beta_hat_T <- NA
      ssq_beta_hat_T <- NA
      pval <- NA
      g <- NA
      if (discordant_viabele & concordant_viabele) {
        model <- Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit)
        beta_hat_T <- model$betaT
        ssq_beta_hat_T <- model$ssq_beta_T
        pval <- model$pval
        g = model$g
      }
      res <- rbind(res, data.frame(
        n = n,
        p = p,
        beta_T = beta_T,
        true_funtion = true_funtion,
        regress_on_X = regress_on_X,
        inference = paste0("bayesian_", prior_type, "_", concordant_fit),
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval,
        g = g
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
    p = p,
    beta_T = beta_T,
    true_funtion = true_funtion,
    regress_on_X = regress_on_X,
    inference = "GLMM",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval,
    g = NA
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
    p = p,
    beta_T = beta_T,
    true_funtion = true_funtion,
    regress_on_X = regress_on_X,
    inference = "GEE",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval,
    g = NA
  ))


  rownames(res) <- NULL
  return(res)
}

# simulate the data
Run_sim <- function(p, beta_T, n) {
  # holds the results of all sub simulations
  BIG_res <- data.frame()
  
  y <- array(NA, n)
  probs <- array(NA, n)

  # See paper for full explanation of the data preparation
  X <- matrix(runif((n / 2) * p, min = -1, max = 1), ncol = p)
  X_plus_eps <- X + matrix(rnorm((n / 2) * p, 0, 0.05), ncol = p)
  combined <- rbind(X, X_plus_eps)
  ids <- order(c(1:(n / 2), 1:(n / 2)))
  X <- combined[ids, ]
  X[, 1] <- runif(n, min = -1, max = 1)
  rm(X_plus_eps, combined, ids)
    
  w <- c(rbind(replicate(n / 2, sample(c(0, 1)), simplify = TRUE)))
  strat <- rep(1:(n / 2), each = 2)


  for (true_funtion in true_funtions) {
    if (true_funtion == "linear") {
      beta_X_value <- 1.25
      beta_X <- rep(beta_X_value, p)
      beta_0 <- -0.5
      probs <- 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_X) + beta_T * w)))
    } else {
      f_x <- sin(pi * X[, 1] * X[, 2]) + X[, 3]^3 + X[, 4]^2 + X[, 5]^2
      probs <- 1 / (1 + exp(-(f_x + beta_T * w)))
    }
    y <- rbinom(n, 1, probs)

    for (regress_on_X in regress_on_Xs) {
      if (regress_on_X == "one") {
        X_run <- X[, 1, drop = FALSE]
      } else if (regress_on_X == "two") {
        X_run <- X[, c(1, 2), drop = FALSE]
      } else {
        X_run <- X
      }
      one_res <- Do_Inference(y, X_run, w, strat, p, beta_T, n, true_funtion, regress_on_X)
      BIG_res <- rbind(BIG_res, one_res)
    }
  }
  return(BIG_res)
}

############### SIM ###############

# set up the simulation
handlers(global = TRUE)
handlers("txtprogressbar")
registerDoFuture()
plan(multisession, workers = num_cores)


# This simulation uses iterative results dumping to minimize RAM use, and protect results
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
      res <- tryCatch(
        {
          out <- Run_sim(p, beta_T, n)
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
  write.csv(results, file = paste0("C:/temp/bclogit_simulation_results/", Nsim, "_", e_nsim, ".csv"), row.names = FALSE)
  rm(results)
  gc()
}

plan(sequential)


############### COMPILE RESULTS ###############

# retrieve the results
results <- read.csv("C:/temp/bclogit_simulation_results/100_1.csv")
for (i in 2:external_nsim) {
  file_path <- paste0("C:/temp/bclogit_simulation_results/100_", i, "_g_prior.csv")
  if (file.exists(file_path)) {
    message("Reading file ", i)
    temp <- read.csv(file_path)
    results <- rbind(results, temp)
  } else {
    message("Skipping missing file ", i)
  }
}

# compile results
res_mod <- results %>%
  mutate(
    lower_ci = beta_hat_T - (1.96 * ssq_beta_hat_T),
    upper_ci = beta_hat_T + (1.96 * ssq_beta_hat_T),
    covered = (lower_ci <= beta_T) & (upper_ci >= beta_T),
    sq_err = (beta_hat_T - beta_T)^2,
    rej = pval < 0.05
  ) %>%
  group_by(p, beta_T, true_funtion, regress_on_X, n, inference) %>%
  summarize(
    num_na = sum(is.na(pval)),
    num_real = sum(!is.na(pval)),
    mse = mean(sq_err, na.rm = TRUE, trim = 0.001),
    med_mse = median(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    coverage = mean(covered, na.rm = TRUE),
    mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_beta_hat_T = mean(ssq_beta_hat_T, trim = 0.001, na.rm = TRUE),
    significant = binom.test(sum(rej, na.rm = TRUE), n = (n() - num_na),  p = 0.05)$p.value,
    coverage_significant = binom.test(x = sum(covered, na.rm = TRUE), n = sum(!is.na(covered)), p = 0.95)$p.value,
    mean_g = mean(g, na.rm = TRUE),
    .groups = "drop"
  )


write.csv(res_mod, file = "C:/temp/bclogit_simulation_results/combined.csv", row.names = FALSE)



############### REAL DATA ###############

library(riskCommunicator)
data("framingham")
D=data.table(framingham)
D=D[!is.na(CIGPDAY)]
D=D[!is.na(BMI)]
D=D[!is.na(HEARTRTE)]
D=D[!is.na(TOTCHOL)]
D=D[!is.na(SYSBP)]
D=D[!is.na(DIABP)]
D=D[!is.na(CURSMOKE)]
D=D[!is.na(DIABETES)]
D=D[!is.na(BPMEDS)]

Dba = D[PERIOD %in% c(1,3)]
Dba[, num_periods_per_id := .N, by = RANDID]
Dba = Dba[num_periods_per_id == 2]
Dba[, num_periods_per_id := NULL]
rm(framingham, D)

strat = Dba$RANDID
w = ifelse(Dba$PERIOD == 3, 1, 0)
y = Dba$PREVCHD 
X = Dba[, c("TOTCHOL", "SYSBP", "DIABP", "CURSMOKE", "CIGPDAY", "BMI", "DIABETES", "BPMEDS", "HEARTRTE")]


res <- data.frame(
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
strat_con = strat[setdiff(1:length(strat), dis_idx)]

# CLOGIT  #
beta_hat_T <- NA
ssq_beta_hat_T <- NA
pval <- NA

y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)
model <- summary(glm(y_dis_0_1 ~ 0 + w_dis + X_dis, family = "binomial"))$coefficients[1, c(1, 2)]
beta_hat_T <- model[1]
ssq_beta_hat_T <- model[2]
pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))

res <- rbind(res, data.frame(
  inference = "clogit",
  beta_hat_T = beta_hat_T,
  ssq_beta_hat_T = ssq_beta_hat_T,
  pval = pval
))

# BAYESIAN  #
# This references the Bayesian_Clogit function from the CODE section
for (prior_type in c("naive", "G prior", "PMP", "Hybrid")) {
  for (concordant_fit in c("GLM", "GEE", "GLMM")) {
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    model <- Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit)
    beta_hat_T <- model$betaT
    ssq_beta_hat_T <- model$ssq_beta_T
    pval <- model$pval
    beta_hat_T ; ssq_beta_hat_T ; pval
    g = model$g
    res <- rbind(res, data.frame(
      inference = paste0("bayesian_", prior_type, "_", concordant_fit),
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
  }
}











