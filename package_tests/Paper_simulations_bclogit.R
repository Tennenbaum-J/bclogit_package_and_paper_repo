pacman::p_load(dplyr, tidyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, ggplot2, geepack, glmmTMB, rstan, binaryMM, rstanarm, ggdist)
rm(list = ls())
options(error = recover)

if (!require("bclogit", character.only = TRUE)) {
    remotes::install_local("bclogit", dependencies = FALSE, force = TRUE, upgrade = "never")
    library(bclogit)
}

############### VARIABLES ###############

num_cores <- availableCores() - 2
ps <- c(6) # number of covariates
Nsim <- 100 # number of simulations per block of simulations
external_nsim <- 100 # number of blocks of simulations
ns <- c(100, 250, 500) # n values
beta_Ts <- c(0, 0.5) # treatment effects
true_functions <- c("linear", "non-linear") # true underlying function either linear or nonlinear
regress_on_Xs <- c("one", "two") # number of covariates the model gets to see

params <- expand.grid(
    nsim = 1:Nsim,
    p = ps,
    beta_T = beta_Ts,
    n = ns
)
params <- params %>%
    arrange(nsim, p, beta_T, n)

############### CODE ###############

# Runs all considered types of inference on the simulated data
# called from `Run_sim`
Do_Inference <- function(y, X, w, strat, p, beta_T, n, true_function, regress_on_X) {
    # result data.frame
    res <- data.frame(
        n = numeric(),
        p = numeric(),
        beta_T = numeric(),
        true_function = character(),
        regress_on_X = character(),
        inference = character(),
        beta_hat_T = numeric(),
        pval = numeric(),
        g = numeric()
    )

    # Check viable conditions
    matched_data <- bclogit:::process_matched_pairs_cpp(
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

    discordant_viable <- if (length(y_dis) > ncol(X) + 7) TRUE else FALSE
    concordant_viable <- if (length(y_con) > ncol(X) + 7) TRUE else FALSE

    ########################### CLOGIT  ###########################
    beta_hat_T <- NA
    ssq_beta_hat_T <- NA
    pval <- NA
    if (discordant_viable) {
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
        true_function = true_function,
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
        true_function = true_function,
        regress_on_X = regress_on_X,
        inference = "GLM",
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval,
        g = NA
    ))

    ########################### BAYESIAN ###########################
    for (prior_type in c("Naive", "G prior", "PMP", "Hybrid")) {
        for (concordant_fit in c("GLM", "GEE", "GLMM")) {
            beta_hat_T <- NA
            ssq_beta_hat_T <- NA
            pval_hdr <- NA
            pval_cr <- NA
            g <- NA
            if (discordant_viable & concordant_viable) {
                fit <- tryCatch(
                    {
                        bclogit::bclogit(
                            y ~ X,
                            data = data.frame(y = y, X = X),
                            treatment = w,
                            strata = strat,
                            concordant_method = concordant_fit,
                            prior_type = prior_type,
                            chains = 1,
                            stan_refresh = 0
                        )
                    },
                    error = function(e) {
                        warning(sprintf("bclogit failed: %s", e$message))
                        NULL
                    }
                )

                if (!is.null(fit) && !is.null(fit$model)) {
                    beta_hat_T <- fit$coefficients[1]
                    ssq_beta_hat_T <- fit$var[1, 1]

                    # Compute P-value based on 95% HDR (Disjoint / Multimodal)
                    hpd_res <- confint(fit, parm = 1, type = "HPD_many", level = 0.95)
                    if (!is.null(hpd_res) && nrow(hpd_res) > 0) {
                        zero_in_hdi <- any(hpd_res[, 1] < 0 & hpd_res[, 2] > 0)
                        reject <- !zero_in_hdi
                        pval_hdr <- if (reject) 0 else 1
                    }

                    # Compute P-value based on 95% Credible Interval
                    mod_summ <- rstan::summary(fit$model)$summary
                    target_param <- ifelse(prior_type %in% c("Naive", "G prior"), "beta[1]", "beta_w")
                    if (target_param %in% rownames(mod_summ)) {
                        ci_lower <- mod_summ[target_param, "2.5%"]
                        ci_upper <- mod_summ[target_param, "97.5%"]
                        reject <- !(ci_lower < 0 & ci_upper > 0)
                        pval_cr <- if (reject) 0 else 1
                    }

                    if (prior_type %in% c("G prior", "Hybrid") && "g" %in% rownames(mod_summ)) {
                        g <- mod_summ["g", "mean"]
                    }
                }
            }
            res <- rbind(res, data.frame(
                n = n, p = p, beta_T = beta_T, true_function = true_function,
                regress_on_X = regress_on_X,
                inference = paste0("bayesian_", prior_type, "_", concordant_fit, "_HDR"),
                beta_hat_T = beta_hat_T, ssq_beta_hat_T = ssq_beta_hat_T, pval = pval_hdr, g = g
            ))
            res <- rbind(res, data.frame(
                n = n, p = p, beta_T = beta_T, true_function = true_function,
                regress_on_X = regress_on_X,
                inference = paste0("bayesian_", prior_type, "_", concordant_fit, "_CR"),
                beta_hat_T = beta_hat_T, ssq_beta_hat_T = ssq_beta_hat_T, pval = pval_cr, g = g
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
        true_function = true_function,
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
        true_function = true_function,
        regress_on_X = regress_on_X,
        inference = "GEE",
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval,
        g = NA
    ))
}

# simulate the data
Run_sim <- function(p, beta_T, n) {
    BIG_res <- data.frame()

    y <- array(NA, n)
    probs <- array(NA, n)

    X <- matrix(runif((n / 2) * p, min = -1, max = 1), ncol = p)
    X_plus_eps <- X + matrix(rnorm((n / 2) * p, 0, 0.05), ncol = p)
    combined <- rbind(X, X_plus_eps)
    ids <- order(c(1:(n / 2), 1:(n / 2)))
    X <- combined[ids, ]
    X[, 1] <- runif(n, min = -1, max = 1)
    rm(X_plus_eps, combined, ids)

    w <- c(rbind(replicate(n / 2, sample(c(0, 1)), simplify = TRUE)))
    strat <- rep(1:(n / 2), each = 2)

    for (true_function in true_functions) {
        if (true_function == "linear") {
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
            one_res <- Do_Inference(y, X_run, w, strat, p, beta_T, n, true_function, regress_on_X)
            BIG_res <- rbind(BIG_res, one_res)
        }
    }
    return(BIG_res)
}



for (j in 1:120) {
    cat("################", j, "################\n")
    p <- params[j, ]$p
    beta_T <- params[j, ]$beta_T
    n <- params[j, ]$n
    X_style <- params[j, ]$X_style
    print(Run_sim(p = p, beta_T = beta_T, n = n))
    cat("\n")
}

############### SIM ###############

handlers(global = TRUE)
handlers("txtprogressbar")
registerDoFuture()
plan(multisession, workers = num_cores)

for (e_nsim in 1:external_nsim) {
    with_progress({
        prog <- progressor(along = 1:nrow(params))

        results <- foreach(
            row = iter(params, by = "row"),
            .combine = rbind,
            .packages = c(
                "bclogit", "nbpMatching", "data.table",
                "dplyr", "MASS", "Rcpp", "rstanarm" # removed rstan from inside foreach to avoid issues
            )
        ) %dorng% {
            nsim <- row$nsim
            p <- row$p
            beta_T <- row$beta_T
            n <- row$n
            res <- tryCatch(
                {
                    out <- Run_sim(p, beta_T, n)
                    prog()
                    out
                },
                error = function(e) {
                    cat(glue::glue("Error in nsim={nsim}: {e$message}"), "\n")
                    prog()
                    NULL
                }
            )
        }
    })

    if (!dir.exists("C:/temp/simulations")) {
        dir.create("C:/temp/simulations", recursive = TRUE)
    }
    write.csv(results, file = paste0("C:/temp/simulations/", Nsim, "_", e_nsim, ".csv"), row.names = FALSE)
    rm(results)
    gc()
}

plan(sequential)

############### COMPILE RESULTS ###############

results <- read.csv("C:/temp/simulations/100_1.csv")
for (i in 2:101) {
    file_path <- paste0("C:/temp/simulations/100_", i, ".csv")
    if (file.exists(file_path)) {
        message("Reading file ", i)
        temp <- read.csv(file_path)
        results <- rbind(results, temp)
    }
}

res_mod <- results %>%
    mutate(
        lower_ci = beta_hat_T - (1.96 * ssq_beta_hat_T),
        upper_ci = beta_hat_T + (1.96 * ssq_beta_hat_T),
        covered = (lower_ci <= beta_T) & (upper_ci >= beta_T),
        sq_err = (beta_hat_T - beta_T)^2,
        rej = pval < 0.05
    ) %>%
    group_by(p, beta_T, true_function, regress_on_X, n, inference) %>%
    summarize(
        num_na = sum(is.na(pval)),
        num_real = sum(!is.na(pval)),
        mse = mean(sq_err, na.rm = TRUE, trim = 0.001),
        med_mse = median(sq_err, na.rm = TRUE),
        percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
        coverage = mean(covered, na.rm = TRUE),
        mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
        mean_sq_beta_hat_T = mean(ssq_beta_hat_T, trim = 0.001, na.rm = TRUE),
        significant = binom.test(sum(rej, na.rm = TRUE), n = (n() - num_na), p = 0.05)$p.value,
        coverage_significant = binom.test(x = sum(covered, na.rm = TRUE), n = sum(!is.na(covered)), p = 0.95)$p.value,
        mean_g = mean(g, na.rm = TRUE),
        .groups = "drop"
    )

write.csv(res_mod, file = "C:/temp/simulations/combined.csv", row.names = FALSE)

############### REAL DATA ###############

library(riskCommunicator)
data("framingham")
D <- data.table(framingham)
cols_to_check <- c("CIGPDAY", "BMI", "HEARTRTE", "TOTCHOL", "SYSBP", "DIABP", "CURSMOKE", "DIABETES", "BPMEDS")
for (col in cols_to_check) {
    D <- D[!is.na(get(col))]
}
D <- D[!is.na(PREVCHD)]

Dba <- D[PERIOD %in% c(1, 3)]
Dba[, num_periods_per_id := .N, by = RANDID]
Dba <- Dba[num_periods_per_id == 2]
Dba[, num_periods_per_id := NULL]
rm(framingham, D)

strat <- Dba$RANDID
w <- ifelse(Dba$PERIOD == 3, 1, 0)
y <- Dba$PREVCHD
X <- Dba[, ..cols_to_check]

res <- data.frame(
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric()
)

# CLOGIT
# Process matched pairs for clogit specifically as bclogit pkg wraps it
matched_data_real <- bclogit:::process_matched_pairs_cpp(
    strata = strat, y = y, X = data.matrix(X), treatment = w
)
X_dis_real <- matched_data_real$X_diffs_discordant
y_dis_real <- matched_data_real$y_diffs_discordant
w_dis_real <- matched_data_real$treatment_diffs_discordant

y_dis_0_1_real <- ifelse(y_dis_real == -1, 0, 1)
model <- summary(glm(y_dis_0_1_real ~ 0 + w_dis_real + X_dis_real, family = "binomial"))$coefficients[1, c(1, 2)]
beta_hat_T <- model[1]
ssq_beta_hat_T <- model[2]
pval <- 2 * pnorm(min(c(-1, 1) * (beta_hat_T / ssq_beta_hat_T)))

res <- rbind(res, data.frame(
    inference = "clogit",
    beta_hat_T = beta_hat_T,
    ssq_beta_hat_T = ssq_beta_hat_T,
    pval = pval
))

# BAYESIAN
for (prior_type in c("Naive", "G prior", "PMP", "Hybrid")) {
    for (concordant_fit in c("GLM", "GEE", "GLMM")) {
        beta_hat_T <- NA
        ssq_beta_hat_T <- NA
        pval_hdr <- NA
        pval_cr <- NA

        fit <- tryCatch(
            {
                bclogit::bclogit(
                    y = y,
                    X = data.frame(X),
                    treatment = w,
                    strata = strat,
                    concordant_method = concordant_fit,
                    prior_type = prior_type,
                    chains = 1,
                    stan_refresh = 0
                )
            },
            error = function(e) NULL
        )

        if (!is.null(fit) && !is.null(fit$model)) {
            beta_hat_T <- fit$coefficients[1]
            ssq_beta_hat_T <- fit$var[1, 1]

            hpd_res <- confint(fit, parm = 1, type = "HPD_many", level = 0.95)
            if (!is.null(hpd_res) && nrow(hpd_res) > 0) {
                zero_in_hdi <- any(hpd_res[, 1] < 0 & hpd_res[, 2] > 0)
                reject <- !zero_in_hdi
                pval_hdr <- if (reject) 0 else 1
            }

            # Compute P-value based on 95% Credible Interval
            mod_summ <- rstan::summary(fit$model)$summary
            target_param <- ifelse(prior_type %in% c("Naive", "G prior"), "beta[1]", "beta_w")
            if (target_param %in% rownames(mod_summ)) {
                ci_lower <- mod_summ[target_param, "2.5%"]
                ci_upper <- mod_summ[target_param, "97.5%"]
                reject <- !(ci_lower < 0 & ci_upper > 0)
                pval_cr <- if (reject) 0 else 1
            }
        }

        res <- rbind(res, data.frame(
            inference = paste0("bayesian_", prior_type, "_", concordant_fit, "_HDR"),
            beta_hat_T = beta_hat_T,
            ssq_beta_hat_T = ssq_beta_hat_T,
            pval = pval_hdr
        ))
        res <- rbind(res, data.frame(
            inference = paste0("bayesian_", prior_type, "_", concordant_fit, "_CR"),
            beta_hat_T = beta_hat_T,
            ssq_beta_hat_T = ssq_beta_hat_T,
            pval = pval_cr
        ))
        
        print(res)
    }
}
