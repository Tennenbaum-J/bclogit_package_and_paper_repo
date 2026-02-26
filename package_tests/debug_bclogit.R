pacman::p_load(dplyr, tidyr, data.table, ggplot2, geepack, glmmTMB, rstan, binaryMM, rstanarm, ggdist)
rm(list = ls())

if (!require("bclogit", character.only = TRUE)) {
    remotes::install_local("bclogit", dependencies = FALSE, force = TRUE, upgrade = "never")
    library(bclogit)
}

source("C:/Users/Jacob/bclogit_package_and_paper_repo/package_tests/Paper_simulations_bclogit.R")

cat("\n--- Running a single sequential simulation ---\n")
p <- 6
beta_T <- 0
n <- 100

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

beta_X_value <- 1.25
beta_X <- rep(beta_X_value, p)
beta_0 <- -0.5
probs <- 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_X) + beta_T * w)))
y <- rbinom(n, 1, probs)

X_run <- X[, 1, drop = FALSE]

cat("Starting inference runs explicitly...\n")

matched_data <- bclogit:::process_matched_pairs_cpp(
    strata = strat,
    y = y,
    X = data.matrix(X_run),
    treatment = w
)

w_dis <- matched_data$treatment_diffs_discordant
X_dis <- matched_data$X_diffs_discordant
y_dis <- matched_data$y_diffs_discordant
y_con <- matched_data$y_concordant

y_dis_0_1 <- ifelse(y_dis == -1, 0, 1)

for (prior_type in c("Naive", "G prior", "PMP", "Hybrid")) {
    for (concordant_fit in c("GLM", "GEE", "GLMM")) {
        cat("\n--------------------------\n")
        cat("Running prior:", prior_type, "concordant_fit:", concordant_fit, "\n")
        fit <- tryCatch(
            {
                bclogit::bclogit(
                    y = y,
                    X = data.frame(X_run),
                    treatment = w,
                    strata = strat,
                    concordant_method = concordant_fit,
                    prior_type = prior_type,
                    chains = 1,
                    stan_refresh = 0
                )
            },
            error = function(e) {
                cat("bclogit failed with error:", e$message, "\n")
                NULL
            }
        )

        if (!is.null(fit) && !is.null(fit$model)) {
            cat("beta_hat_T:", fit$coefficients[1], "\n")
            cat("ssq_beta:", fit$var[1, 1], "\n")

            hpd_res <- confint(fit, parm = 1, type = "HPD_many", level = 0.95)
            cat("hpd is valid?\n")
            print(!is.null(hpd_res))
        }
    }
}
cat("\nDone with loop.\n")
