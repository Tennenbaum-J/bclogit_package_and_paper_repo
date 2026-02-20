#' Summarize a bclogit model
#'
#' @param object A `bclogit` object.
#' @param level Confidence level for credible intervals (default 0.95).
#' @param inference_method Method used for both the displayed confidence set bounds and the
#'   p-value (computed via bisection over alpha). Options are:
#'   \code{"HPD_one"} (default) uses unimodal HPD intervals (C++) with 20 bisection iterations,
#'   \code{"HPD_many"} uses \code{ggdist::hdi} which supports disjoint (multimodal) HPD regions,
#'     with 50 bisection iterations (requires the \pkg{ggdist} package).
#'     Confidence set bounds are shown when the HPD is a single interval; if disjoint,
#'     they are set to \code{NA} (use \code{\link{confint.bclogit}} with \code{type = "HPD_many"}
#'     to see all intervals),
#'   \code{"CR"} uses equal-tailed credible intervals (quantile-based, C++) with 20 bisection iterations.
#' @param ... Additional arguments (not used).
#' @return A list of class \code{"summary.bclogit"} containing:
#'   \item{call}{The original function call.}
#'   \item{coefficients}{A matrix with one row per parameter and columns for the posterior mean
#'     estimate, posterior median estimate, standard error, lower and upper credible interval bounds,
#'     optionally \code{Rhat} and \code{n_eff} convergence diagnostics (when available from Stan),
#'     and \code{Pr(tx!=0)} (the Bayesian p-value).}
#'   \item{num_discordant}{Number of discordant pairs used for fitting.}
#'   \item{num_concordant}{Number of concordant pairs used for the prior.}
#'   \item{level}{The credible interval level used.}
#'   \item{inference_method}{The inference method used for interval and p-value computation.}
#'   \item{prior_info}{A list with elements \code{mu} and \code{Sigma} describing the prior
#'     derived from the concordant pairs model.}
#'   \item{treatment_name}{Name of the treatment variable.}
#' @export
summary.bclogit <- function(object, level = 0.95, inference_method = "HPD_one", ...) {
    assertClass(object, "bclogit")
    assertNumber(level, lower = 0, upper = 1)
    assertChoice(inference_method, c("HPD_one", "HPD_many", "CR"))

    if (is.null(object$model)) {
        cat("No discordant model available.\n")
        return(invisible(object))
    }

    sims <- rstan::extract(object$model)
    stan_summary <- rstan::summary(object$model)$summary

    # Reconstruct beta_post matrix
    if (!is.null(sims$beta)) {
        beta_post <- sims$beta
    } else {
        beta_post <- cbind(as.vector(sims$beta_w), sims$beta_nuis)
    }

    # Calculate summary stats
    est <- colMeans(beta_post)
    median_val <- apply(beta_post, 2, median)
    se <- apply(beta_post, 2, sd)

    # Confidence set bounds
    alpha <- (1 - level) / 2
    ci_low_label <- paste0(format(alpha * 100, digits = 3, trim = TRUE), "%")
    ci_high_label <- paste0(format((1 - alpha) * 100, digits = 3, trim = TRUE), "%")

    if (inference_method == "CR") {
        ci_low <- apply(beta_post, 2, quantile, probs = alpha)
        ci_high <- apply(beta_post, 2, quantile, probs = 1 - alpha)
    } else if (inference_method == "HPD_one") {
        mcmc_obj <- coda::mcmc(beta_post)
        hpd_ci <- coda::HPDinterval(mcmc_obj, prob = level)
        ci_low <- hpd_ci[, "lower"]
        ci_high <- hpd_ci[, "upper"]
    } else {
        # HPD_many via ggdist - show bounds when unimodal, NA when disjoint
        if (!requireNamespace("ggdist", quietly = TRUE)) {
            stop("Package 'ggdist' is required for 'HPD_many'. Please install it or use 'HPD_one'.")
        }
        # Helper to normalize ggdist::hdi output to a data frame with lower/upper columns
        normalize_hdi <- function(h) {
            h_df <- as.data.frame(h)
            if (ncol(h_df) >= 2) {
                colnames(h_df)[1:2] <- c("lower", "upper")
            }
            h_df
        }
        n_coefs <- ncol(beta_post)
        ci_low <- numeric(n_coefs)
        ci_high <- numeric(n_coefs)
        for (j in seq_len(n_coefs)) {
            h_df <- normalize_hdi(ggdist::hdi(beta_post[, j], .width = level))
            if (nrow(h_df) == 1) {
                ci_low[j] <- h_df$lower
                ci_high[j] <- h_df$upper
            } else {
                # Disjoint intervals - cannot represent as single bounds
                ci_low[j] <- NA
                ci_high[j] <- NA
            }
        }
    }

    # P-value computation via bisection
    if (inference_method == "HPD_many") {
        # R-based bisection using ggdist::hdi which supports disjoint intervals
        n_iter <- 50L
        n_coefs <- ncol(beta_post)
        pvals <- numeric(n_coefs)
        for (j in seq_len(n_coefs)) {
            samps <- beta_post[, j]
            # If all samples are on one side of zero, p-value is 0
            if (min(samps) > 0 || max(samps) < 0) {
                pvals[j] <- 0
                next
            }
            low <- 0
            high <- 1
            for (i in seq_len(n_iter)) {
                mid <- (low + high) / 2
                mid <- max(0.001, min(0.999, mid))
                level_bis <- 1 - mid
                h_df <- normalize_hdi(ggdist::hdi(samps, .width = level_bis))
                in_interval <- any(h_df$lower <= 0 & h_df$upper >= 0)
                if (in_interval) {
                    low <- mid
                } else {
                    high <- mid
                }
            }
            pvals[j] <- low
        }
    } else {
        n_iter <- 20L
        interval_type <- if (inference_method == "CR") 1L else 0L  # 0 = HPD, 1 = CR
        pvals <- calc_bisection_pvals_cpp(beta_post, n_iter, interval_type)
    }

    # Extract diagnostics if available
    rhat <- NULL
    ess <- NULL

    if (!is.null(stan_summary)) {
        all_row_names <- rownames(stan_summary)

        if (!is.null(sims$beta)) {
            param_names <- paste0("beta[", seq_len(ncol(beta_post)), "]")
        } else {
            if (ncol(beta_post) > 1) {
                param_names <- c("beta_w", paste0("beta_nuis[", seq_len(ncol(beta_post) - 1), "]"))
            } else {
                param_names <- "beta_w"
            }
        }

        if (all(param_names %in% all_row_names)) {
            rhat <- stan_summary[param_names, "Rhat"]
            ess <- stan_summary[param_names, "n_eff"]
        }
    }

    coef_mat <- cbind(
        `mean estimate` = est,
        `median estimate` = median_val,
        `Std. Error` = se,
        ci_low,
        ci_high
    )
    colnames(coef_mat)[4:5] <- c(ci_low_label, ci_high_label)

    if (!is.null(rhat)) {
        coef_mat <- cbind(coef_mat, Rhat = rhat, n_eff = ess)
    }

    coef_mat <- cbind(coef_mat, `Pr(tx!=0)` = pvals)

    rownames(coef_mat) <- names(object$coefficients)

    res <- list(
        call = object$call,
        coefficients = coef_mat,
        num_discordant = object$num_discordant,
        num_concordant = object$num_concordant,
        level = level,
        inference_method = inference_method,
        prior_info = object$prior_info,
        treatment_name = object$treatment_name
    )

    class(res) <- "summary.bclogit"
    res
}
