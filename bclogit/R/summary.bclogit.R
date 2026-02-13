#' Summarize a bclogit model
#'
#' @param object A `bclogit` object.
#' @param conf.level Confidence level for credible intervals (default 0.95).
#' @param ... Additional arguments (not used).
#' @return A `summary.bclogit` object containing coefficients, standard errors, and credible intervals.
#' @export
summary.bclogit <- function(object, conf.level = 0.95, ...) {
    if (is.null(object$model)) {
        # Fallback if no model was fit / converged
        cat("No discordant model available.\n")
        return(invisible(object))
    }

    # Extract posterior samples again or use stored summary if optimization needed
    # We'll extract from standard 'stanfit' object
    if (!is.null(object$model)) {
        sims <- rstan::extract(object$model)
        stan_summary <- rstan::summary(object$model)$summary
    } else {
        sims <- list()
        stan_summary <- NULL
    }

    # Reconstruct beta_post matrix
    if (!is.null(sims$beta)) {
        beta_post <- sims$beta
        # Columns in stan_summary for beta
        # Usually they are named beta[1], beta[2], etc.
        # We need to map them correctly.
        # rstan summary has row names like "beta[1]", "beta[2]"...
        # and "beta_w[1]" etc for PMP.
    } else {
        beta_post <- cbind(as.vector(sims$beta_w), sims$beta_nuis)
    }

    # Calculate summary stats
    est <- colMeans(beta_post)
    median_val <- apply(beta_post, 2, median)
    se <- apply(beta_post, 2, sd)

    alpha <- (1 - conf.level) / 2
    q_low <- apply(beta_post, 2, quantile, probs = alpha)
    q_high <- apply(beta_post, 2, quantile, probs = 1 - alpha)

    # Calculate "p-value" based on HPD interval
    # We find the smallest alpha such that the (1-alpha) HDI includes 0 on the boundary.
    # Algorithm: Binary search for alpha using ggdist::hdi.

    calc_hpd_pval <- function(x) {
        if (all(x > 0) || all(x < 0)) {
            return(0)
        } # Shortcut if all samples are same sign (0 not in range)

        low <- 0
        high <- 1
        for (i in 1:20) { # Precision approx 1e-6
            mid <- (low + high) / 2
            # ggdist::hdi returns a data frame-like structure or vector
            # .width = 1 - mid.
            # Handle cases where width is 0 or 1 safely?
            if (mid > 0.999) mid <- 0.999
            if (mid < 0.001) mid <- 0.001

            # ggdist::hdi on a vector returns a tidy data frame
            # We need to extract lower/upper.
            h <- tryCatch(ggdist::hdi(x, .width = 1 - mid), error = function(e) NULL)

            if (is.null(h)) {
                break
            }

            # Check if 0 is in interval(s).
            # ggdist::hdi returns a matrix with columns [lower, upper]
            # It may return multiple rows if multimodal (but usually 1 row).
            # We check if 0 is contained in ANY of the intervals.

            in_interval <- FALSE
            if (length(h) > 0) {
                # Ensure it's treated as matrix even if 1 row / vector returned
                if (!is.matrix(h)) {
                    h <- matrix(h, ncol = 2)
                }

                # Check if 0 is within [lower, upper] for any row
                # Column 1 is lower, Column 2 is upper
                in_interval <- any(h[, 1] <= 0 & h[, 2] >= 0)
            }

            if (in_interval) {
                # 0 is inside (1-mid) interval.
                # We can try a narrower interval (higher alpha/mid).
                low <- mid
            } else {
                # 0 is outside. Need wider interval (lower alpha/mid).
                high <- mid
            }
        }
        return(low)
    }

    hpd_pvals <- apply(beta_post, 2, calc_hpd_pval)

    # Extract diagnostics if available
    rhat <- NULL
    ess <- NULL

    if (!is.null(stan_summary)) {
        # Match parameter names from stan summary to our coefficients
        all_row_names <- rownames(stan_summary)

        # Construct expected stan names
        if (!is.null(sims$beta)) {
            # Naive / G-prior
            param_names <- paste0("beta[", seq_len(ncol(beta_post)), "]")
        } else {
            # PMP / Hybrid
            # beta_w is scalar? check code.
            # "beta_w_post <- as.vector(sims$beta_w)" implies it might be scalar or vector.
            # In code: "parameters { vector[P] beta; }" for naive.
            # For PMP: "parameters { real beta_w; vector[P-1] beta_nuis; }"

            if (ncol(beta_post) > 1) {
                param_names <- c("beta_w", paste0("beta_nuis[", seq_len(ncol(beta_post) - 1), "]"))
            } else {
                param_names <- "beta_w"
            }
        }

        # Check if these exist in summary
        if (all(param_names %in% all_row_names)) {
            rhat <- stan_summary[param_names, "Rhat"]
            ess <- stan_summary[param_names, "n_eff"]
        }
    }

    coef_mat <- cbind(
        Estimate = est,
        Median = median_val,
        `Std. Error` = se,
        `Q_low` = q_low,
        `Q_high` = q_high
    )

    if (!is.null(rhat)) {
        coef_mat <- cbind(coef_mat, Rhat = rhat, n_eff = ess)
    }

    # Add P-value last so printCoefmat recognizes it for stars
    coef_mat <- cbind(coef_mat, `Pr(!=0)` = hpd_pvals)

    # Set row names
    rownames(coef_mat) <- names(object$coefficients)

    # Rename Quantile columns nicely to reflect actual quantiles
    # Columns 4 and 5 correspond to Q_low and Q_high
    pct_low <- paste0(format(alpha * 100, digits = 3, trim = TRUE), "%")
    pct_high <- paste0(format((1 - alpha) * 100, digits = 3, trim = TRUE), "%")
    colnames(coef_mat)[4:5] <- c(pct_low, pct_high)

    res <- list(
        call = object$call,
        coefficients = coef_mat,
        num_discordant = object$num_discordant,
        num_concordant = object$num_concordant,
        conf.level = conf.level,
        prior_info = object$prior_info,
        treatment_name = object$treatment_name
    )

    class(res) <- "summary.bclogit"
    res
}
