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

    # Probability of being positive/negative (for hypothesis testing)
    # Pr(beta > 0)
    prob_pos <- colMeans(beta_post > 0)

    # Extract diagnostics if available
    rhat <- NULL
    ess <- NULL

    if (!is.null(stan_summary)) {
        # Match parameter names from stan summary to our coefficients
        # This can be tricky because stan names are generic (beta[1], etc)
        # while our coefficients have meaningful names.
        # However, the order in beta_post should match the order we extracted.

        # We'll just take the corresponding rows from the summary if we can identify them.
        # For naive/G-prior: "beta[1]", "beta[2]"...
        # For PMP/Hybrid: "beta_w", "beta_nuis[1]", ...

        # Simplified approach: assume order in summary (excluding lp__) matches?
        # No, stan summary includes many things.

        # Let's try to map by index since we know the structure
        # But `rstan::extract` flattens things differently than `summary`.

        # Let's just calculate ESS and Rhat manually from chains if needed, or use `monitor`.
        # `rstan::monitor` takes an array of simulations.
        # But `object$model` is a stanfit object.

        # Let's try to get them from summary by matching names.
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
        `Est.Error` = se,
        `Q_low` = q_low,
        `Q_high` = q_high,
        `Pr(>0)` = prob_pos
    )

    if (!is.null(rhat)) {
        coef_mat <- cbind(coef_mat, Rhat = rhat, n_eff = ess)
    }

    # Set row names
    rownames(coef_mat) <- names(object$coefficients)

    # Rename Quantile columns nicely
    # Note: Column indices shifted by 1 due to Median insertion
    colnames(coef_mat)[4:5] <- paste0(c("L", "U"), format(conf.level * 100), "%")

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
