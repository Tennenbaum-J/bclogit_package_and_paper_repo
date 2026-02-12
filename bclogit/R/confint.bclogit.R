#' Credible Intervals for bclogit Parameters
#'
#' Computes Bayesian credible intervals for the model parameters.
#'
#' @param object A `bclogit` object.
#' @param parm A specification of which parameters to be given credible intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required (default 0.95).
#' @param type Type of interval to compute: "quantile" (default), "HPD_one" (unimodal HPD interval via coda), "HPD_many" (multimodal HPD interval via ggdist).
#' @param ... Additional arguments.
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.bclogit <- function(object, parm, level = 0.95, type = "quantile", ...) {
    type <- match.arg(type)

    if (is.null(object$model)) {
        warning("No discordant model available for confidence intervals.")
        return(NULL)
    }

    sims <- rstan::extract(object$model)
    if (!is.null(sims$beta)) {
        beta_post <- sims$beta
    } else {
        beta_post <- cbind(as.vector(sims$beta_w), sims$beta_nuis)
    }

    if (ncol(beta_post) == length(object$coefficients)) {
        colnames(beta_post) <- names(object$coefficients)
    }

    if (missing(parm)) {
        parm <- names(object$coefficients)
    } else if (is.numeric(parm)) {
        parm <- names(object$coefficients)[parm]
    }

    present_parms <- intersect(parm, colnames(beta_post))
    if (length(present_parms) == 0) {
        stop("No parameters found matching the 'parm' argument.")
    }

    beta_subset <- beta_post[, present_parms, drop = FALSE]

    if (type == "quantile") {
        alpha <- (1 - level) / 2
        probs <- c(alpha, 1 - alpha)
        ci <- apply(beta_subset, 2, quantile, probs = probs)
        return(t(ci))
    } else if (type == "HPD_one") {
        # Standard Unimodal HPD Interval (coda)
        mcmc_obj <- coda::mcmc(beta_subset)
        ci <- coda::HPDinterval(mcmc_obj, prob = level)
        return(ci)
    } else {
        # "HPD_many" (ggdist) - supports multimodal / multiple intervals
        # Requires transforming to long format for ggdist::hdi
        # We return a tibble/data.frame here as the structure might be complex

        # Convert to tidy format
        beta_df <- as.data.frame(beta_subset)
        beta_long <- do.call(rbind, lapply(names(beta_df), function(p) {
            data.frame(Parameter = p, value = beta_df[[p]])
        }))

        # Calculate HDI
        # ggdist::hdi expects a vector of values. We can group by Parameter.
        # We can use dplyr if available or base split.

        # Using base R split + lapply to avoid heavy dplyr dependency inside the function body if possible,
        # but we imported tibble/ggdist so we likely can use them.
        # However, ggdist::hdi works on vectors.

        res_list <- lapply(split(beta_long, beta_long$Parameter), function(sub_df) {
            # Loop over levels to handle vector input safely
            do.call(rbind, lapply(level, function(l) {
                h <- ggdist::hdi(sub_df$value, .width = l)
                h_df <- as.data.frame(h)
                # ggdist returns matrix with unnamed columns or V1, V2
                # We standardize to lower/upper
                if (ncol(h_df) == 2) {
                    colnames(h_df) <- c("lower", "upper")
                } else {
                    # Fallback if structure is different
                    names(h_df) <- paste0("col", 1:ncol(h_df))
                }

                # Add metadata
                h_df$Parameter <- sub_df$Parameter[1]
                h_df$.width <- l
                h_df
            }))
        })

        res <- do.call(rbind, res_list)
        # Reorder columns to put Parameter first
        res <- res[, c("Parameter", setdiff(names(res), "Parameter"))]
        return(res)
    }
}
