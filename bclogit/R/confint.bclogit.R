#' Credible Intervals for bclogit Parameters
#'
#' Computes Bayesian credible intervals for the model parameters.
#'
#' @param object A `bclogit` object.
#' @param parm A specification of which parameters to be given credible intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required (default 0.95).
#' @param type Type of interval to compute: "HPD_one" (default unimodal HPD interval via coda), "CR" (equal-tailed credible region), "HPD_many" (multimodal HPD interval via ggdist).
#' @param ... Additional arguments.
#' @return A matrix with columns \code{lower} and \code{upper}.
#'   For \code{"HPD_many"}, a parameter may appear on multiple rows when the interval is disjoint.
#'   The matrix has a \code{Probability} attribute.
#' @seealso \code{\link{summary.bclogit}}, \code{\link{coef.bclogit}}
#' @export
confint.bclogit <- function(object, parm, level = 0.95, type = c("HPD_one", "CR", "HPD_many"), ...) {
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

    # Helper to normalize ggdist::hdi output column names
    normalize_hdi <- function(h) {
        h_df <- as.data.frame(h)
        if (ncol(h_df) >= 2) {
            colnames(h_df)[1:2] <- c("lower", "upper")
        }
        h_df
    }

    if (type == "CR") {
        alpha <- (1 - level) / 2
        probs <- c(alpha, 1 - alpha)
        res <- t(apply(beta_subset, 2, quantile, probs = probs))
        colnames(res) <- c("lower", "upper")
    } else if (type == "HPD_one") {
        # Standard Unimodal HPD Interval (coda)
        mcmc_obj <- coda::mcmc(beta_subset)
        res <- coda::HPDinterval(mcmc_obj, prob = level)
    } else {
        # "HPD_many" (ggdist) - supports multimodal / disjoint intervals
        if (!requireNamespace("ggdist", quietly = TRUE)) {
            stop("Package 'ggdist' is required for 'HPD_many' intervals. Please install it or use 'HPD_one'.")
        }
        res_list <- lapply(present_parms, function(p) {
            h_df <- normalize_hdi(ggdist::hdi(beta_subset[, p], .width = level))
            mat <- as.matrix(h_df[, c("lower", "upper"), drop = FALSE])
            rownames(mat) <- rep(p, nrow(mat))
            mat
        })
        res <- do.call(rbind, res_list)
    }
    attr(res, "Probability") <- level
    return(res)
}
