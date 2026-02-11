#' Credible Intervals for bclogit Parameters
#'
#' Computes Bayesian credible intervals for the model parameters.
#'
#' @param object A `bclogit` object.
#' @param parm A specification of which parameters to be given credible intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required (default 0.95).
#' @param type Type of interval to compute: "quantile" (default) or "HPD" (Highest Posterior Density).
#' @param ... Additional arguments.
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.bclogit <- function(object, parm, level = 0.95, type = c("quantile", "HPD"), ...) {
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
    } else {
        # HPD Interval
        # Using coda::HPDinterval if valid mcmc object, but here we have a matrix
        # We can implement a simple HPD for unimodal posterior samples or default to coda if available.
        # To avoid extra dependencies if not needed, we can use the `coda` package if the user has it,
        # or implement a basic version.
        # Given this is a package, it's safer to rely on `coda` or `HDInterval` if we add them to Imports,
        # but for now I will implement a basic HPD calculation to avoid dependency hell unless requested.

        # Actually, `rstan` usually suggests `coda` or `loo`.
        # Let's check `coda` availability or use a simple method.

        compute_hpd <- function(x, prob = 0.95) {
            x <- sort(x)
            n <- length(x)
            gap <- max(1, min(n - 1, round(n * prob)))
            init <- 1:(n - gap)
            inds <- which.min(x[init + gap] - x[init])
            c(lower = x[init[inds]], upper = x[init[inds] + gap])
        }

        ci <- apply(beta_subset, 2, compute_hpd, prob = level)
        rownames(ci) <- c(paste0("lower (", level * 100, "%)"), paste0("upper (", level * 100, "%)"))
        return(t(ci))
    }
}
