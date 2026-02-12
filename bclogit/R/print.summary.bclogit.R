#' Print summary of a bclogit model
#'
#' @param x A `summary.bclogit` object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#' @export
print.summary.bclogit <- function(x, digits = 4, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Discordant Pairs Bayesian Model\n")
    cat("-------------------------------\n")
    cat("Discordant pairs: ", x$num_discordant, "\n")
    cat("Concordant pairs: ", x$num_concordant, " (used for prior)\n")

    if (!is.null(x$prior_info$mu)) {
        cat("\nPrior Information (from concordant data):\n")

        prior_mu <- x$prior_info$mu
        if (length(prior_mu) == nrow(x$coefficients)) {
            names(prior_mu) <- rownames(x$coefficients)
            print(prior_mu, digits = digits)
        } else {
            cat("Mean vector length:", length(prior_mu), "\n")
            print(prior_mu, digits = digits)
        }
    }
    cat("\n")

    cat("Coefficients (Posterior Mean and Quantiles):\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, P.values = FALSE, has.Pvalue = FALSE)

    cat("\nColumn 'Pr(>0)' is the posterior probability that the coefficient is positive.\n")
    if ("Rhat" %in% colnames(x$coefficients)) {
        cat("Columns 'Rhat' and 'n_eff' are convergence diagnostics (R-hat should be < 1.05).\n")
    }

    invisible(x)
}
