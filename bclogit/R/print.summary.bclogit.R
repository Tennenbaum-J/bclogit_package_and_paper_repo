#' Print summary of a bclogit model
#'
#' @param x A `summary.bclogit` object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#' @export
print.summary.bclogit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Discordant Pairs Bayesian Model\n")
    cat("-------------------------------\n")
    cat("Discordant pairs: ", x$num_discordant, "\n")
    cat("Concordant pairs: ", x$num_concordant, " (used for prior)\n")

    if (!is.null(x$prior_info$mu)) {
        cat("\nPrior Estimates (from concordant data):\n")

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
    # Identify coefficient columns for special formatting (Estimate, Median, Std. Error, Quantiles)
    # Estimate=1, Median=2, Std. Error=3, Q_low=4, Q_high=5
    # Rhat=6, n_eff=7 (if present), Pr(!=0)=Last

    # We pass cs.ind to ensure n_eff (large integer) doesn't mess up significant digits of coefficients
    # Check if Rhat/n_eff exist to determine indices
    n_cols <- ncol(x$coefficients)
    has_rhat <- "Rhat" %in% colnames(x$coefficients)

    # Assuming standard structure: Est, Med, SE, QL, QH, [Rhat, n_eff], Pval
    # Coefficients are 1:5.

    printCoefmat(x$coefficients,
        digits = digits, signif.stars = TRUE,
        P.values = TRUE, has.Pvalue = TRUE, cs.ind = 1:5
    )

    cat("\nColumn 'Pr(!=0)' is the probability that the coefficient is not zero (based on HPD interval).\n")
    if ("Rhat" %in% colnames(x$coefficients)) {
        cat("Columns 'Rhat' and 'n_eff' are convergence diagnostics (R-hat should be < 1.05).\n")
    }

    invisible(x)
}
