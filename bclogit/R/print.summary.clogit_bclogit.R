#' Print summary of a clogit_bclogit model
#'
#' @param x A \code{summary.clogit_bclogit} object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#' @examples
#' \dontrun{
#' n <- 200
#' dat <- data.frame(
#'   y = rbinom(n, 1, 0.5), x1 = rnorm(n),
#'   treatment = rep(c(0, 1), n / 2),
#'   strata = rep(1:(n / 2), each = 2)
#' )
#' fit <- clogit(y ~ x1, data = dat, treatment = treatment, strata = strata)
#' print(summary(fit))
#' }
#' @export
print.summary.clogit_bclogit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    assertClass(x, "summary.clogit_bclogit")
    assertCount(digits, positive = TRUE)
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Conditional Logistic Regression (Frequentist)\n")
    cat("----------------------------------------------\n")
    cat("Total observations:  ", x$n, "\n")
    cat("Discordant pairs:    ", x$num_discordant, "\n")
    cat("Concordant pairs:    ", x$num_concordant, "\n")
    cat("\nCoefficients:\n")

    printCoefmat(x$coefficients,
        digits = digits, signif.stars = TRUE,
        P.values = TRUE, has.Pvalue = TRUE
    )

    invisible(x)
}
