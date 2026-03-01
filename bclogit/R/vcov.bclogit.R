#' Extract variance-covariance matrix from a bclogit model
#'
#' @param object A `bclogit` object.
#' @param ... Additional arguments.
#' @return A square matrix of the posterior covariance of the coefficients,
#'   derived from the MCMC samples.
#' @seealso \code{\link{coef.bclogit}}, \code{\link{summary.bclogit}}
#' @examples
#' \donttest{
#' data("fhs")
#' fit <- bclogit(PREVHYP ~ TOTCHOL + BMI, data = fhs,
#'                treatment = PERIOD, strata = RANDID)
#' vcov(fit)
#' }
#' @export
vcov.bclogit <- function(object, ...) {
    object$var
}
