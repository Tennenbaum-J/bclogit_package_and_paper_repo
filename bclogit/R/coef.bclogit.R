#' Extract coefficients from a bclogit model
#'
#' @param object A `bclogit` object.
#' @param ... Additional arguments.
#' @return Named numeric vector of posterior mean coefficients.
#' @seealso \code{\link{vcov.bclogit}}, \code{\link{summary.bclogit}}
#' @examples
#' \donttest{
#' data("fhs")
#' fit <- bclogit(PREVHYP ~ TOTCHOL + BMI, data = fhs,
#'                treatment = PERIOD, strata = RANDID)
#' coef(fit)
#' }
#' @export
coef.bclogit <- function(object, ...) {
    object$coefficients
}
