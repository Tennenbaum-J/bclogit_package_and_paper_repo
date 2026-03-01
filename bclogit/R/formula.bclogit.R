#' Extract model formula
#'
#' @param x A `bclogit` object.
#' @param ... Additional arguments.
#' @return The \code{formula} used in the model, or \code{NULL} when the model was
#'   fitted via the default (matrix) interface.
#' @seealso \code{\link{bclogit}}
#' @examples
#' \donttest{
#' data("fhs")
#' fit <- bclogit(PREVHYP ~ TOTCHOL + BMI, data = fhs,
#'                treatment = PERIOD, strata = RANDID)
#' formula(fit)
#' }
#' @export
formula.bclogit <- function(x, ...) {
    if (!is.null(x$terms)) {
        formula(x$terms)
    } else {
        NULL
    }
}
