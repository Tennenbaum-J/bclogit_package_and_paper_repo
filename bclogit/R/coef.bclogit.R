#' Extract coefficients from a bclogit model
#'
#' @param object A `bclogit` object.
#' @param ... Additional arguments.
#' @return Numeric vector of coefficients.
#' @export
coef.bclogit <- function(object, ...) {
    object$coefficients
}
