#' Extract variance-covariance matrix from a bclogit model
#'
#' @param object A `bclogit` object.
#' @param ... Additional arguments.
#' @return A matrix of the estimated covariance of the coefficients.
#' @export
vcov.bclogit <- function(object, ...) {
    object$var
}
