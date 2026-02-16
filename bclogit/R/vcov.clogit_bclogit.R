#' Extract variance-covariance matrix from a clogit_bclogit model
#'
#' @param object A \code{clogit_bclogit} object.
#' @param ... Additional arguments.
#' @return A matrix of the estimated covariance of the coefficients.
#' @export
vcov.clogit_bclogit <- function(object, ...) {
    assertClass(object, "clogit_bclogit")
    object$var
}
