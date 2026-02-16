#' Extract coefficients from a clogit_bclogit model
#'
#' @param object A \code{clogit_bclogit} object.
#' @param ... Additional arguments.
#' @return Numeric vector of coefficients.
#' @examples
#' \dontrun{
#' n <- 200
#' dat <- data.frame(
#'   y = rbinom(n, 1, 0.5), x1 = rnorm(n),
#'   treatment = rep(c(0, 1), n / 2),
#'   strata = rep(1:(n / 2), each = 2)
#' )
#' fit <- clogit(y ~ x1, data = dat, treatment = treatment, strata = strata)
#' coef(fit)
#' }
#' @export
coef.clogit_bclogit <- function(object, ...) {
    assertClass(object, "clogit_bclogit")
    object$coefficients
}
