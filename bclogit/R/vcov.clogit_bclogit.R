#' Extract variance-covariance matrix from a clogit_bclogit model
#'
#' @param object A \code{clogit_bclogit} object.
#' @param ... Additional arguments.
#' @return A diagonal matrix of the estimated coefficient variances (squared standard
#'   errors), or \code{NULL} when \code{do_inference_on_var} is not \code{"all"}.
#' @seealso \code{\link{coef.clogit_bclogit}}, \code{\link{summary.clogit_bclogit}}
#' @examples
#' n <- 200
#' dat <- data.frame(
#'   y = rbinom(n, 1, 0.5),
#'   x1 = rnorm(n),
#'   treatment = rep(c(0, 1), n / 2),
#'   strata = rep(seq_len(n / 2), each = 2)
#' )
#' fit <- clogit(y ~ x1, data = dat, treatment = treatment, strata = strata)
#' vcov(fit)
#' @export
vcov.clogit_bclogit <- function(object, ...) {
    assertClass(object, "clogit_bclogit")
    object$var
}
