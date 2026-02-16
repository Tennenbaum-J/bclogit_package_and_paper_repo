#' Frequentist Conditional Logistic Regression
#'
#' Fits a conditional logistic regression model for matched pairs using
#' the discordant-pair GLM trick. This is a fast frequentist alternative
#' to \code{\link{bclogit}}.
#'
#' @param formula For the formula method, a symbolic description of the model.
#' @param y For the default method, a binary (0,1) response vector.
#' @param X A data.frame, data.table, or model.matrix containing the variables.
#' @param treatment Vector specifying the treatment variable.
#' @param strata Vector specifying the strata (matched pairs).
#' @param data A data frame containing the variables (for formula method).
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action A function which indicates what should happen when the data contain NAs. 
#' @param treatment_name Optional string name for the treatment variable.
#' @param call Optional call object to store in the result.
#' @param ... Additional arguments passed to methods.
#' @return An object of class \code{"clogit_bclogit"}.
#' @seealso \code{\link{bclogit}}, \code{\link{summary.clogit_bclogit}}
#' @examples
#' \dontrun{
#' n <- 200
#' dat <- data.frame(
#'   y = rbinom(n, 1, 0.5),
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   treatment = rep(c(0, 1), n / 2),
#'   strata = rep(1:(n / 2), each = 2)
#' )
#' fit <- clogit(y ~ x1 + x2, data = dat, treatment = treatment, strata = strata)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' }
#' @export
clogit <- function(formula, data, treatment = NULL, strata = NULL, 
                   subset = NULL, na.action = NULL,
                   do_inference_on_var = "all", ...) {
    UseMethod("clogit")
}
