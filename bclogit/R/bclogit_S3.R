#' Initialize a new bclogit model
#'
#' This function fits a Bayesian conditional logistic regression model, incorporating
#' information from concordant pairs to improve estimation.
#'
#' @param formula A symbolic description of the model to be fitted (for formula method).
#' @param data A data.frame, data.table, or model.matrix containing the variables (optional for formula method).
#' @param treatment Optional vector specifying the treatment variable (required for default method, or can be specified in formula method).
#' @param strata Vector specifying the strata (matched pairs).
#' @param y A binary (0,1) vector containing the response of each subject (for default method).
#' @param X A data.frame, data.table, or model.matrix containing the variables (optional for formula method, required for default method).
#' @param concordant_method The method to use for fitting the concordant pairs and reservoir. Options are "GLM", "GEE", and "GLMM".
#' @param prior_type The type of prior to use for the discordant pairs. Options are "Naive", "G prior", "PMP", and "hybrid".
#' @param chains Number of chains for Stan sampling. Default is 4.
#' @param treatment_name Optional string name for the treatment variable.
#' @param ... Additional arguments passed to the default method.
#' @return A list of class `bclogit` containing:
#'   \item{coefficients}{Estimated coefficients (posterior means).}
#'   \item{var}{Variance-covariance matrix of coefficients.}
#'   \item{model}{The fitted Stan model object.}
#'   \item{concordant_model}{The fitted model object for the concordant pairs/reservoir (GLM/GEE/GLMM).}
#'   \item{prior_info}{Information about the prior derived from concordant pairs.}
#'   \item{call}{The function call.}
#'   \item{terms}{The model terms.}
#'   \item{num_discordant}{Number of discordant pairs used.}
#'   \item{num_concordant}{Number of concordant pairs/reservoir entries used.}
#' @seealso \code{\link{summary.bclogit}}, \code{\link{confint.bclogit}},
#' \code{\link{vcov.bclogit}}, \code{\link{coef.bclogit}}
#' @examples
#' \dontrun{
#' # Example usage
#' fit <- bclogit(y ~ x1 + x2, data = mydata, treatment = trt, strata = id)
#' summary(fit)
#' }
#' @export
bclogit <- function(x, ...) {
    UseMethod("bclogit")
}
