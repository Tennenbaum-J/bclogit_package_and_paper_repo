#' Initialize a new bclogit model
#'
#' This function fits a Bayesian conditional logistic regression model, incorporating
#' information from concordant pairs to improve estimation.
#'
#' @param formula For the formula method, a symbolic description of the model to be fitted. 
#' @param y For the default method, a binary (0,1) vector containing the response of each subject.
#' @param data A data.frame, data.table, or model.matrix containing the variables (optional for formula method).
#' @param treatment Optional vector specifying the treatment variable (required for default method, or can be specified in formula method).
#' @param strata Vector specifying the strata (matched pairs).
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action A function which indicates what should happen when the data contain NAs. 
#' @param X A data.frame, data.table, or model.matrix containing the variables (optional for formula method, required for default method).
#' @param concordant_method The method to use for fitting the concordant pairs and reservoir. Options are "GLM", "GEE", and "GLMM".
#' @param prior_type The type of prior to use for the discordant pairs. Options are "Naive", "G prior", "PMP", and "Hybrid".
#' @param chains Number of chains for Stan sampling. Default is 4.
#' @param treatment_name Optional string name for the treatment variable.
#' @param return_raw_stan_output Logical; if \code{TRUE}, the raw Stan posterior samples
#'   (iterations x chains x parameters) are stored in the returned object. Default \code{FALSE}.
#' @param prior_variance_treatment Prior variance for the treatment coefficient in the
#'   covariance matrix \code{Sigma_con}. Default is 100.
#' @param stan_refresh How often Stan reports sampling progress (in iterations).
#'   Default is 0 (silent). Set to a positive integer (e.g., 1 or 100) to see progress.
#' @param call Optional call object to store in the result.
#' @param ... Additional arguments passed to \code{rstan::sampling} (e.g., \code{iter}, \code{warmup}, \code{thin}, \code{seed}, \code{control}).
#' @return A list of class `bclogit` containing:
#'   \item{coefficients}{Estimated coefficients (posterior means).}
#'   \item{var}{Variance-covariance matrix of coefficients.}
#'   \item{model}{The fitted Stan model object.}
#'   \item{posterior_samples}{Raw posterior samples as a 3D array (iterations x chains x parameters) from \code{rstan::extract(model, permuted = FALSE)}. Only populated when \code{return_raw_stan_output = TRUE}; \code{NULL} otherwise.}
#'   \item{concordant_model}{The fitted model object for the concordant pairs/reservoir (GLM/GEE/GLMM).}
#'   \item{matched_data}{The processed matched pairs data from the premodeling step.}
#'   \item{prior_info}{Information about the prior derived from concordant pairs.}
#'   \item{call}{The function call.}
#'   \item{terms}{The model terms.}
#'   \item{num_discordant}{Number of discordant pairs used.}
#'   \item{num_concordant}{Number of concordant pairs/reservoir entries used.}
#' @seealso \code{\link{summary.bclogit}}, \code{\link{confint.bclogit}},
#' \code{\link{vcov.bclogit}}, \code{\link{coef.bclogit}}
#' @examples
#' \donttest{
#' # Example usage
#' data("fhs")
#' fit <- bclogit(PREVHYP ~ TOTCHOL + CIGPDAY + BMI + HEARTRTE, 
#'   data = fhs, treatment = PERIOD, strata = RANDID)
#' summary(fit)
#' }
#' @export
bclogit <- function(formula, data, treatment = NULL, strata = NULL,
                    subset = NULL, na.action = NULL,
                    concordant_method = "GLM", prior_type = "Naive",
                    chains = 4, return_raw_stan_output = FALSE,
                    prior_variance_treatment = 100, stan_refresh = 0, ...) {
    UseMethod("bclogit")
}
