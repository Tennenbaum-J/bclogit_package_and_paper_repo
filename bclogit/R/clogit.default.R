#' @return A list of class \code{"clogit_bclogit"} containing:
#'   \item{coefficients}{Estimated coefficients (posterior means).}
#'   \item{var}{Variance-covariance matrix of the coefficients (diagonal, built from standard errors).
#'     \code{NULL} when \code{do_inference_on_var} is not \code{"all"}.}
#'   \item{flr_model}{The fitted fast logistic regression model object returned by
#'     \code{fastLogisticRegressionWrap::fast_logistic_regression}.}
#'   \item{call}{The function call.}
#'   \item{terms}{The model terms.}
#'   \item{n}{Total number of observations.}
#'   \item{num_discordant}{Number of discordant pairs used for fitting.}
#'   \item{num_concordant}{Number of concordant pairs.}
#'   \item{X_model_matrix_col_names}{Column names of the covariate model matrix.}
#'   \item{treatment_name}{Name of the treatment variable.}
#'   \item{se}{Standard errors of the coefficients.}
#'   \item{z}{Z-statistics for each coefficient.}
#'   \item{pval}{Approximate p-values for each coefficient.}
#'   \item{do_inference_on_var}{The value of the \code{do_inference_on_var} argument.}
#' @describeIn clogit Default method for matrix/data input.
#' @examples
#' data("fhs")
#' fit <- clogit(PREVHYP ~ TOTCHOL + CIGPDAY + BMI + HEARTRTE, 
#'   data = fhs, treatment = PERIOD, strata = RANDID)
#' summary(fit)
#' @param subset An optional vector specifying a subset of observations.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param do_inference_on_var Which variable(s) to compute standard errors for.
#'   \code{"all"} (default) computes SEs for all coefficients. An integer \code{j}
#'   computes the SE only for the \code{j}th coefficient (1 = treatment, then covariates
#'   in order). \code{"none"} skips inference entirely.
#' @export
clogit.default <- function(formula = NULL,
                           data = NULL,
                           treatment = NULL,
                           strata = NULL,
                           subset = NULL,
                           na.action = NULL,
                           do_inference_on_var = "all",
                           ...,
                           y = NULL,
                           X = NULL,
                           treatment_name = NULL,
                           call = NULL) {

  # --------------------------------------------------------------------------
  # 1. Input Validation
  # --------------------------------------------------------------------------
  if (missing(treatment)) stop("The 'treatment' argument is required.")
  assertVector(subset, null.ok = TRUE)
  assertFunction(na.action, null.ok = TRUE)
  assertString(treatment_name, null.ok = TRUE)
  assert(
      checkString(do_inference_on_var),
      checkCount(do_inference_on_var, positive = TRUE)
  )

  if (test_data_frame(X, types = c("numeric", "integer", "factor", "logical"))) {
    data_mat <- model.matrix(~ 0 + ., data = X)
  } else if (test_matrix(X, mode = "numeric")) {
    if (sd(X[, 1]) == 0) {
      X[, 1] <- NULL
    }
    data_mat <- as.matrix(X)
  } else {
    assert(
      check_data_frame(X),
      check_matrix(X, mode = "numeric")
    )
  }

  n <- nrow(data_mat)

  assertNumeric(y, lower = 0, upper = 1, any.missing = FALSE, len = n)
  assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)
  assertNumeric(strata, any.missing = FALSE, len = n)

  if (!all(treatment %in% c(0, 1))) {
    stop("Treatment must be binary 0 or 1.")
  }
  if (is.null(treatment_name)) {
    treatment_name <- deparse(substitute(treatment))
  }

  model_terms <- tryCatch(terms(y ~ ., data = as.data.frame(data_mat)), error = function(e) NULL)

  # --------------------------------------------------------------------------
  # 2. Data Preparation via C++ matched pair processing
  # --------------------------------------------------------------------------
  X <- data_mat
  X_model_matrix_col_names <- colnames(X)
  if (is.null(X_model_matrix_col_names)) {
    X_model_matrix_col_names <- paste0("X", 1:ncol(X))
  }

  matched_data <- process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )

  X_diffs <- matched_data$X_diffs_discordant
  y_diffs <- matched_data$y_diffs_discordant
  treatment_diffs <- matched_data$treatment_diffs_discordant

  num_concordant <- length(matched_data$y_concordant) / 2
  num_discordant <- length(y_diffs)

  if (num_discordant < ncol(X_diffs) + 5) {
    stop("There are not enough discordant pairs. The model will not be fit.")
  }

  # --------------------------------------------------------------------------
  # 3. Fit logistic regression on discordant pair differences (no intercept)
  # --------------------------------------------------------------------------
  y_01 <- ifelse(y_diffs == -1, 0, 1)
  X_full <- cbind(treatment_diffs, X_diffs)
  col_names_final <- c(treatment_name, X_model_matrix_col_names)
  colnames(X_full) <- col_names_final

  flr_fit <- fastLogisticRegressionWrap::fast_logistic_regression(
    Xmm = X_full,
    ybin = y_01,
    do_inference_on_var = do_inference_on_var,
    ...
  )

  # --------------------------------------------------------------------------
  # 4. Result Construction
  # --------------------------------------------------------------------------
  coefficients <- flr_fit$coefficients
  names(coefficients) <- col_names_final

  # Build vcov from SEs (diagonal; off-diagonals only available with "all")
  se_vec <- flr_fit$se
  if (do_inference_on_var == "all") {
    # Full vcov: use se to build diagonal (off-diagonal not provided by fastLR)
    var_cov <- diag(se_vec^2, nrow = length(se_vec))
    rownames(var_cov) <- col_names_final
    colnames(var_cov) <- col_names_final
  } else {
    var_cov <- NULL
  }

  res <- list(
    coefficients = coefficients,
    var = var_cov,
    flr_model = flr_fit,
    call = if (!is.null(call)) call else match.call(),
    terms = model_terms,
    n = n,
    num_discordant = num_discordant,
    num_concordant = num_concordant,
    X_model_matrix_col_names = X_model_matrix_col_names,
    treatment_name = treatment_name,
    se = se_vec,
    z = flr_fit$z,
    pval = flr_fit$approx_pval,
    do_inference_on_var = do_inference_on_var
  )

  class(res) <- c("clogit_bclogit", "list")
  return(res)
}
