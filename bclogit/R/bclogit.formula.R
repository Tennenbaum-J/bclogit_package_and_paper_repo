#' @return An object of class \code{"bclogit"}.
#' @export
#' @describeIn bclogit Formula method
bclogit.formula <- function(formula, data, treatment = NULL, strata = NULL,
                            subset = NULL, na.action = NULL,
                            concordant_method = "GLM", prior_type = "Naive",
                            chains = 4, return_raw_stan_output = FALSE,
                            prior_variance_treatment = 100, stan_refresh = 0, ...) {
    assertClass(formula, "formula")
    # subset and na.action are handled by model.frame; do not touch subset here to avoid premature evaluation
    assertFunction(na.action, null.ok = TRUE)
    assertChoice(concordant_method, c("GLM", "GEE", "GLMM"))
    assertChoice(prior_type, c("Naive", "G prior", "PMP", "Hybrid"))

    cl <- match.call(expand.dots = TRUE)
    m <- match(c("formula", "data", "subset", "na.action", "treatment", "strata"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(stats::model.frame)
    mf <- eval(cl, parent.frame())
    
    response <- model.response(mf, "numeric")
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)

    if (colnames(X)[1] == "(Intercept)") {
        X <- X[, -1, drop = FALSE]
    }

    trt_vec <- mf[["(treatment)"]]
    strata_vec <- mf[["(strata)"]]

    if (is.null(trt_vec)) stop("The 'treatment' argument is required.")
    if (is.null(strata_vec)) stop("The 'strata' argument is required.")

    # Extract treatment name from call
    full_cl <- match.call()
    if ("treatment" %in% names(full_cl)) {
        t_name <- deparse(full_cl$treatment)
    } else {
        t_name <- "treatment"
    }

    bclogit.default(
        y = response,
        X = X,
        treatment = trt_vec,
        strata = strata_vec,
        treatment_name = t_name,
        call = full_cl,
        concordant_method = concordant_method,
        prior_type = prior_type,
        chains = chains,
        return_raw_stan_output = return_raw_stan_output,
        prior_variance_treatment = prior_variance_treatment,
        stan_refresh = stan_refresh,
        ...
    )
}
