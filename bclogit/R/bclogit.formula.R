#' @param formula A symbolic description of the model to be fitted.
#' @param data An optional data frame containing the variables in the model.
#' @param treatment Optional vector specifying the treatment variable.
#' @param strata Vector specifying the strata (matched pairs).
#' @param ... Additional arguments passed to the default method.
#' @export
#' @describeIn bclogit Formula method
bclogit.formula <- function(formula, data, treatment, strata = NULL, ...) {
    cl <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "drop.unused.levels"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(stats::model.frame)
    mf <- eval(cl, parent.frame())

    response <- model.response(mf, "numeric")
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)

    if (colnames(X)[1] == "(Intercept)") {
        X <- X[, -1, drop = FALSE]
    }

    if (missing(data)) data <- environment(formula)

    if (missing(treatment)) stop("The 'treatment' argument is required.")

    trt_vec <- eval(substitute(treatment), data, parent.frame())
    strata_vec <- eval(substitute(strata), data, parent.frame())

    # Extract treatment name from call
    # We use the full match.call() which has all arguments
    full_cl <- match.call()
    if ("treatment" %in% names(full_cl)) {
        t_name <- deparse(full_cl$treatment)
        # Clean up if it's long or complex? usually deparse is fine.
    } else {
        t_name <- "treatment"
    }

    bclogit.default(response = response, data = X, treatment = trt_vec, strata = strata_vec, treatment_name = t_name, ...)
}
