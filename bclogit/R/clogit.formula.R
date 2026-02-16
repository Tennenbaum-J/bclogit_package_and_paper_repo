#' @export
#' @describeIn clogit Formula method
#' @examples
#' \dontrun{
#' n <- 200
#' dat <- data.frame(
#'   y = rbinom(n, 1, 0.5),
#'   x1 = rnorm(n),
#'   treatment = rep(c(0, 1), n / 2),
#'   strata = rep(1:(n / 2), each = 2)
#' )
#' fit <- clogit(y ~ x1, data = dat, treatment = treatment, strata = strata)
#' # Inference on treatment only (faster):
#' fit2 <- clogit(y ~ x1, data = dat, treatment = treatment, strata = strata,
#'                do_inference_on_var = 1)
#' }
clogit.formula <- function(formula, data, treatment = NULL, strata = NULL, 
                           subset = NULL, na.action = NULL,
                           do_inference_on_var = "all", ...) {
    assertClass(formula, "formula")
    assertFunction(na.action, null.ok = TRUE)
    assert(
        checkString(do_inference_on_var),
        checkCount(do_inference_on_var, positive = TRUE)
    )

    cl <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "drop.unused.levels", "treatment", "strata"), names(cl), 0L)
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

    full_cl <- match.call()
    if ("treatment" %in% names(full_cl)) {
        t_name <- deparse(full_cl$treatment)
    } else {
        t_name <- "treatment"
    }

    clogit.default(
        y = response,
        X = X,
        treatment = trt_vec,
        strata = strata_vec,
        treatment_name = t_name,
        call = full_cl,
        do_inference_on_var = do_inference_on_var,
        ...
    )
}
