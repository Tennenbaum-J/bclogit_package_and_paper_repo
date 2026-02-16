#' Summarize a clogit_bclogit model
#'
#' @param object A \code{clogit_bclogit} object.
#' @param ... Additional arguments (not used).
#' @return A \code{summary.clogit_bclogit} object.
#' @export
summary.clogit_bclogit <- function(object, ...) {
    assertClass(object, "clogit_bclogit")
    coefs <- object$coefficients
    se <- object$se
    z_val <- object$z
    p_val <- object$pval

    coef_mat <- cbind(
        Estimate = coefs,
        `Std. Error` = se,
        `z value` = z_val,
        `Pr(>|z|)` = p_val
    )
    rownames(coef_mat) <- names(coefs)

    res <- list(
        call = object$call,
        coefficients = coef_mat,
        num_discordant = object$num_discordant,
        num_concordant = object$num_concordant,
        n = object$n,
        treatment_name = object$treatment_name,
        do_inference_on_var = object$do_inference_on_var
    )

    class(res) <- "summary.clogit_bclogit"
    res
}
