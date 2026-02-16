#' Print summary of a bclogit model
#'
#' @param x A `summary.bclogit` object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments.
#' @export
print.summary.bclogit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    cat("Discordant Pairs Bayesian Model\n")
    cat("-------------------------------\n")
    cat("Discordant pairs: ", x$num_discordant, "\n")
    cat("Concordant pairs: ", x$num_concordant, " (used for prior)\n")

    if ("Rhat" %in% colnames(x$coefficients)) {
        cat("\nConvergence Diagnostics:\n")
        diag_cols <- c("Rhat", "n_eff")
        print(x$coefficients[, diag_cols, drop = FALSE], digits = digits)
        cat("Note: R-hat should be < 1.05.\n")
    }

    inference_desc <- switch(x$inference_method,
        "HPD_one" = "unimodal HPD interval (with pval binary search 20 iters)",
        "HPD_many" = "HPD interval with disjoint support via ggdist (with pval binary search 50 iters)",
        "CR" = "equal-tailed credible region (with pval binary search 20 iters)",
        "unknown method"
    )
    cat(sprintf("\nInterval method: %s.\n", inference_desc))

    cat("\nCoefficients:\n")
    
    # Select columns to display: exclude diagnostics (Rhat, n_eff) but keep everything else.
    # Pr(!=0) is always the last column.
    exclude_cols <- c("Rhat", "n_eff")
    all_cols <- colnames(x$coefficients)
    main_cols <- which(!all_cols %in% exclude_cols)
    main_mat <- x$coefficients[, main_cols, drop = FALSE]
    
    nc <- ncol(main_mat)
    # cs.ind and tst.ind must cover columns 1:(nc-1); the p-value column (nc)
    # is handled implicitly by has.Pvalue = TRUE
    printCoefmat(main_mat,
        digits = digits, signif.stars = TRUE,
        P.values = TRUE, has.Pvalue = TRUE,
        cs.ind = seq_len(nc - 1), tst.ind = integer(0)
    )

    if (!is.null(x$prior_info$mu)) {
        cat("\nPrior Estimates (from concordant data):\n")

        prior_mu <- x$prior_info$mu
        prior_sigma <- x$prior_info$Sigma
        
        if (length(prior_mu) == nrow(x$coefficients)) {
            p_names <- rownames(x$coefficients)
            # Remove the treatment variable (strata variable) which is the first element
            if (length(prior_mu) > 1) {
                p_mu <- prior_mu[-1]
                p_names <- p_names[-1]
                
                if (!is.null(prior_sigma)) {
                    p_se <- sqrt(diag(prior_sigma))[-1]
                    # Create a data frame for nice alignment
                    p_df <- data.frame(
                        Estimate = p_mu,
                        `Std. Error` = p_se,
                        row.names = p_names,
                        check.names = FALSE
                    )
                    print(p_df, digits = digits)
                } else {
                    names(p_mu) <- p_names
                    print(p_mu, digits = digits)
                }
            } else {
                cat("No covariate priors to display.\n")
            }
        } else {
            cat("Mean vector length:", length(prior_mu), "\n")
            print(prior_mu, digits = digits)
        }
    }

    invisible(x)
}
