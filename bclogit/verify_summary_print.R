source("R/summary.bclogit.R")
source("R/print.summary.bclogit.R")

# Mock a bclogit object
mock_coefs <- c(treatment = 0.5, x1 = 0.2, x2 = -0.1)
mock_bclogit <- list(
    coefficients = mock_coefs,
    call = quote(bclogit(y ~ x1 + x2, data = data, treatment = trt)),
    num_discordant = 100,
    num_concordant = 50,
    prior_info = list(mu = c(0, 0, 0), Sigma = diag(3)),
    treatment_name = "treatment",
    model = NULL # simulating no stan model for basic check, or we can mock summary behavior
)
class(mock_bclogit) <- "bclogit"

print("Testing summary.bclogit with no model...")
summ <- summary.bclogit(mock_bclogit)
print(class(summ))

# Now let's try to mock the internal behavior of summary.bclogit regarding stan extraction
# valid summary.bclogit calls rstan::extract and rstan::summary.
# We can't easily mock rstan namespace calls without pkgload or mockery.
# But we can verify that print.summary.bclogit works on a constructed summary object.

print("Testing print.summary.bclogit with constructed object...")
# Construct a summary object manually that mimics what summary.bclogit returns
coef_mat <- matrix(
    c(
        0.5, 0.1, 0.3, 0.7, 0.99, 1.001, 1000,
        0.2, 0.1, 0.0, 0.4, 0.95, 1.000, 1000,
        -0.1, 0.1, -0.3, 0.1, 0.10, 1.002, 1000
    ),
    nrow = 3, byrow = TRUE
)
colnames(coef_mat) <- c("Estimate", "Est.Error", "L95%", "U95%", "Pr(>0)", "Rhat", "n_eff")
rownames(coef_mat) <- names(mock_coefs)

mock_summary <- list(
    call = mock_bclogit$call,
    coefficients = coef_mat,
    num_discordant = 100,
    num_concordant = 50,
    prior_info = mock_bclogit$prior_info
)
class(mock_summary) <- "summary.bclogit"

out <- capture.output(print.summary.bclogit(mock_summary))
writeLines(out, "summary_clean.txt")
print("Verification script finished.")
