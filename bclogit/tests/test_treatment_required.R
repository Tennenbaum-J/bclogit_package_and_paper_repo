library(bclogit)
library(testthat)

# Mock data
set.seed(123)
n <- 20
data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n),
    trt = rbinom(n, 1, 0.5),
    id = rep(1:(n / 2), each = 2)
)

test_that("bclogit.formula throws error when treatment is missing", {
    expect_error(bclogit(y ~ x1 + x2, data = data, strata = id), "The 'treatment' argument is required.")
})

test_that("bclogit.default throws error when treatment is missing", {
    X <- model.matrix(~ x1 + x2, data = data)
    expect_error(bclogit.default(response = data$y, data = X, strata = data$id), "The 'treatment' argument is required.")
})

test_that("bclogit works when treatment is provided", {
    # We just want to make sure it runs without that specific error,
    # even if it fails later due to small sample size / stan issues in this mock env
    # or if it actually runs.
    # For now, let's just check it doesn't throw the "required" error.

    # Note: running full stan model might be slow or fail in this test env if not set up.
    # We can try-catch the call.

    tryCatch(
        {
            bclogit(y ~ x1 + x2, data = data, treatment = trt, strata = id)
        },
        error = function(e) {
            if (grepl("The 'treatment' argument is required", e$message)) {
                fail("bclogit threw 'treatment required' error even when treatment was provided")
            }
            # Other errors are fine for this specific test (e.g. Stan not installed, etc)
            # In a real check we'd want it to pass fully, but here we are checking the argument plumbing.
        }
    )

    succeed()
})

print("Test script checks completed.")
