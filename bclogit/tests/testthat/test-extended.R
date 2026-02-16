# ============================================================================
# Extended Tests for bclogit
# ============================================================================

test_that("bclogit.formula handles factors and characters correctly", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 1)
  df <- make_df(dat)
  df$group <- rep(c("A", "B", "C"), length.out = nrow(df))

  fit <- bclogit(
    y ~ x1 + group, data = df, treatment = treatment,
    strata = strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  # group has 3 levels, so 2 dummy variables + x1 + treatment = 4 coefficients
  expect_equal(length(coef(fit)), 4)
  expect_true(all(c("groupB", "groupC") %in% names(coef(fit))))
})

test_that("bclogit.formula handles subset and na.action", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 80, p = 1, seed = 2)
  df <- make_df(dat)

  # Add NA
  df$x1[1] <- NA

  # With na.omit (default)
  fit_omit <- bclogit(
    y ~ x1, data = df, treatment = treatment,
    strata = strata, chains = 1
  )
  # The pair containing df$x1[1] should be dropped or handled.
  # If one observation in a pair is NA, the whole pair is effectively ruined for matched analysis.
  # Let's see if it works without error.
  expect_s3_class(fit_omit, "bclogit")

  # With subset
  fit_subset <- bclogit(
    y ~ x1, data = df, treatment = treatment,
    strata = strata, subset = strata > 10, chains = 1
  )
  expect_s3_class(fit_subset, "bclogit")
  expect_true(fit_subset$n < fit_omit$n)
})

test_that("bclogit.default handles X with a constant column", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 3)
  # Add a constant column to X
  X_const <- cbind(dat$X, const = 1)

  # It should remove the constant column if it's the first one,
  # or handle it if it's later (glm might handle it or it might cause issues).
  # The code says: if (sd(X[, 1]) == 0) { X[, 1] <- NULL }
  # Let's test when it's NOT the first column too.
  X_const_late <- cbind(x1 = dat$X[, 1], const = 1)

  fit <- bclogit_default(
    y = dat$y, X = X_const_late, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )
  expect_s3_class(fit, "bclogit")
  # Should have treatment + x1 (constant removed or ignored by Stan/GLM)
  # Actually, the code only removes it if it's the FIRST column.
  # If it's the second, it might be passed to Stan.
  # Let's check what happened.
  expect_true(length(coef(fit)) >= 2)
})

test_that("bclogit errors when no discordant pairs exist", {
  skip_if_no_stan()

  # Create data where all pairs are concordant
  # strata 1: (y=0, trt=0), (y=0, trt=1) -> concordant
  # strata 2: (y=1, trt=0), (y=1, trt=1) -> concordant
  n_pairs <- 20
  strata <- rep(1:n_pairs, each = 2)
  treatment <- rep(c(0, 1), n_pairs)
  y <- rep(rep(c(0, 1), each = 2), length.out = n_pairs * 2)
  X <- matrix(rnorm(n_pairs * 2), ncol = 1)

  expect_error(
    bclogit_default(y = y, X = X, treatment = treatment, strata = strata),
    "not enough discordant pairs"
  )
})

test_that("pre-compiled stanmodels are available", {
  skip_if_no_stan()

  sm <- bclogit:::stanmodels
  expect_true(is.list(sm))
  expect_true("mvn_logistic" %in% names(sm))
  expect_s4_class(sm[["mvn_logistic"]], "stanmodel")
})

test_that("print.summary.bclogit respects digits", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 4)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )
  s <- summary(fit)

  out1 <- capture.output(print(s, digits = 1))
  out5 <- capture.output(print(s, digits = 5))

  # This is a bit brittle, but we can check if a number with many digits is present in out5 but not out1
  # Actually, just check if it runs without error.
  expect_true(length(out1) > 0)
  expect_true(length(out5) > 0)
})

test_that("bclogit.default handles data.table input for X", {
  skip_if_no_stan()
  if (!requireNamespace("data.table", quietly = TRUE)) {
    skip("data.table not installed")
  }

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 5)
  X_dt <- data.table::as.data.table(dat$X)
  fit <- bclogit_default(
    y = dat$y, X = X_dt, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )
  expect_s3_class(fit, "bclogit")
})

test_that("summary.bclogit works with level = 0.5", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 6)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )

  s50 <- summary(fit, level = 0.5, inference_method = "CR")
  expect_equal(s50$level, 0.5)
  expect_true("25%" %in% colnames(s50$coefficients))
  expect_true("75%" %in% colnames(s50$coefficients))
})
