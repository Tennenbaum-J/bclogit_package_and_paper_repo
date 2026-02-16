# ============================================================================
# More Tests for bclogit
# ============================================================================

test_that("bclogit.formula handles interaction terms", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 2, seed = 101)
  df <- make_df(dat)

  fit <- bclogit(
    y ~ x1 * x2, data = df, treatment = treatment,
    strata = strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  # treatment, x1, x2, x1:x2 = 4 coefficients
  expect_equal(length(coef(fit)), 4)
  expect_true("x1:x2" %in% names(coef(fit)))
})

test_that("bclogit handles many covariates", {
  skip_if_no_stan()

  p <- 10
  dat <- generate_test_data(n_pairs = 100, p = p, seed = 102)
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  expect_equal(length(coef(fit)), p + 1)
})

test_that("bclogit.default handles multiple constant columns in X", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 103)
  X_const <- cbind(dat$X, c1 = 1, c2 = 2)

  # The first column is not constant, but c1 and c2 are.
  # bclogit.default currently only removes the FIRST column if it's constant.
  # So c1 and c2 will be passed to model.matrix or as is.
  # Stan might have issues if they are exactly collinear, but GLM might drop them.
  # Let's see how it behaves.
  fit <- bclogit_default(
    y = dat$y, X = X_const, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  expect_true(length(coef(fit)) >= 2)
})

test_that("bclogit handles collinear columns in X", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 104)
  X_coll <- cbind(dat$X, x2 = dat$X[, 1] * 2)

  fit <- bclogit_default(
    y = dat$y, X = X_coll, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  # Collinearity might lead to NAs in concordant model,
  # but bclogit.default tries to handle them.
  expect_true(!is.null(coef(fit)))
})

test_that("bclogit respects chains argument", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 40, p = 1, seed = 105)
  # Use 2 chains
  fit <- bclogit_default(
    y = dat$y, X = dat$X, treatment = dat$treatment,
    strata = dat$strata, chains = 2
  )

  expect_equal(fit$model@sim$chains, 2)
})

test_that("bclogit handles factors in strata and treatment", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 106)
  df <- make_df(dat)
  df$treatment <- factor(df$treatment, levels = c(0, 1))
  df$strata <- factor(df$strata)

  # bclogit.formula should handle these
  fit <- bclogit(
    y ~ x1, data = df, treatment = treatment,
    strata = strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  expect_equal(length(coef(fit)), 2)
})

test_that("bclogit.default handles X as data.frame with factors", {
  skip_if_no_stan()

  dat <- generate_test_data(n_pairs = 60, p = 1, seed = 107)
  X_df <- data.frame(x1 = dat$X[, 1], f = factor(rep(c("A", "B"), length.out = nrow(dat$X))))

  fit <- bclogit_default(
    y = dat$y, X = X_df, treatment = dat$treatment,
    strata = dat$strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  # treatment, x1, fB = 3 coefficients
  expect_equal(length(coef(fit)), 3)
})

test_that("bclogit handles very few discordant pairs (just enough)", {
  skip_if_no_stan()

  # ncol(X) + 2 = 1 + 2 = 3 discordant pairs needed
  n_pairs <- 20
  strata <- rep(1:n_pairs, each = 2)
  treatment <- rep(c(0, 1), n_pairs)
  # Make 3 discordant pairs, rest concordant
  y <- rep(0, n_pairs * 2)
  # Discordant pairs: 1-3
  for (i in 1:3) {
    y[2*i] <- 1 # (0, 1) -> discordant
  }
  # Concordant pairs: 4-20
  # already (0, 0)

  X <- matrix(rnorm(n_pairs * 2), ncol = 1)

  fit <- bclogit_default(
    y = y, X = X, treatment = treatment, strata = strata, chains = 1
  )

  expect_s3_class(fit, "bclogit")
  expect_equal(fit$num_discordant, 3)
})

test_that("bclogit handles missing values in treatment or strata", {
  dat <- generate_test_data(n_pairs = 20, p = 1, seed = 108)
  df <- make_df(dat)

  df$treatment[1] <- NA
  expect_error(
    bclogit(y ~ x1, data = df, treatment = treatment, strata = strata, na.action = na.fail),
    "missing values"
  )

  df <- make_df(dat)
  df$strata[1] <- NA
  expect_error(
    bclogit(y ~ x1, data = df, treatment = treatment, strata = strata, na.action = na.fail),
    "missing values"
  )
})
