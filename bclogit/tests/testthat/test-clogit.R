test_that("clogit matches survival::clogit", {
  skip_if_not_installed("survival")
  library(survival)
  
  set.seed(42)
  n <- 200
  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n),
    treatment = rep(c(0, 1), n / 2),
    strata = rep(1:(n / 2), each = 2)
  )
  
  fit_surv <- survival::clogit(y ~ treatment + x1 + x2 + strata(strata), data = dat)
  fit_bcl  <- bclogit::clogit(y ~ x1 + x2, data = dat, treatment = treatment, strata = strata)
  
  expect_equal(unname(coef(fit_bcl)), unname(coef(fit_surv)), tolerance = 1e-5)
  # fast_logistic_regression currently only returns diagonal vcov (SEs)
  expect_equal(unname(diag(vcov(fit_bcl))), unname(diag(vcov(fit_surv))), tolerance = 1e-5)
})

test_that("clogit summary and print work", {
  set.seed(42)
  n <- 100
  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    treatment = rep(c(0, 1), n / 2),
    strata = rep(1:(n / 2), each = 2)
  )
  
  fit <- bclogit::clogit(y ~ x1, data = dat, treatment = treatment, strata = strata)
  s <- summary(fit)
  expect_s3_class(s, "summary.clogit_bclogit")
  expect_output(print(s), "Conditional Logistic Regression")
})

test_that("clogit works with do_inference_on_var", {
  set.seed(42)
  n <- 200
  dat <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n),
    treatment = rep(c(0, 1), n / 2),
    strata = rep(1:(n / 2), each = 2)
  )
  
  # Inference on treatment (variable 1 in wX)
  fit <- bclogit::clogit(y ~ x1 + x2, data = dat, treatment = treatment, strata = strata, do_inference_on_var = 1)
  expect_equal(length(coef(fit)), 3)
  expect_true(!is.na(fit$se[1]))
})
