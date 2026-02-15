# bclogit

**Conditional Logistic Regression with Optional Concordant Pairs**

## Overview

`bclogit` is an R package for fitting conditional logistic regression models that can optionally incorporate information from concordant pairs (and a reservoir of additional controls) to improve estimation of discordant pairs. This approach effectively uses a Bayesian prior derived from the concordant data to inform the discordant pairs analysis, suitable for matched case-control studies where data sparsity might be an issue.

## Installation

You can install the development version of `bclogit` from source:

```r
# install.packages("devtools")
devtools::install()
```

or use the command line

```
R CMD INSTALL bclogit
```

## Usage

Here is a basic example of how to use `bclogit`:

```r
library(bclogit)

# Simulate some matched pair data
set.seed(123)
n_pairs <- 50
x1 <- rnorm(n_pairs * 2)
x2 <- rbinom(n_pairs * 2, 1, 0.5)
strata <- rep(1:n_pairs, each = 2)
treatment <- rep(c(0, 1), n_pairs) # Case-Control
# ... (simulate response y) ...

# Fit the model
# Using Formula Interface
fit <- bclogit(y ~ x1 + x2, data = mydata, treatment = treatment, strata = strata)

# Summary of results
summary(fit)

# Coefficients
coef(fit)

# Confidence Intervals
confint(fit)
```

## Features

- **Standard Formula Interface**: Works like `glm` or `lm`.
- **Bayesian Prior Integration**: Uses `rstan` to fit models with informative priors.
- **Multiple Prior Types**: Supports "naive", "G prior", "PMP", and "hybrid".
- **Concordant Pairs Method**: Supports "GLM", "GEE", and "GLMM" for the concordant pairs step.

## License

This package is licensed under the MIT License.
