data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0, upper=1> y[N]; // Binary outcome
  vector[N] xw;                // Parameter of interest (orthogonalized)
  matrix[N, P] X;              // Nuisance predictors

  // Informative Prior from Dataset A
  vector[P] mu_A;
  matrix[P, P] Sigma_A;
}

transformed data {
  matrix[P, P] L_A = cholesky_decompose(Sigma_A);
}

parameters {
  real beta_w;
  vector[P] z;                 // Helper: Standard Normal(0,1) for nuisance
  real<lower=0> g;             // The hyper-parameter for nuisance strength
}

transformed parameters {
  // Non-centered parameterization: beta_nuis = mu_A + sqrt(g) * L_A * z
  // Mathematically equivalent to beta_nuis ~ N(mu_A, g * Sigma_A)
  vector[P] beta_nuis = mu_A + sqrt(g) * (L_A * z);
}

model {
  g ~ inv_gamma(0.5, N / 2.0);
  z ~ std_normal();

  // 1. Calculate the PMP Adjustment (log square root of Fisher Info)
  vector[N] eta = xw * beta_w + X * beta_nuis;
  vector[N] p = inv_logit(eta);
  real I_ww = 0;
  for (i in 1:N) {
    // Weighting term for logistic information: p*(1-p)
    I_ww += square(xw[i]) * p[i] * (1 - p[i]);
  }

  // Add log(sqrt(I_ww)) to the target log-density
  target += 0.5 * log(I_ww);

  // 2. Likelihood
  y ~ bernoulli_logit(eta);
}
