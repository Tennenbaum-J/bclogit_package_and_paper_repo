data {
  int<lower=1> N;                      // number of observations
  int<lower=1> K;                      // number of predictors (treatment + covariates)
  matrix[N, K] X;                      // design matrix
  int<lower=0,upper=1> y[N];           // binary outcomes
  vector[K] mu;                        // prior mean vector
  matrix[K, K] Sigma;                  // prior covariance matrix
}

transformed data {
  // Separate treatment (index 1) from nuisance (indices 2:K)
  real mu_trt = mu[1];
  real sigma_trt = sqrt(Sigma[1, 1]);

  // Extract nuisance block
  int Kn = K - 1;
  vector[Kn] mu_nuis = mu[2:K];
  matrix[Kn, Kn] Sigma_nuis = Sigma[2:K, 2:K];
  matrix[Kn, Kn] L_nuis = cholesky_decompose(Sigma_nuis);
}

parameters {
  real beta_trt;                       // treatment coefficient (fixed prior, not scaled by g)
  vector[Kn] z;                        // Helper: Standard Normal(0,1) for nuisance
  real<lower=0> g;                     // The hyper-parameter g (scales nuisance only)
}

transformed parameters {
  // Assemble full beta vector
  vector[K] beta;
  beta[1] = beta_trt;
  // Nuisance: beta_nuis = mu_nuis + sqrt(g) * L_nuis * z
  // This is mathematically beta_nuis ~ N(mu_nuis, g * Sigma_nuis)
  beta[2:K] = mu_nuis + sqrt(g) * (L_nuis * z);
}

model {
  // Treatment: fixed prior (not scaled by g)
  beta_trt ~ normal(mu_trt, sigma_trt);
  // Hyper-prior: Zellner-Siow (applied to nuisance only)
  g ~ inv_gamma(0.5, N / 2.0);
  // Prior for the helper (leads to the g-prior for nuisance beta)
  z ~ std_normal();
  // Likelihood
  y ~ bernoulli_logit(X * beta);
}
