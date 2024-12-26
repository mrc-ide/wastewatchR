// logistic_regression_probabilities.stan
data {
  int<lower=1> N;              // number of data points
  real x[N];                   // Infections
  real<lower=0, upper=1> y[N]; // Observed probabilities
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  // Priors (you can adjust based on domain knowledge)
  alpha ~ normal(0, 5);
  beta  ~ normal(0, 5);
  sigma ~ exponential(1);  // positivity-constrained scale

  // Likelihood
  for (n in 1:N) {
    real mu = inv_logit(alpha + beta * x[n]);  // logistic curve
    y[n] ~ normal(mu, sigma);
  }
}
