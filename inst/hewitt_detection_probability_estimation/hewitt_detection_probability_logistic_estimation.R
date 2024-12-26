library(rstan)

# Observed data
data <- read.csv("inst/hewitt_detection_probability_estimation/hewitt_detection_probability_data.csv")

Infections <- data$Infections
Probability <- data$Probability
N <- length(Infections)

# Put into a list for Stan
stan_data <- list(
  N = N,
  x = Infections,
  y = Probability
)

# Compile & sample
fit <- stan(
  file = "inst/hewitt_detection_probability_estimation/hewitt_detection_stan_model.stan",
  data = stan_data,
  iter = 2000,      # total iterations
  warmup = 1000,    # burn-in
  chains = 4,
  cores = 2,
  seed = 123
)

# Print a summary of posterior estimates
print(fit, pars = c("alpha", "beta", "sigma"))

posterior_draws <- rstan::extract(fit, pars = c("alpha", "beta"))
alpha_mean <- mean(posterior_draws$alpha)
beta_mean  <- mean(posterior_draws$beta)

check <- plogis(alpha_mean + beta_mean * Infections)
plot(Infections, check, type = "l")
points(Infections, Probability, pch = 20)
