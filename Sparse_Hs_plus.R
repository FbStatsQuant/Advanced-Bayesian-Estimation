# Load required libraries
library(MASS)       # For multivariate normal sampling
library(statmod)    # For inverse-Gaussian sampling
library(Matrix)

# Gibbs sampler settings
num_samples <- 5000  # Number of iterations
burn_in <- 1000      # Burn-in period
# n <- 2000             # Number of observations
# p <- 150              # Number of predictors
# 
# # Generate synthetic data
# set.seed(42)
# X <- matrix(rnorm(n * p), n, p)
# true_beta <- rep(0, p)  # Initialize all zeros
# nonzero_indices <- sample(1:p, size = round(0.1 * p), replace = FALSE)  # Select 10% nonzero indices
# true_beta[nonzero_indices] <- rnorm(length(nonzero_indices), mean = 2, sd = 10)  # Sample from N(0,5)
# sigma_true <- 1
# y <- X %*% true_beta + rnorm(n, 0, sigma_true)

X <- B
y <- y
n <- nrow(X)
p <- ncol(X)


# Initialize parameters
beta <- rep(0, p)
lambda_j <- rep(1, p)
phi_j <- rep(1, p)
tau <- 1
sigma2 <- 1  # Initial noise variance

# Storage for posterior samples
beta_samples <- matrix(0, num_samples, p)
tau_samples <- rep(0, num_samples)
sigma2_samples <- rep(0, num_samples)

# Gibbs Sampling Loop
for (i in 1:num_samples) {
  
  # Sample beta
  D_inv <- Diagonal(p, 1 / (lambda_j^2 * phi_j^2 * tau^2 * sigma2))
  Sigma_beta <- solve(t(X) %*% X + solve(D))
  mu_beta <- Sigma_beta %*% t(X) %*% y
  beta <- mvrnorm(1, mu_beta, Sigma_beta)  # Multivariate normal sample
  
  # Sample lambda_j
  z_j <- beta / (phi_j * tau * sqrt(sigma2))
  lambda_j <- 1 / sqrt(rinvgauss(p, mean = 1 / abs(z_j), shape = 1))
  
  # Sample phi_j
  z_j <- beta / (lambda_j * tau * sqrt(sigma2))
  phi_j <- 1 / sqrt(rinvgauss(p, mean = 1 / abs(z_j), shape = 1))
  
  # Sample tau
  xi <- rgamma(1, shape = 1, rate = 1)
  tau <- sqrt(rgamma(1, shape = (p + 1) / 2, 
                     rate = sum((beta / (lambda_j * phi_j * sqrt(sigma2)))^2) / 2 + 1 / xi))
  
  # Sample sigma^2 using inverse-chi-squared distribution (Jeffreys prior)
  residuals <- y - X %*% beta
  S <- sum(residuals^2) / n
  sigma2 <- 1 / rgamma(1, shape = n / 2, rate = S / 2)  # Inverse-Chi-Squared
  
  # Store samples
  beta_samples[i, ] <- beta
  tau_samples[i] <- tau
  sigma2_samples[i] <- sigma2
}

# Remove burn-in samples
beta_samples <- beta_samples[(burn_in + 1):num_samples, ]
tau_samples <- tau_samples[(burn_in + 1):num_samples]
sigma2_samples <- sigma2_samples[(burn_in + 1):num_samples]

# Plot posterior distributions of tau and sigma^2
par(mfrow = c(1, 2))
hist(tau_samples, breaks = 30, probability = TRUE, col = "lightblue", main = "Posterior of Tau (Global Shrinkage)")
hist(sigma2_samples, breaks = 30, probability = TRUE, col = "lightblue", main = "Posterior of Sigma^2 (Noise Variance)")

# Compute summary statistics of beta estimates
beta_summary <- data.frame(
  Mean = colMeans(beta_samples),
  Median = apply(beta_samples, 2, median),
  Std_Dev = apply(beta_samples, 2, sd),
  True_Beta = true_beta
)

# Print summary
print(beta_summary)

mean(abs(colMeans(beta_samples)-true_beta))
