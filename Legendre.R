library(horseshoe)
library(MLmetrics)
library(glmnet)

# Nonparametric Regression with Legendre Polynomial Basis (n=1000, J=20)
# ------------------------------------------------------------

# Generate synthetic data
set.seed(123)
n <- 1000
x <- seq(0, 1, length.out = n)

# True function (nonlinear example)
true_f <- function(x) {
  0.5 * sin(4 * pi * x) + 2 * exp(-20 * (x - 0.5)^2) + 3*sqrt(x)
}
f_true <- true_f(x)

# Add Gaussian noise
sigma <- 0.7
y <- f_true + rnorm(n, 0, sigma)

# Rescale x to [-1, 1] for Legendre polynomials
x_scaled <- 2 * x - 1

# Function to compute Legendre polynomials up to degree J-1
generate_legendre_basis <- function(x, J) {
  X <- matrix(0, nrow = length(x), ncol = J)
  
  # Recurrence relation for Legendre polynomials
  X[, 1] <- 1                     # P0(x) = 1
  if(J >= 2) X[, 2] <- x          # P1(x) = x
  
  for (k in 3:J) {
    n <- k - 2  # Degree index (0-based)
    # P_{n}(x) = [(2n-1)x P_{n-1}(x) - (n-1) P_{n-2}(x)] / n
    X[, k] <- ((2*n - 1) * x * X[, k-1] - (n - 1) * X[, k-2]) / n
  }
  
  return(X)
}

# Create design matrix with J=20 Legendre polynomials
J <- 200
X <- generate_legendre_basis(x_scaled, J)

# Bayesian regression with normal prior (ridge regression)
lambda <- cv.glmnet(X,f_true, alpha=0, nfolds = 20)$lambda.min
sigma_known <- sigma  # Assumed noise level

# Compute posterior mean coefficients
XtX <- t(X) %*% X
Xty <- t(X) %*% y
lambda_I <- diag(J) * (sigma_known^2 / lambda)
beta_post <- solve(XtX + lambda_I) %*% Xty

# Predictions
f_hat <- X %*% beta_post

# Compute MSE
mse <- mean((f_true - f_hat)^2)
cat("MSE:", round(mse, 5), "\n")

# Plot results
plot(x, y, col = "gray", main = "Legendre Polynomial Regression",
     xlab = "x", ylab = "y", cex = 0.5)
lines(x, f_true, col = "black", lwd = 2, lty = 2)
lines(x, f_hat, col = "red", lwd = 2)
legend("topright", legend = c("True Function", "Legendre Fit"),
       col = c("black", "red"), lwd = 2, lty = c(2, 1))



Hs1 <- horseshoe(y, X, method.tau = c("halfCauchy"), Sigma2 = (sigma)^2, tau=1,
                 method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
beta_hat1 <- as.vector(unlist(Hs1[1]))

# Predictions
f_hat <- X %*% beta_hat1

# Compute MSE
mse <- mean((f_true - f_hat)^2)
cat("MSE:", round(mse, 5), "\n")

# Plot results
plot(x, y, col = "gray", main = "Legendre Polynomial Regression",
     xlab = "x", ylab = "y", cex = 0.5)
lines(x, f_true, col = "black", lwd = 2, lty = 2)
lines(x, f_hat, col = "red", lwd = 2)
legend("topright", legend = c("True Function", "Legendre Fit"),
       col = c("black", "red"), lwd = 2, lty = c(2, 1))

beta_hat1

