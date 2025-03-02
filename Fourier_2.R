# 1. Define a function to generate Fourier basis functions
fourier_basis <- function(x, K) {
  # Check that x is within [0, 1]
  if (any(x < 0 | x > 1)) stop("x must be in [0,1]")
  
  n <- length(x)
  # Create a matrix that starts with a constant (intercept) column.
  # Then, for each k=1,...,K, add sin and cos terms.
  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
  
  for (k in 1:K) {
    basis[, 2 * k]   <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }
  
  return(basis)
}

# 2. Simulate data
set.seed(42)
n <- 200
x <- seq(0, 1, length.out = n)
# Define a true function that is a sum of sine waves (periodic)
f_true <- function(x) sin(2 * pi * x) + 0.5 * sin(4 * pi * x)
# Generate noisy observations from the true function
y <- f_true(x) + rnorm(n, mean = 0, sd = 0.3)

# 3. Generate Fourier basis expansion and fit the regression model
K <- 3  # Number of Fourier terms to use
X <- fourier_basis(x, K)
# Fit a linear model using the Fourier basis.
# We remove the intercept (using "-1") since our basis already includes a constant.
fit <- lm(y ~ X - 1)

# 4. Make predictions using the fitted model
y_pred <- predict(fit)

# 5. Plot the results
plot(x, y, pch = 16, col = "gray",
     main = "Nonparametric Regression using Fourier Basis",
     xlab = "x", ylab = "y")
lines(x, f_true(x), col = "black", lwd = 2, lty = 2)  # True function (dashed)
lines(x, y_pred, col = "blue", lwd = 2)  # Fitted Fourier regression curve
legend("topright", legend = c("True Function", "Fourier Fit"),
       col = c("black", "blue"), lty = c(2, 1), lwd = 2)
