library(ggplot2)
library(splines)
library(MLmetrics)
library(Hmisc)
library(Matrix)
library(glmnet)

num_seeds <- 1
random_seeds <- sample(0:2000, num_seeds)


n <- 1000
x <- seq(0, 1, length.out = n)
f <- function(x) exp(-2 * x^2)
z <- f(x)
sd <- 0.3
sigma_sq <- sd^2
y <- z + rnorm(n, mean = 0, sd = sd)
mean(z^2)/mean((z-y)^2)
J <- 100


data <- data.frame(x = x, y = y, z = z)

ggplot(data, aes(x = x)) +
  geom_line(aes(y = z), color = "blue") +
  labs(x = "x", y = "y")

ggplot(data, aes(x = x)) +
  geom_point(aes(y = y), color = "lightblue") +
  geom_line(aes(y = z), color = "black") +
  labs(title = "Function vs Random Points", x = "x", y = "y") +
  theme(plot.title = element_text(size = 16)) +
  scale_color_manual(values = c("lightblue", "black")) +
  theme_minimal()


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



B <- fourier_basis(x, floor(J/2)-1)

Hs1 <- horseshoe(y, B, method.tau = c("halfCauchy"), Sigma2 = (sd)^2, tau=1,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
beta_hat1 <- as.vector(unlist(Hs1[1]))


# Create a data frame
plot_data <- data.frame(
  x = x,
  True_Function = z,
  Estimated = as.numeric(B%*%beta_hat1)
)

# Plot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = True_Function), color = "red", linewidth = 0.8) +  # True function (line)
  geom_point(aes(y = Estimated), color = "blue", size = 0.5, alpha = 0.5) +  # Estimates (points)
  labs(
    title = "True Function vs. Posterior Estimates",
    x = "x",
    y = "Value"
  ) +
  theme_minimal()



k <- cv.glmnet(B,z, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (1/((sigma_sq)/k) * diag(ncol(B)) + (1/sigma_sq) * crossprod(B))
theta <- solve(precision_matrix, (1/sigma_sq) * crossprod(B, y))

MSE(z,as.numeric(B%*%beta_hat1))
MSE(z, as.numeric(B%*%theta))


J_values <- c(seq(8,100, by = 4))
k_values <- rep(0,length(J_values))
mse_values_2_1 <- numeric(length(J_values))
mse_values_2_2 <- numeric(length(J_values))

for (i in seq_along(J_values)) {
  
  B <- fourier_basis(x, floor(J_values[i]/2)-1)
  Hs1 <- horseshoe(y, B, method.tau = c("halfCauchy"), Sigma2 = (sd)^2, tau=1,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
  
  beta_hat1 <- as.vector(unlist(Hs1[1]))
  
  k_values[i] <- cv.glmnet(B,z, alpha=0, nfolds = 20)$lambda.min
  
  precision_matrix <- (1/((sigma_sq)/k_values[i]) * diag(ncol(B)) + (1/sigma_sq) * crossprod(B))
  
  beta_hat2 <- solve(precision_matrix, (1/sigma_sq) * crossprod(B, y))
  
  mse_values_2_1[i] <- MSE(z, as.numeric(B %*% beta_hat1))
  
  mse_values_2_2[i] <- MSE(z, as.numeric(B %*% beta_hat2))
  
  cat(sprintf("J=%d | k=%.2f | MSE_1=%.4f, | MSE_2=%.4f\n", J_values[i], k_values[i], 
              
              mse_values_2_1[i], mse_values_2_2[i]))
}

results <- data.frame(
  J = J_values,
  k = k_values,
  MSE_1 = mse_values_2_1,
  MSE_2 = mse_values_2_2
)

ggplot(results, aes(x = J_values)) +
  geom_line(aes(y = MSE_1), color = "red", linewidth = 0.8) +  # True function (line)
  geom_line(aes(y = MSE_2), color = "blue", linewidth = 0.8) +  # True function (line)
  labs(
    title = "Hs vs Normal",
    x = "Knots",
    y = "MSE"
  ) +
  theme_minimal()


