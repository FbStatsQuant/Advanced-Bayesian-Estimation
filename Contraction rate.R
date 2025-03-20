library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)

n_values <- seq(1000, 20000, by = 1000)
MSE <- rep(0, length(n_values))
MSE_r <- rep(0, length(n_values))
b <- 0.7

s <- function(m, x) {
  sqrt(2) * sum(seq(m)^(-1.5) * sin(seq(m)) * cos(((seq(m)) - 0.5) * pi * x))
}

sd <- 0.3

for (i in seq_along(n_values)) {
  n <- n_values[i]
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = floor(n^b)+1),1.0001, 1.0002)
  Sp_order <- 3
  B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  Hs <- horseshoesp(y, B_s, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                    Sigma2 = sd^2)
  
  f_hat_hs <- as.numeric(B_s %*% as.vector(unlist(Hs$BetaHat)))
  MSE[i] <- mean((z - f_hat_hs)^2)
}

print(MSE)

library(ggplot2)

# Create data frame
my_data <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE)
)

# Fit linear model
fit_hs_bs <- lm(log_MSE ~ log_n, data = my_data)
summary_fit <- summary(fit_hs_bs)

# Extract slope
slope <- round(coef(fit_hs_bs)[2], 3)

# Create plot
ggplot(my_data, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): B-splines / Horseshoe",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data$log_n) + 0.1, y = min(my_data$log_MSE) + 0.1,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()





fourier_basis <- function(x, K) {
  
  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
  
  for (m in 1:K) {
    basis[, 2 * m]   <- sin(2 * pi * m * x)
    basis[, 2 * m + 1] <- cos(2 * pi * m * x)
  }
  
  return(basis)
}

for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- 2*floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  B_f <- fourier_basis(x, floor(J/2))
  y <- z + rnorm(n, mean = 0, sd = sd)
  lambda_s <- cv.glmnet(B_f,z, alpha=0, nfolds = 10)$lambda.min
  lambda_s <- 0.5
  precision_matrix <- (lambda_s/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
  theta_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
  f_hat_n <- B_f%*%theta_n
  MSE_r[i] <- mean((z - f_hat_n)^2)
}

print(MSE_r)
 
my_data <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_r)
)


fit <- lm(log_MSE ~ log_n, data = my_data)
summary(fit)

ggplot(my_data, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n) for ridge regression under best lambda",
       x = "log(n)", y = "log(MSE)") +
  theme_minimal()





my_data <- data.frame(
  log_n = log(n_values), # New predictor log(log(n))
  log_MSE = log(MSE)
)

# Fit multiple regression model
fit_hs_bs <- lm(log_MSE ~ log_n, data = my_data)
summary(fit_hs_bs)



for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- 2*floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  x_scaled <- 2 * x - 1
  B_l <- matrix(0, nrow = length(x), ncol = J+1)
  B_l[, 1] <- 1
  if(J >= 2) B_l[, 2] <- x_scaled
  for (k in 3:(J+1)) {
    m <- k - 2  
    B_l[, k] <- ((2*m - 1) * x_scaled * B_l[, k-1] - (m - 1) * B_l[, k-2]) / m
  }
  z <- sapply(x, function(xi) s(1000, xi))
  y <- z + rnorm(n, mean = 0, sd = sd)
  lambda_s <- cv.glmnet(B_l,z, alpha=0, nfolds = 5)$lambda.min
  #lambda_s <- 0.05
  precision_matrix <- (lambda_s/((sd^2)) * diag(ncol(B_l)) + (1/(sd^2)) * crossprod(B_l))
  theta_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_l, y))
  f_hat_n <- B_l%*%theta_n
  MSE_r[i] <- mean((z - f_hat_n)^2)
}

print(MSE_r)

my_data <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_r)
)


fit <- lm(log_MSE ~ log_n, data = my_data)
summary(fit)

ggplot(my_data, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n) for ridge regression under best lambda",
       x = "log(n)", y = "log(MSE)") +
  theme_minimal()




