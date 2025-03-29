library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)

set.seed(1991)

n <- 2000
b <- 0.6
J <- floor(n^b)
tau <- n^(-b-0.5)
x <- seq(0,1, length.out = n)

s <- function(m, x) {
  sqrt(2) * sum(seq(m)^(-1.5) * sin(seq(m)) * cos(((seq(m)) - 0.5) * pi * x))
}

z <- sapply(x, function(xi) s(1000, xi))
sd <- 0.3
y <- z + rnorm(n, mean = 0, sd = sd)

sum(z^2) / sum((z-y)^2)

data <- data.frame(x=x, z=z, y=y)

ggplot(data, aes(x=x)) +
  geom_point(aes(y=y), color = "lightblue", size = 1, alpha = 1) +
  geom_line(aes(y=z), color = "black", linewidth = 0.8) +
  labs(
    title = "True functions and random points",
    x = "x",
    y = "value"
  ) +
  theme_minimal()



# Bsplines



t <- c(-0.0002, -0.0001, seq(0,1, length.out = J+1),1.0001, 1.0002)
Sp_order <- 3

## Horseshoe

B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = FALSE)

Hs_s <- horseshoe(y, B_s, method.tau = c("fixed"), method.sigma = c("fixed"), tau = tau,
                  Sigma2 = sd^2)

theta_s_hs <- unlist(Hs_s$BetaHat)
f_hat_s_hs <- B_s%*%theta_s_hs

## Ridge

lambda_s <- cv.glmnet(B_s,z, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_s/((sd^2)) * diag(ncol(B_s)) + (1/(sd^2)) * crossprod(B_s))
theta_s_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_s, y))
f_hat_s_n <- B_s%*%theta_s_n

## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_s_hs, ridge = f_hat_s_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe Prior"), size = 1, alpha = 1) +  
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1, alpha = 1) +    
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +        
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (B-Splines)",
    x = "x",
    y = "Value",
    caption = "Black: True Function | Light Blue: Horseshoe Prior | Blue: Ridge Regression"
  ) +
  scale_color_manual(values = c("True Function" = "black", 
                                "Horseshoe Prior" = "lightblue", 
                                "Ridge Regression" = "blue")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0.5, vjust = 2, size = 10, face = "italic")  # Centered caption below the plot
  )






## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_s_hs, ridge = f_hat_s_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe Prior"), size = 1, alpha = 1) +  
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +        
  labs(
    x = "x",
    y = "Value"
  ) +
  scale_color_manual(values = c("True Function" = "black", 
                                "Horseshoe Prior" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0.5, vjust = 2, size = 10, face = "italic")  # Centered caption below the plot
  )


# Fourier



fourier_basis <- function(x, K) {
  
  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
  
  for (k in 1:K) {
    basis[, 2 * k]   <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }
  
  return(basis)
}

B_f <- fourier_basis(x, floor(J/2))

## Hs

Hs_f <- horseshoe(y, B_f, method.tau = c("fixed"), method.sigma = c("fixed"), tau = tau,
                  Sigma2 = sd^2)
theta_f_hs <- unlist(Hs_f$BetaHat)
f_hat_f_hs <- B_f%*%theta_f_hs

## Ridge

lambda_f <- cv.glmnet(B_f,z, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_f/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
theta_f_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
f_hat_f_n <- B_f%*%theta_f_n

## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_f_hs, ridge = f_hat_f_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe Prior"), size = 1, alpha = 1) +  
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1, alpha = 1) +    
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +        
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (Fourier)",
    x = "x",
    y = "Value",
    caption = "Black: True Function | Light Blue: Horseshoe Prior | Blue: Ridge Regression"
  ) +
  scale_color_manual(values = c("True Function" = "black", 
                                "Horseshoe Prior" = "lightblue", 
                                "Ridge Regression" = "blue")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0.5, vjust = 2, size = 10, face = "italic")  # Centered caption below the plot
  )



# Legendre



x_scaled <- 2 * x - 1

generate_legendre_basis <- function(x, J) {
  X <- matrix(0, nrow = length(x), ncol = J)
  
  X[, 1] <- 1                     # P0(x) = 1
  if(J >= 2) X[, 2] <- x          # P1(x) = x
  
  for (k in 3:J) {
    n <- k - 2  
    # P_{n}(x) = [(2n-1)x P_{n-1}(x) - (n-1) P_{n-2}(x)] / n
    X[, k] <- ((2*n - 1) * x * X[, k-1] - (n - 1) * X[, k-2]) / n
  }
  
  return(X)
}

B_l <- generate_legendre_basis(x_scaled, J+1)

## Horseshoe

Hs_l <- horseshoe(y, B_l, method.tau = c("fixed"), method.sigma = c("fixed"), tau=tau,
                  Sigma2 = sd^2)
theta_l_hs <- unlist(Hs_l$BetaHat)
f_hat_l_hs <- B_l%*%theta_l_hs

## Ridge

lambda_l <- cv.glmnet(B_l,z, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_l/((sd^2)) * diag(ncol(B_l)) + (1/(sd^2)) * crossprod(B_l))
theta_l_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_l, y))
f_hat_l_n <- B_l%*%theta_l_n

## Plot 

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_l_hs, ridge = f_hat_l_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe Prior"), size = 1, alpha = 1) +  
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1, alpha = 1) +    
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +        
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (Legendre)",
    x = "x",
    y = "Value",
    caption = "Black: True Function | Light Blue: Horseshoe Prior | Blue: Ridge Regression"
  ) +
  scale_color_manual(values = c("True Function" = "black", 
                                "Horseshoe Prior" = "lightblue", 
                                "Ridge Regression" = "blue")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0.5, vjust = 2, size = 10, face = "italic")  # Centered caption below the plot
  )


subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x[floor(n*0.05):floor(n*0.95)], true = z[floor(n*0.05):floor(n*0.95)], 
                   hs_s = f_hat_s_hs[floor(n*0.05):floor(n*0.95)], hs_f = f_hat_f_hs[floor(n*0.05):floor(n*0.95)], 
                   hs_l = f_hat_l_hs[floor(n*0.05):floor(n*0.95)])

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs_s, color = "Horseshoe Splines"), size = 1, alpha = 1) +  
  geom_point(aes(y = hs_l, color = "Horseshoe Fourier"), size = 1, alpha = 1) +
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +        
  labs(
    title = "True Function, Horseshoe, n = 5000 and J = 100",
    x = "x",
    y = "Value"
  ) +
  scale_color_manual(values = c("True Function" = "black", 
                                "Horseshoe Splines" = "blue", 
                                "Horseshoe Fourier" = "red")) +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0.5, vjust = 2, size = 10, face = "italic")  # Centered caption below the plot
  )



# MSE

mse_s_hs <- mean((z - f_hat_s_hs)^2)
mse_s_ridge <- mean((z - f_hat_s_n)^2)
mse_f_hs <- mean((z - f_hat_f_hs)^2)
mse_f_ridge <- mean((z - f_hat_f_n)^2)
mse_l_hs <- mean((z - f_hat_l_hs)^2)
mse_l_ridge <- mean((z - f_hat_l_n)^2)

cat(sprintf(
  "MSE - Bs / Horseshoe: %.5f\nMSE - Bs / Ridge: %.5f\nMSE - F / Horseshoe: %.5f\nMSE - F / Ridge: %.5f\nMSE - Leg / Horseshoe: %.5f\nMSE - Leg / Ridge: %.5f\n",
  mse_s_hs, mse_s_ridge, mse_f_hs, mse_f_ridge, mse_l_hs, mse_l_ridge
))
