library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)

set.seed(1991)


n <- 1000
b <- 0.6
J <- floor(n^b)
x <- seq(0,1, length.out = n)

f <- function(x) abs(x-0.5)
z <- f(x)
sd <- 0.15
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

t <- seq(0.0001,0.9999, length.out = J)
Sp_order <- 3
B_s <- bSpline(x, knots = t, degree = Sp_order, intercept = FALSE)

#Hs

Hs_s <- horseshoe(y, B_s, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                  Sigma2 = sd^2)

theta_s_hs <- unlist(Hs_s$BetaHat)
f_hat_s_hs <- B_s%*%theta_s_hs

#Ridge

lambda_s <- cv.glmnet(B_s,y, alpha=0, nfolds = 20)$lambda.min
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
  ) +
  scale_color_manual(values = c("Horseshoe Prior" = "blue",
                                "Ridge Regression" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.11),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    #plot.caption = element_text(hjust = 0.1, vjust = 1, size = 10, face = "italic")  # Centered caption below the plot
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

Hs_f <- horseshoe(y, B_f, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                  Sigma2 = sd^2)
theta_f_hs <- unlist(Hs_f$BetaHat)
f_hat_f_hs <- B_f%*%theta_f_hs

## Ridge

lambda_f <- cv.glmnet(B_f,y, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_f/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
theta_f_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
f_hat_f_n <- B_f%*%theta_f_n

## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_f_hs, ridge = f_hat_f_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe Prior"), size = 1.3, alpha = 1) +
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1.7, alpha = 1) +
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (Fourier)",
    x = "x",
    y = "Value",
  ) +
  scale_color_manual(values = c("Horseshoe Prior" = "blue",
                                "Ridge Regression" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.11),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
    #plot.caption = element_text(hjust = 0.1, vjust = 1, size = 10, face = "italic")  # Centered caption below the plot
  )



(mse_s_hs <- mean((z - f_hat_s_hs)^2))
(mse_s_ridge <- mean((z - f_hat_s_n)^2))
(mse_f_hs <- mean((z - f_hat_f_hs)^2))
(mse_f_ridge <- mean((z - f_hat_f_n)^2))
