library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)
library(splines2)

set.seed(1991)


n <- 5000
b <- 0.6
J <- floor(n^b)
x <- seq(0,1, length.out = n)

f <- function(x) 0.5*(abs(x-0.2))^5 + 2*(abs(x-0.5))^5 -0.5*(abs(x-0.9))^5
z <- f(x)
sd <- 0.1
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

t <- seq(0.0001, 0.9999, length.out = J)
Sp_order <- 3
B_s <- bSpline(x, knots = t, degree = Sp_order, intercept = FALSE)

#Hs (halfCauchy tau)

Hs_s <- horseshoe(y, B_s, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                  Sigma2 = sd^2, nmc = 5000, burn = 1000)

theta_s_hs <- drop(Hs_s$BetaHat)
f_hat_s_hs <- drop(B_s%*%theta_s_hs)

#Hs (fixed tau = 1/J)

Hs_s_fixed <- horseshoe(y, B_s, method.tau = c("fixed"), tau = 1/J, method.sigma = c("fixed"),
                        Sigma2 = sd^2, nmc = 5000, burn = 1000)

theta_s_hs_fixed <- drop(Hs_s_fixed$BetaHat)
f_hat_s_hs_fixed <- drop(B_s%*%theta_s_hs_fixed)

#Ridge

lambda_s <- cv.glmnet(B_s,y, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_s/((sd^2)) * diag(ncol(B_s)) + (1/(sd^2)) * crossprod(B_s))
theta_s_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_s, y))
f_hat_s_n <- B_s%*%theta_s_n

## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_s_hs, hs_fixed = f_hat_s_hs_fixed, ridge = f_hat_s_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe (halfCauchy)"), size = 1, alpha = 1) +
  geom_point(aes(y = hs_fixed, color = "Horseshoe (fixed tau)"), size = 1, alpha = 1) +
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1, alpha = 1) +
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (B-Splines)",
    x = "x",
    y = "Value",
  ) +
  scale_color_manual(values = c("Horseshoe (halfCauchy)" = "blue",
                                "Horseshoe (fixed tau)" = "darkgreen",
                                "Ridge Regression" = "lightblue",
                                "True Function" = "black")) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.11),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
  )

# Fourier (full sin+cos basis)

fourier_basis <- function(x, K) {

  basis <- matrix(1, nrow = length(x), ncol = 2 * K + 1)

  for (k in 1:K) {
    basis[, 2 * k]     <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }

  return(basis)
}

B_f <- fourier_basis(x, num_fourier)

## Hs (fixed tau = 1/J)

Hs_f <- horseshoe(y, B_f, method.tau = c("fixed"), tau = 1/(J), method.sigma = c("fixed"),
                  Sigma2 = sd^2, nmc = 20000, burn = 4000, thin = 10)
theta_f_hs <- drop(Hs_f$BetaHat)
f_hat_f_hs <- drop(B_f%*%theta_f_hs)

## Hs (halfCauchy tau)

Hs_f_hc <- horseshoe(y, B_f, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                     Sigma2 = sd^2, nmc = 20000, burn = 4000)
theta_f_hs_hc <- drop(Hs_f_hc$BetaHat)
f_hat_f_hs_hc <- drop(B_f%*%theta_f_hs_hc)

## Ridge

lambda_f <- cv.glmnet(B_f,y, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (lambda_f/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
theta_f_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
f_hat_f_n <- B_f%*%theta_f_n

## Plot

subset_indices <- seq(1, n, by = 2)
data <- data.frame(x = x, true = z, hs = f_hat_f_hs, hs_hc = f_hat_f_hs_hc, ridge = f_hat_f_n)[subset_indices, ]

ggplot(data, aes(x = x)) +
  geom_point(aes(y = hs, color = "Horseshoe (fixed tau)"), size = 1.3, alpha = 1) +
  geom_point(aes(y = hs_hc, color = "Horseshoe (halfCauchy)"), size = 1.3, alpha = 1) +
  geom_point(aes(y = ridge, color = "Ridge Regression"), size = 1.7, alpha = 1) +
  geom_line(aes(y = true, color = "True Function"), linewidth = 0.8) +
  labs(
    title = "True Function, Horseshoe, and Ridge Estimates (Fourier)",
    x = "x",
    y = "Value",
  ) +
  scale_color_manual(values = c("Horseshoe (fixed tau)" = "blue",
                                "Horseshoe (halfCauchy)" = "red",
                                "Ridge Regression" = "lightblue",
                                "True Function" = "black")) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.11),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_blank(),
  )



(mse_s_hs <- mean((z - f_hat_s_hs)^2))
(mse_s_hs_fixed <- mean((z - f_hat_s_hs_fixed)^2))
(mse_s_ridge <- mean((z - f_hat_s_n)^2))
(mse_f_hs <- mean((z - f_hat_f_hs)^2))
(mse_f_hs_hc <- mean((z - f_hat_f_hs_hc)^2))
(mse_f_ridge <- mean((z - f_hat_f_n)^2))

# Histogram of absolute horseshoe coefficients

threshold <- 1 / J

abs_coef_f <- data.frame(abs_coef = abs(theta_f_hs),
                         index    = seq_along(theta_f_hs))
pct_f <- mean(abs_coef_f$abs_coef < threshold) * 100

ggplot(abs_coef_f, aes(x = abs_coef)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("%.1f%% < 1/J", pct_f), size = 4) +
  labs(
    title = "Histogram of |Horseshoe Coefficients| (Fourier, fixed tau)",
    x = "|Coefficient|",
    y = "Count"
  ) +
  theme_minimal()

abs_coef_f_hc <- data.frame(abs_coef = abs(theta_f_hs_hc),
                             index    = seq_along(theta_f_hs_hc))
pct_f_hc <- mean(abs_coef_f_hc$abs_coef < threshold) * 100

ggplot(abs_coef_f_hc, aes(x = abs_coef)) +
  geom_histogram(bins = 100, fill = "red", color = "white") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("%.1f%% < 1/J", pct_f_hc), size = 4) +
  labs(
    title = "Histogram of |Horseshoe Coefficients| (Fourier, halfCauchy tau)",
    x = "|Coefficient|",
    y = "Count"
  ) +
  theme_minimal()

abs_coef_s_hs <- data.frame(abs_coef = abs(theta_s_hs),
                             index    = seq_along(theta_s_hs))
pct_s_hs <- mean(abs_coef_s_hs$abs_coef < threshold) * 100

ggplot(abs_coef_s_hs, aes(x = abs_coef)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("%.1f%% < 1/J", pct_s_hs), size = 4) +
  labs(
    title = "Histogram of |Horseshoe Coefficients| (B-Splines, halfCauchy tau)",
    x = "|Coefficient|",
    y = "Count"
  ) +
  theme_minimal()

abs_coef_s_fixed <- data.frame(abs_coef = abs(theta_s_hs_fixed),
                                index    = seq_along(theta_s_hs_fixed))
pct_s_fixed <- mean(abs_coef_s_fixed$abs_coef < threshold) * 100

ggplot(abs_coef_s_fixed, aes(x = abs_coef)) +
  geom_histogram(bins = 100, fill = "darkgreen", color = "white") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("%.1f%% < 1/J", pct_s_fixed), size = 4) +
  labs(
    title = "Histogram of |Horseshoe Coefficients| (B-Splines, fixed tau)",
    x = "|Coefficient|",
    y = "Count"
  ) +
  theme_minimal()
