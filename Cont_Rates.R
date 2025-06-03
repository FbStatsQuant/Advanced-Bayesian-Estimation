library(horseshoe)
library(splines)
library(ggplot2)

set.seed(2024)

# Settings
f <- function(x) 1.5*(abs(x-0.1))^3 - 5*(abs(x-0.4))^3
sigma <- 0.08
b <- 0.5
Sp_order <- 3
n_seq <- seq(1000, 10000, by = 1000)

# Initialize storage
mse_bs <- numeric(length(n_seq))
mse_fourier <- numeric(length(n_seq))

# Main loop
for (i in seq_along(n_seq)) {
  n <- n_seq[i]
  J <- floor(n^b)
  x <- seq(0, 1, length.out = n)
  z <- f(x)
  y <- z + rnorm(n, 0, sigma)

  # B-spline basis
  t <- seq(0.0001, 0.9999, length.out = J)
  B_s <- bSpline(x, knots = t, degree = Sp_order, intercept = FALSE)
  Hs_s <- horseshoe(y, B_s, method.tau = "halfCauchy", method.sigma = "fixed", Sigma2 = sigma^2)
  f_hat_s <- B_s %*% Hs_s$BetaHat
  mse_bs[i] <- mean((z - f_hat_s)^2)

  # Fourier basis
  fourier_basis <- function(x, K) {
    mat <- matrix(1, nrow = length(x), ncol = 2 * K + 1)
    for (k in 1:K) {
      mat[, 2 * k]     <- sin(2 * pi * k * x)
      mat[, 2 * k + 1] <- cos(2 * pi * k * x)
    }
    mat
  }

  B_f <- fourier_basis(x, floor(J / 2))
  Hs_f <- horseshoe(y, B_f, method.tau = "halfCauchy", method.sigma = "fixed", Sigma2 = sigma^2)
  f_hat_f <- B_f %*% Hs_f$BetaHat
  mse_fourier[i] <- mean((z - f_hat_f)^2)
}

# Prepare data
df <- data.frame(
  n = n_seq,
  logn = log(n_seq),
  log_mse_bs = log(mse_bs),
  log_mse_fourier = log(mse_fourier)
)

# Fit models
model_bs <- lm(log_mse_bs ~ logn, data = df)
model_fourier <- lm(log_mse_fourier ~ logn, data = df)

summary(model_bs)
summary(model_fourier)

cat("B-splines regression slope:", coef(model_bs)[2], "\n")
cat("Fourier regression slope:", coef(model_fourier)[2], "\n")
slope_bs <- round(coef(model_bs)[2], 3)
slope_fourier <- round(coef(model_fourier)[2], 3)

# Plot: B-splines
ggplot(df, aes(x = logn, y = log_mse_bs)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = paste0("log(MSE) vs log(n) — Horseshoe with B-splines (slope = ", slope_bs, ")"),
    x = "log(n)", y = "log(MSE)"
  ) +
  theme_minimal()

# Plot: Fourier
ggplot(df, aes(x = logn, y = log_mse_fourier)) +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = paste0("log(MSE) vs log(n) — Horseshoe with Fourier (slope = ", slope_fourier, ")"),
    x = "log(n)", y = "log(MSE)"
  ) +
  theme_minimal()
