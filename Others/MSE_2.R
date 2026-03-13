library(splines2)
library(MLmetrics)
library(Hmisc)
library(horseshoe)
library(ggplot2)
set.seed(2024)

n <- 1000
x <- seq(0,1, length.out = n)
f <- function(x) 1.5*(abs(x-0.1))^3 - 5*(abs(x-0.4))^3
z <- f(x)
sd <- 0.07
y <- z + rnorm(n, mean = 0, sd = sd)


# Fourier Matrix
fourier_basis <- function(x, K) {

  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)

  for (k in 1:K) {
    basis[, 2 * k]   <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }

  return(basis)
}


# Range of J values to test
J_seq <- seq(300, n, by = 50)
rmse_fourier <- numeric(length(J_seq))
mse_fourier <- numeric(length(J_seq))
mae_fourier <- numeric(length(J_seq))

# Test each J value
for(i in seq_along(J_seq)) {
  J <- J_seq[i]

  # B-fourier with Horseshoe
  t <- c(min(x) - 0.02, min(x) - 0.01, seq(min(x), max(x), length.out = J-4), max(x) + 0.01, max(x) + 0.02)
  Sp_order <- 3
  B_f <- fourier_basis(x, floor(J/2))
  Hs_f <- horseshoe(y, B_f, method.tau = "halfCauchy", method.sigma = "fixed", Sigma2 = sd^2,
                    burn = 1000, nmc = 3000, thin = 1)
  beta_f_hat <- as.vector(unlist(Hs_f[1]))
  f_hat_f <- B_f %*% beta_f_hat

  # Calculate error metrics
  mse_fourier[i] <- mean((z - f_hat_f)^2)
  rmse_fourier[i] <- sqrt(mse_fourier[i])
  mae_fourier[i] <- mean(abs(z - f_hat_f))

  # Print progress
  cat("Completed J =", J, "(", i, "of", length(J_seq), ")\n")
  cat("  MSE:", mse_fourier[i], "\n")
  cat("  RMSE:", rmse_fourier[i], "\n")
  cat("  MAE:", mae_fourier[i], "\n\n")
}

# Create results dataframe
results_df <- data.frame(
  J = J_seq,
  MSE = mse_fourier,
  RMSE = rmse_fourier,
  MAE = mae_fourier
)

# Create plot of RMSE vs J
ggplot(results_df, aes(x = J, y = MSE)) +
  geom_point(color = "red", size = 3) +
  geom_line(color = "red") +
  labs(
    title = "MSE for Fourier with Horseshoe Prior",
    x = "Number of Basis Functions (J)",
    y = "Mean Squared Error"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  ylim(c(0.0005,0.0007))  # This line fixes the y-axis


