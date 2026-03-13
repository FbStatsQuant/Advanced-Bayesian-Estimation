library(readxl)
library(horseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)
library(ggplot2)
set.seed(2024)

# Read data
df <- read.csv("SP500.csv", header = TRUE, sep = ",")
y <- rev(as.vector(unlist(df['Close'])))
n <- length(y)
x_orig <- 1:n                      # original time index
x_scaled <- x_orig / n             # rescaled to (0,1)

# Fourier basis function
fourier_basis <- function(x, K) {
  basis <- matrix(1, nrow = length(x), ncol = 2 * K + 1)
  for (k in 1:K) {
    basis[, 2 * k]     <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }
  basis
}

# Range of J values to test
J_seq <- seq(50, 500, by = 50)
mse_values <- numeric(length(J_seq))
rmse_values <- numeric(length(J_seq))
mae_values <- numeric(length(J_seq))
mape_values <- numeric(length(J_seq))

# Test each J value
for(i in seq_along(J_seq)) {
  J <- J_seq[i]

  # Create Fourier basis
  B_f <- fourier_basis(x_scaled, floor(J / 2))

  # Horseshoe estimation
  Hs_f <- horseshoe(y, B_f, method.tau = "halfCauchy", method.sigma = "Jeffreys")
  theta_f_hs <- unlist(Hs_f$BetaHat)
  f_hat_f_hs <- B_f %*% theta_f_hs

  # Calculate error metrics
  mse_values[i] <- mean((y - f_hat_f_hs)^2)
  rmse_values[i] <- sqrt(mse_values[i])
  mae_values[i] <- mean(abs(y - f_hat_f_hs))
  mape_values[i] <- mean(abs((y - f_hat_f_hs)/y)) * 100

  # Print progress
  cat("Completed J =", J, "(", i, "of", length(J_seq), ")\n")
}

# Create results dataframe
results_df <- data.frame(
  J = J_seq,
  MSE = mse_values,
  RMSE = rmse_values,
  MAE = mae_values,
  MAPE = mape_values
)

# Create plot of MSE vs J
ggplot(results_df, aes(x = J, y = RMSE)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue") +
  labs(
    title = "RMSE vs Number of Basis Functions (S&P 500)",
    x = "Number of Basis Functions (J)",
    y = "Mean Squared Error"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )






# Find optimal J
optimal_idx <- which.min(mse_values)
optimal_J <- J_seq[optimal_idx]
cat("Optimal number of basis functions (J):", optimal_J,
    "\nMSE:", mse_values[optimal_idx],
    "\nRMSE:", rmse_values[optimal_idx],
    "\nMAE:", mae_values[optimal_idx],
    "\nMAPE:", mape_values[optimal_idx], "%\n")

# Generate plot with optimal J
J <- optimal_J
B_f <- fourier_basis(x_scaled, floor(J / 2))
Hs_f <- horseshoe(y, B_f, method.tau = "halfCauchy", method.sigma = "Jeffreys")
theta_f_hs <- unlist(Hs_f$BetaHat)
f_hat_f_hs <- B_f %*% theta_f_hs

# Exclude first and last 10 points
cut <- 10
data <- data.frame(
  Time = x_orig[(cut+1):(n - cut)],
  True = y[(cut+1):(n - cut)],
  Estimated = as.vector(f_hat_f_hs)[(cut+1):(n - cut)]
)

# Subset every 10th point for estimated values
scatter_data <- data[seq(1, nrow(data), by = 10), ]

# Create final plot with optimal J
ggplot() +
  geom_line(data = data, aes(x = Time, y = True, color = "True (Observed)"), size = 0.8) +
  geom_point(data = scatter_data, aes(x = Time, y = Estimated, color = "Estimated (Horseshoe + Fourier)"),
             size = 1.1, alpha = 0.7) +
  scale_color_manual(values = c("True (Observed)" = "lightblue",
                                "Estimated (Horseshoe + Fourier)" = "blue")) +
  labs(
    title = paste("True vs Estimated S&P 500 (Horseshoe + Fourier, J =", optimal_J, ")"),
    x = "Time Index",
    y = "Closing Price",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(hjust = 0.5)
  )

# Save results if needed
# write.csv(results_df, "sp500_fourier_horseshoe_results.csv", row.names = FALSE)
