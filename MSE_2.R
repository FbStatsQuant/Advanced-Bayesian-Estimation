library(splines)
library(MLmetrics)
library(Hmisc)
library(horseshoe)
library(ggplot2)
set.seed(2024)

# Read data
df <- read.csv("CSIRO_Recons_gmsl_yr_2019.csv", header = TRUE, sep = ",")
x <- unlist(df["Time"])
y <- unlist(df["GMSL"])
n <- length(y)

# Range of J values to test
J_seq <- seq(20, 180, by = 20)
rmse_splines <- numeric(length(J_seq))
mse_splines <- numeric(length(J_seq))
mae_splines <- numeric(length(J_seq))

# Test each J value
for(i in seq_along(J_seq)) {
  J <- J_seq[i]

  # B-Splines with Horseshoe
  t <- c(min(x) - 0.02, min(x) - 0.01, seq(min(x), max(x), length.out = J-4), max(x) + 0.01, max(x) + 0.02)
  Sp_order <- 3
  B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = FALSE)
  Hs_s <- horseshoe(y, B_s, method.tau = "halfCauchy", method.sigma = "Jeffreys",
                    burn = 500, nmc = 2000, thin = 1)
  beta_s_hat <- as.vector(unlist(Hs_s[1]))
  f_hat_s <- B_s %*% beta_s_hat

  # Calculate error metrics
  mse_splines[i] <- mean((y - f_hat_s)^2)
  rmse_splines[i] <- sqrt(mse_splines[i])
  mae_splines[i] <- mean(abs(y - f_hat_s))

  # Print progress
  cat("Completed J =", J, "(", i, "of", length(J_seq), ")\n")
  cat("  MSE:", mse_splines[i], "\n")
  cat("  RMSE:", rmse_splines[i], "\n")
  cat("  MAE:", mae_splines[i], "\n\n")
}

# Create results dataframe
results_df <- data.frame(
  J = J_seq,
  MSE = mse_splines,
  RMSE = rmse_splines,
  MAE = mae_splines
)

# Create plot of RMSE vs J
ggplot(results_df, aes(x = J, y = RMSE)) +
  geom_point(color = "red", size = 3) +
  geom_line(color = "red") +
  labs(
    title = "RMSE for B-Splines with Horseshoe Prior (GMSL)",
    x = "Number of Basis Functions (J)",
    y = "Root Mean Squared Error"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Find optimal J
optimal_idx <- which.min(rmse_splines)
optimal_J <- J_seq[optimal_idx]

cat("Optimal results for B-Splines:\n")
cat("J =", optimal_J, "\n")
cat("MSE =", mse_splines[optimal_idx], "\n")
cat("RMSE =", rmse_splines[optimal_idx], "\n")
cat("MAE =", mae_splines[optimal_idx], "\n")

# Generate optimal fit
t_opt <- c(min(x) - 0.02, min(x) - 0.01, seq(min(x), max(x), length.out = optimal_J-4), max(x) + 0.01, max(x) + 0.02)
B_s_opt <- splineDesign(t_opt, x, ord = 3, outer.ok = TRUE, sparse = FALSE)
Hs_s_opt <- horseshoe(y, B_s_opt, method.tau = "halfCauchy", method.sigma = "Jeffreys",
                      burn = 1000, nmc = 5000, thin = 1)
beta_s_opt <- as.vector(unlist(Hs_s_opt[1]))
f_hat_s_opt <- B_s_opt %*% beta_s_opt

# Create final plot with optimal fit
df_plot <- data.frame(
  x = x,
  y = y,
  fit = f_hat_s_opt
)

ggplot(df_plot, aes(x = x)) +
  geom_point(aes(y = y, color = "Observed"), size = 1, alpha = 0.7) +
  geom_line(aes(y = fit, color = "B-Spline Estimate"), linewidth = 0.8) +
  scale_color_manual(values = c("B-Spline Estimate" = "red", "Observed" = "black")) +
  labs(
    title = "CSIRO Global Mean Sea Level Estimation",
    subtitle = paste("B-Splines with J =", optimal_J, "| RMSE =", round(rmse_splines[optimal_idx], 4)),
    x = "Years",
    y = "Millimeters",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save results if needed
# write.csv(results_df, "gmsl_bsplines_results.csv", row.names = FALSE)
