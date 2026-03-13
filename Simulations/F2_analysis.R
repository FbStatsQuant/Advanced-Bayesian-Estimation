# ============================================================
#  F1_analysis.R
#  Load chains from F1_chains_trial.rds and compute:
#    - True f(x), fitted f(x) by mean and median of betas
#    - RMSE (mean and median estimator) -- matches empirical L2 norm in proof
#    - a_n: number of |beta_j| > 1/J_n
#    - Regression: log(RMSE) ~ log(n), log(a_n) ~ log(n)
#  Exports summary table to Excel
# ============================================================

rds_path <- "C:/Users/felip/OneDrive - Rice University/Overcomplete random series priors/Advanced-Bayesian-Estimation/Simulations/F2_chains_trial.rds"

results <- readRDS(rds_path)

# True function
f <- function(x) 1.5*(abs(x-0.1))^3 - 5*(abs(x-0.4))^3

fourier_basis <- function(x, K) {
  basis <- matrix(1, nrow = length(x), ncol = K + 1)
  for (k in 1:K) basis[, k + 1] <- cos(pi * k * x)
  basis
}

# ============================================================
#  Storage
# ============================================================
n_seq     <- sapply(results, function(r) r$n)
J_seq     <- sapply(results, function(r) r$J)
RMSE_mean <- numeric(length(results))
RMSE_med  <- numeric(length(results))
a_n       <- numeric(length(results))

# ============================================================
#  Main loop
# ============================================================
for (i in seq_along(results)) {
  
  r    <- results[[i]]
  n    <- r$n
  J    <- r$J
  x    <- r$x
  y    <- r$y
  samp <- r$BetaSamples   # J x 4000
  
  B    <- fourier_basis(x, J)
  
  beta_mean <- rowMeans(samp)
  beta_med  <- apply(samp, 1, median)
  
  f_hat_mean <- B %*% beta_mean
  f_hat_med  <- B %*% beta_med
  
  f_true <- f(x)
  
  # RMSE: matches empirical L2 norm used in the proof
  RMSE_mean[i] <- sqrt(mean((f_hat_mean - f_true)^2))
  RMSE_med[i]  <- sqrt(mean((f_hat_med  - f_true)^2))
  
  # a_n: number of |beta_j| > 1/J_n (posterior median)
  a_n[i] <- sum(abs(beta_med) > 1 / n)
  
  cat(sprintf("n = %4d | J = %3d | a_n = %3d (%.1f%% of J) | RMSE_mean = %.6f | RMSE_med = %.6f\n",
              n, J, a_n[i], 100 * a_n[i] / J, RMSE_mean[i], RMSE_med[i]))
}


val <- 10
# ============================================================
#  Regressions
# ============================================================
cat("\n--- Regression: log(RMSE_mean) ~ log(n) ---\n")
fit_rmse_mean <- lm(log(RMSE_mean[-(1:val)]) ~ log(n_seq[-(1:val)]))
print(summary(fit_rmse_mean))
cat(sprintf("Estimated convergence rate (mean): %.4f\n", coef(fit_rmse_mean)[2]))

cat("\n--- Regression: log(RMSE_med) ~ log(n) ---\n")
fit_rmse_med <- lm(log(RMSE_med[-(1:val)]) ~ log(n_seq[-(1:val)]))
print(summary(fit_rmse_med))
cat(sprintf("Estimated convergence rate (median): %.4f\n", coef(fit_rmse_med)[2]))

cat("\n--- Regression: log(a_n) ~ log(n) ---\n")
fit_an <- lm(log(a_n[-(1:val)]) ~ log(n_seq[-(1:val)]))
print(summary(fit_an))
cat(sprintf("Estimated sparsity order: %.4f\n", coef(fit_an)[2]))

# ============================================================
#  Export summary table to Excel
# ============================================================
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
library(writexl)

summary_table <- data.frame(
  n           = n_seq,
  J           = J_seq,
  RMSE_mean   = RMSE_mean,
  RMSE_median = RMSE_med,
  a_n         = a_n,
  pct_nonzero = round(100 * a_n / J_seq, 2)
)

excel_path <- "C:/Users/felip/OneDrive - Rice University/Overcomplete random series priors/Advanced-Bayesian-Estimation/Simulations/F2_summary.xlsx"
write_xlsx(summary_table, excel_path)
cat(sprintf("\nSummary table saved to: %s\n", excel_path))

# ============================================================
#  Plot 1: log(n) vs log(RMSE_median)
# ============================================================
plot(log(n_seq), log(RMSE_med),
     xlab = "log(n)", ylab = "log(RMSE)",
     main = "log(n) vs log(RMSE) -- posterior median",
     pch = 16, col = "darkgreen")
abline(fit_rmse_med, col = "darkgreen", lwd = 2)
legend("topright",
       legend = sprintf("slope = %.4f", coef(fit_rmse_med)[2]),
       col = "darkgreen", lwd = 2, bty = "n")

# ============================================================
#  Plot 2: log(n) vs log(RMSE_mean)
# ============================================================
plot(log(n_seq), log(RMSE_mean),
     xlab = "log(n)", ylab = "log(RMSE)",
     main = "log(n) vs log(RMSE) -- posterior mean",
     pch = 16, col = "steelblue")
abline(fit_rmse_mean, col = "steelblue", lwd = 2)
legend("topright",
       legend = sprintf("slope = %.4f", coef(fit_rmse_mean)[2]),
       col = "steelblue", lwd = 2, bty = "n")

# ============================================================
#  Plot 3: log(n) vs log(a_n)
# ============================================================
plot(log(n_seq), log(a_n),
     xlab = "log(n)", ylab = "log(a_n)",
     main = "log(n) vs log(a_n) -- sparsity",
     pch = 16, col = "firebrick")
abline(fit_an, col = "firebrick", lwd = 2)
legend("topleft",
       legend = sprintf("slope = %.4f", coef(fit_an)[2]),
       col = "firebrick", lwd = 2, bty = "n")



### Plots in ggplot

data <- data.frame(
  n = n_seq,
  RMSE_mean = RMSE_mean,
  RMSE_med = RMSE_med,
  a_n = a_n
)

library(ggplot2)

#Mean
ggplot(data, aes(x = log(n), y = log(RMSE_mean))) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  labs(title = "log(n) vs log(RMSE) -- posterior mean",
       x = "log(n)", y = "log(RMSE)") +
  annotate("text", x = min(log(n_seq)), y = max(log(RMSE_mean)),
           label = sprintf("slope = %.4f", coef(fit_rmse_mean)[2]),
           hjust = 0, vjust = 1, color = "steelblue")

#Median
ggplot(data, aes(x = log(n), y = log(RMSE_med))) +  
  geom_point(color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  labs(title = "log(n) vs log(RMSE) -- posterior median",
       x = "log(n)", y = "log(RMSE)") +
  annotate("text", x = min(log(n_seq)), y = max(log(RMSE_med)),
           label = sprintf("slope = %.4f", coef(fit_rmse_med)[2]),
           hjust = 0, vjust = 1, color = "darkgreen")

#Sparsity
ggplot(data, aes(x = log(n), y = log(a_n))) +
  geom_point(color = "firebrick") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(title = "log(n) vs log(a_n) -- sparsity",
       x = "log(n)", y = "log(a_n)") +
  annotate("text", x = min(log(n_seq)), y = max(log(a_n)),
           label = sprintf("slope = %.4f", coef(fit_an)[2]),
           hjust = 0, vjust = 1, color = "firebrick")

