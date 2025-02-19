# Load necessary libraries
library(splines)
library(sparsehorseshoe)
library(MLmetrics)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Parameters for the AR(1) process
n <- 1000  # Length of the time series
phi <- 0.7  # Linear AR(1) coefficient
stationary_ar1 <- numeric(n)
stationary_ar1[1] <- rnorm(1)

# Generate the linear AR(1) process
for (i in 2:n) {
  stationary_ar1[i] <- phi * stationary_ar1[i-1] + rnorm(1)
}

# Plot the linear AR(1) process
plot.ts(stationary_ar1, main="Linear AR(1) Process", ylab="Value", xlab="Time")

# Prepare data for LP (Sparse Horseshoe)
y <- stationary_ar1[2:n]
x <- stationary_ar1[1:(n-1)]

# Set up a grid for J values (e.g., 400 to 1600, step of 200)
J_values <- seq(400, 1600, by=200)

# To store MSE results for each J
mse_results <- numeric(length(J_values))

# Loop over each J to calculate MSE
for (i in 1:length(J_values)) {
  J <- J_values[i]
  
  # Generate B-spline basis functions for each J
  t <- c(0.998 * min(x), 0.999 * min(x), seq(min(x), max(x), length.out = J - 4), max(x) * 1.001, max(x) * 1.002)
  Sp_order <- 3
  B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
  
  # Apply sparse horseshoe method (LP) for Bayesian regression
  Hs4 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), s = 0.5, method.sigma = c("Jeffreys"), 
                     burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
  
  # Get the beta estimates
  beta_hat4 <- as.vector(unlist(Hs4[1]))
  
  # Compute MSE
  mse_results[i] <- sqrt(MSE(y, as.numeric(B %*% beta_hat4)))
}

# Plot J vs sqrt(MSE) using ggplot2
mse_df <- data.frame(J = J_values, sqrt_MSE = mse_results)

ggplot(mse_df, aes(x=J, y=sqrt_MSE)) +
  geom_line(color="blue", size=1) +
  geom_point(color="blue", size=3) +
  labs(title="J vs sqrt(MSE) for Linear AR(1) Process", x="J (Number of B-spline basis functions)", y="sqrt(MSE)") +
  theme_minimal()

# Choose a specific J for fitting the model and plotting predictions
chosen_J <- 800  # You can choose any J value you prefer

# Generate B-spline basis functions for chosen J
t <- c(0.998 * min(x), 0.999 * min(x), seq(min(x), max(x), length.out = chosen_J - 4), max(x) * 1.001, max(x) * 1.002)
Sp_order <- 3
B_chosen <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)

# Apply sparse horseshoe method for chosen J
Hs_chosen <- horseshoesp(y, B_chosen, method.tau = c("truncatedCauchy"), s = 0.5, method.sigma = c("Jeffreys"), 
                         burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)

# Get the beta estimates for chosen J
beta_hat_chosen <- as.vector(unlist(Hs_chosen[1]))

# Compute the fitted values B * beta_hat
fitted_values <- as.numeric(B_chosen %*% beta_hat_chosen)

# Fit a Gaussian Process for comparison
library(kernlab)

# Use Gaussian Process with RBF kernel
gp_model <- gausspr(x = x, y = stationary_ar1[2:n], kernel = "rbfdot", kpar = list(sigma = 0.1))
gp_predictions <- predict(gp_model, x)

# Create a data frame for plotting the AR(1) process, LP Fitted, and GP Predictions
plot_df <- data.frame(x = x, AR1 = stationary_ar1[2:n], LP_Fitted = fitted_values, GP_Predictions = gp_predictions)

# Plot the AR(1) process, LP fitted values, and GP predictions using ggplot2
ggplot(plot_df, aes(x=x)) +
  geom_point(aes(y=AR1), color="lightblue", size=2) +   # AR(1) process as points in light blue
  geom_smooth(aes(y=LP_Fitted), color="red", size=1.2, linetype="dashed", method="loess") +  # LP fitted, smoothed with loess
  geom_smooth(aes(y=GP_Predictions), color="blue", size=1.2, linetype="dotted", method="loess") + # GP predictions, smoothed with loess
  labs(title="Linear AR(1) Process with LP and GP Predictions", x="x_(t-1)", y="x_t") +
  theme_minimal() +
  theme(legend.position="top")
