library(splines)
library(MLmetrics)
library(Hmisc)
library(forecast)
library(ggplot2)

set.seed(123)

n <- 2000
b <- 1.2
J <- floor(n^b)


phi <- 0.7
stationary_ar1 <- numeric(n)
stationary_ar1[1] <- rnorm(1)

# Generate AR(1) process
for (i in 2:n) {
  stationary_ar1[i] <- phi * stationary_ar1[i-1] + rnorm(1)
}

# Plot the AR(1) series
plot.ts(stationary_ar1, main="Weakly Stationary AR(1) Time Series", ylab="Value", xlab="Time")

# #ARIMA
# data <- ts(y)
# best_model <- auto.arima(data)
# summary(best_model)


#With respect to time

# x <- c(1:length(stationary_ar1))
# y <- stationary_ar1
# J <- 500
# t <- c(0.998*min(x), 0.999*min(x), seq(min(x),max(x), length.out = J-4),max(x)*1.001, max(x)*1.002)
# Sp_order <- 3
# B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
#                   sparse = TRUE)
# Hs4 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), k = max(0.5,b/2),
#                    method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
# beta_hat4 <- as.vector(unlist(Hs4[1]))
# lp_prediction <- as.numeric(B%*%beta_hat4)
# 
# subset_indices <- seq(1, length(stationary_ar1), by=4)
# 
# plot_df <- data.frame(
#   Time = 1:length(stationary_ar1),
#   Actual = stationary_ar1,
#   LP_Estimation = lp_prediction
# )
# 
# plot_df_subset <- plot_df[subset_indices, ]
# 
# ggplot() +
#   geom_point(data = plot_df_subset, aes(x = Time, y = Actual, color = "Actual values"), size = 2, alpha = 0.8) +
#   geom_line(data = plot_df, aes(x = Time, y = LP_Estimation, color = "Estimation"), size = 0.8, alpha = 0.7) +
#   labs(x = "Observations",
#        y = "Value") +
#   scale_color_manual(values = c("Actual values" = "blue", "Estimation" = "tomato")) +
#   theme_minimal() +
#   theme(legend.title = element_blank())


y <- stationary_ar1[2:n]
x <- stationary_ar1[1:n-1]
t <- c(0.998*min(x), 0.999*min(x), seq(min(x),max(x), length.out = J-4),max(x)*1.001, max(x)*1.002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)
Hs4 <- horseshoesp(y, B, method.tau = c("fixed"), k=max(0.5,b/3), tau = b/2,   
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat4 <- as.vector(unlist(Hs4[1]))
lp_prediction <- as.numeric(B%*%beta_hat4)

s = sqrt(J)
# 
# gp_model <- gausspr(x = x, y = y, kernel = "rbfdot", kpar = list(sigma = 0.1))
# gp_predictions <- predict(gp_model, x)
# sqrt(MSE(y,gp_predictions))
# 
# rmse_results <- numeric()
# 
# for (J in seq(400, 1800, by=200)) {
#   
#   # Generate B-spline basis functions for each J
#   t <- c(0.998 * min(x), 0.999 * min(x), seq(min(x), max(x), length.out = J - 4), max(x) * 1.001, max(x) * 1.002)
#   Sp_order <- 3
#   B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
#   
#   # Apply sparse horseshoe method (LP) for Bayesian regression
#   Hs4 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), s = 0.5, method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
#   
#   # Get the beta estimates
#   beta_hat4 <- as.vector(unlist(Hs4[1]))
#   
#   # Compute predictions using the fitted model
#   lp_prediction <- as.numeric(B %*% beta_hat4)
#   
#   # Compute RMSE for LP predictions
#   rmse <- sqrt(MSE(y, lp_prediction))
#   
#   # Store the RMSE for the current J value
#   rmse_results <- c(rmse_results, rmse)
# }
# 
# # Plot J vs RMSE using ggplot2
# rmse_df <- data.frame(J = seq(400, 1800, by=200), RMSE = rmse_results)
# 
# ggplot(rmse_df, aes(x=J, y=RMSE)) +
#   geom_line(color="blue", size=1) +
#   geom_point(color="blue", size=3) +
#   labs(x="J (Number of B-spline basis functions)", 
#        y="RMSE") +
#   theme_minimal() +
#   ylim(0.89, 0.94)  # Set y-axis limits from 0.8 to 1.0
# 
# 
# 
# # Fit the underlying AR(1) process (straight line with slope = 0.7)
# ar_line <- 0.7 * x  # Underlying AR(1) process line with slope 0.7
# 
# # Subset the data to show every 10th point
# subset_indices <- seq(1, length(x), by=5)
# 
# # Create a data frame for plotting the results
# plot_df <- data.frame(
#   x = x[subset_indices],
#   AR1 = y[subset_indices],
#   LP_Fitted = lp_prediction[subset_indices],
#   AR_Line = ar_line[subset_indices]
# )
# # Plot the AR(1) process, LP fitted values, and AR line using ggplot2
# ggplot(plot_df, aes(x=x)) +
#   geom_point(aes(y=AR1, color="AR(1) Process"), size=3, alpha=0.7) +   # AR(1) process as scatter points in light blue
#   geom_point(aes(y=LP_Fitted, color="LP Predictions"), size=3, alpha=0.7) +  # LP fitted, red points
#   geom_line(aes(y=AR_Line), color="black", size=1.2, linetype="solid") +  # AR line with slope 0.7
#   labs(x="x_(t-1)", y="x_t") +
#   scale_color_manual(values=c("lightblue", "red", "black")) +  # Define colors for AR(1), LP, and AR line
#   theme_minimal() +
#   theme(legend.position="top") +
#   guides(color=guide_legend(title="Legend"))  



sqrt(MSE(y,lp_prediction))
Hs4$Sigma2Hat
