library(splines)
library(sparsehorseshoe)
library(DiceKriging)
library(ggplot2)
library(tidyr)
library(dplyr)

#Splines 

n <- 2000
b <- 0.6
set.seed(2024)
J <- floor(n^b)
x <- seq(0, 1, length.out = n)

a1 <- 0.5
a2 <- 7
s <- function(m,x) {
  sum(a1^(seq(m))*cos(a2^seq(10)*pi*x))
}

z <- sapply(x, function(j) s(1000, j))

y <- z + rnorm(n, mean = 0, sd = 0.4)

data <- data.frame(x = x, y = y)

t <- c(-0.0002, -0.0001, seq(0, 1, length.out = J - 4), 1.0001, 1.0002)
Sp_order <- 3

set.seed(123)
train_idx <- sample(1:n, size = round(0.8 * n))
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

B_train <- splineDesign(t, train_data$x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
B_test <- splineDesign(t, test_data$x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)

Hs_model <- horseshoesp(
  y = train_data$y, X = B_train, 
  method.tau = "truncatedCauchy", Sigma2 = 0.16, tau = 0.4, s = 0.5,
  method.sigma = "fixed", burn = 1000, nmc = 5000, thin = 1, alpha = 0.05
)

predictions_lp <- B_test %*% Hs_model$BetaHat  

mse_2_lp <- mean((test_data$y - as.numeric(predictions_lp))^2)
print(paste("Mean Squared Error (MSE):", mse_2_lp))

##GPs

gpr_model <- km(
  y ~ 1, design = train_data[, "x", drop = FALSE], response = train_data$y, 
  covtype = "matern5_2", nugget = 0.1
)

predictions_gp <- predict(gpr_model, 
                       newdata = test_data[, "x", drop = FALSE], type = "SK")$mean

mse_2_gp <- mean((test_data$y - predictions_gp)^2)
print(paste("Mean Squared Error (MSE):", mse_2_gp))

# Add predictions to the test data
test_data <- test_data %>%
  mutate(
    pred_spline = as.numeric(predictions_lp),
    pred_gp = predictions_gp
  )

# Reshape data into long format for ggplot
plot_data <- test_data %>%
  select(x, y, pred_spline, pred_gp) %>%
  pivot_longer(cols = c(y, pred_spline, pred_gp), 
               names_to = "source", values_to = "value") %>%
  mutate(source = recode(source, 
                         "y" = "Original Data", 
                         "pred_spline" = "Splines Estimate", 
                         "pred_gp" = "GP Estimate"))

# Create ggplot
ggplot(plot_data, aes(x = x, y = value, color = source)) +
  geom_point(data = filter(plot_data, source == "Original Data"), 
             size = 0.5, alpha = 0.5) +  # Original data as points
  geom_line(data = filter(plot_data, source != "Original Data"), 
            size = 1) +  # Estimates as lines
  labs(title = "Comparison of Original Data, LP and GP Estimates",
       x = "x",
       y = "y",
       color = "Data Source") +
  theme_minimal()

