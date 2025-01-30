library(splines)
library(sparsehorseshoe)
library(DiceKriging)

#Splines 

n <- 2000
b <- 0.7
set.seed(2024)
J <- floor(n^b)

x <- seq(0, 1, length.out = n)

s <- function(m, x) {
  sqrt(2) * sum(seq(m)^(-1.5) * sin(seq(m)) * cos(((seq(m)) - 0.5) * pi * x))
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

predictions <- B_test %*% Hs_model$BetaHat  

mse_sp <- mean((test_data$y - as.numeric(predictions))^2)
print(paste("Mean Squared Error (MSE):", mse_sp))

##GPs

gpr_model <- km(
  y ~ 1, design = train_data[, "x", drop = FALSE], response = train_data$y, 
  covtype = "matern5_2", nugget = 1e-2
)

predictions <- predict(gpr_model, 
                       newdata = test_data[, "x", drop = FALSE], type = "SK")$mean

mse_gp <- mean((test_data$y - predictions)^2)
print(paste("Mean Squared Error (MSE):", mse_gp))
