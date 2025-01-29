library(readxl)
library(splines)
library(sparsehorseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)

set.seed(123)

n <- 1000
b <- 0.95
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

y <- stationary_ar1[2:n]
x <- stationary_ar1[1:n-1]
t <- c(0.998*min(x), 0.999*min(x), seq(1,n, length.out = J-4),max(x)*1.001, max(x)*1.002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)
Hs4 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), s=0.5,  
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat4 <- as.vector(unlist(Hs4[1]))
sqrt(MSE(y,as.numeric(B%*%beta_hat4)))


#ARIMA
data <- ts(y)
best_model <- auto.arima(data)
summary(best_model)


##Grid and MSE
b_grid4 <- seq(0.7, 1.5, by = 0.05)
mse_values4 <- numeric(length(b_grid4))
J_values <- numeric(length(b_grid4))

for (i in seq_along(b_grid4)) {
  J_values[i] <- floor(n^(b_grid4[i]))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J_values[i]-4),1.0001, 1.0002)
  Sp_order <- 3
  B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                    sparse = TRUE)
  Hs4 <- horseshoesp(y, B, 
                     method.tau = "truncatedCauchy",
                     Sigma2 = 0.16,
                     tau = 0.4,
                     s = 0.5,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  beta_hat4 <- as.vector(unlist(Hs4[1]))
  mse_values4[i] <- MSE(y, as.numeric(B %*% beta_hat4))
  cat(sprintf("b=%.2f | J=%d | MSE=%.3f\n", b_grid4[i], J_values[i], mse_values4[i]))
}
results4 <- data.frame(
  J = J_values,
  MSE = mse_values4
)
print(results4)

