library(readxl)
library(splines)
library(horseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)

setwd("C:/Users/felip/OneDrive - Rice University/Overcomplete random series priors/Coding")
df <- read.csv("TSLA2.csv", header = TRUE, sep = ",")
y <- rev(as.vector(unlist(df['Close'])))
set.seed(2024)
length(y)
x <- c(1:length(y))

n <- length(x)
b <- 0.8
J <- floor(n^b)

t <- c(0.98*min(x[1:n]), 0.99*min(x[1:n]),
       seq(min(x[1:n]),max(x[1:n]), length.out = J-4),1.01*max(x[1:n]), 1.02*max(x[1:n]))
Sp_order <- 3
B <- splineDesign(t, x[1:n], ord = Sp_order, outer.ok = TRUE,
                  sparse = FALSE)
Hs <- horseshoe(y[1:n], B, method.tau = c("halfCauchy"),
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat <- as.vector(unlist(Hs[1]))
c(J,sum(HS.var.select(Hs, y[1:n], method = "intervals"))) #selected betas

MSE(y[1:n],B%*%beta_hat)

plot(x[1:n],y[1:n], xlab = "Observarions",
     ylab = "Closing Value", col = "black", cex.main=2)
points(x[1:n],B%*%beta_hat, type = "l",
       col = "red")
legend("bottomright",
       legend = c("Estimation", "Actual values"),
       col = c("red", "black"),
       lty = c(1, 1))



beta_hat_t1 <- vector()
B_t1 <- matrix()
MSE_test1 <- vector()
y_pred <- vector()

m <- 202
for (i in 1:50){
  B_t1 <- splineDesign(t, x[1:(m+i)], ord = Sp_order, outer.ok = TRUE,
                      sparse = FALSE)
  Hs_t <- horseshoe(y[1:(m+i)], B_t1, method.tau = c("halfCauchy"),
                    method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
  beta_hat1_t <- as.vector(unlist(Hs_t[1]))
  y_pred[i] <- B_t1[(m+i),]%*%beta_hat1_t
  MSE_test1[i] <- MSE(y[(m+i)],y_pred[i])
  print(MSE_test1)
  t <- c(min(x[1:(m+i)])*0.98, min(x[1:(m+i)])*0.99, seq(min(x[1:(m+i)]),max(x[1:(m+i)]), length.out = J-4),max(x[1:(m+i)])*1.01,
          max(x[1:(m+i)])*1.02)
}

indices <- (m+1):252

plot(y[indices],y_pred, xlab = 'Actual', ylab = 'Predicted')
max_val <- max(c(y[indices], y_pred[indices]))
min_val <- min(c(y[indices], y_pred[indices]))
abline(a = 0, b = 1, col = "red", lty = 2)  # Identity line (y = x)

sqrt(mean(MSE_test1))


x1 <- y[1:(n-1)]
x2 <- y[2:n]



t1 <- c(min(x1)*0.98, min(x1)*0.99, seq(min(x1),max(x1), length.out = J-4),max(x1)*1.01,
        max(x1)*1.02)
Sp_order <- 3
B1 <- splineDesign(t1, x1, ord = Sp_order, outer.ok = TRUE,
                   sparse = FALSE)
Hs1 <- horseshoe(x2, B1, method.tau = c("halfCauchy"),
                 method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat1 <- as.vector(unlist(Hs1[1]))
c(J,sum(HS.var.select(Hs1,x2, method = c("intervals"), threshold = 0.05))) #selected betas

sqrt(MSE(x2,B1%*%beta_hat1))

#Test


t2 <- c(min(x1[1:m])*0.98, min(x1[1:m])*0.99, seq(min(x1[1:m]),max(x1[1:m]), length.out = J-4),max(x1[1:m])*1.01,
        max(x1[1:m])*1.02)
Sp_order <- 3
B2 <- splineDesign(t2, x1[1:m], ord = Sp_order, outer.ok = TRUE,
                   sparse = FALSE)
Hs1 <- horseshoe(x2[1:m], B2, method.tau = c("halfCauchy"),
                 method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat1 <- as.vector(unlist(Hs1[1]))
c(J,sum(HS.var.select(Hs1,x2, method = c("intervals"), threshold = 0.05))) #selected betas

sqrt(MSE(x2,B1%*%beta_hat1))

beta_hat_t <- vector()
B_t <- matrix()
MSE_test <- vector()
x2_pred <- vector()


for (i in 1:50){
  B_t <- splineDesign(t2, x1[1:(m+i)], ord = Sp_order, outer.ok = TRUE,
                      sparse = FALSE)
  Hs_t <- horseshoe(x2[1:(m+i)], B_t, method.tau = c("halfCauchy"),
                   method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
  beta_hat1_t <- as.vector(unlist(Hs_t[1]))
  x2_pred[i] <- B_t[(m+i),]%*%beta_hat1_t
  MSE_test[i] <- MSE(x2[m+i],x2_pred[i])
  print(MSE_test)
  t2 <- c(min(x1[1:(m+i)])*0.98, min(x1[1:(m+i)])*0.99, seq(min(x1[1:(m+i)]),max(x1[1:(m+i)]), length.out = J-4),max(x1[1:(m+i)])*1.01,
          max(x1[1:(m+i)])*1.02)
}



sqrt(mean(MSE_test))

indices <- (m+1):252

plot(x2[indices],x2_pred, xlab = 'Actual', ylab = 'Predicted')
max_val <- max(c(x2[indices], x2_pred[indices]))
min_val <- min(c(x2[indices], x2_pred[indices]))
abline(a = 0, b = 1, col = "red", lty = 2)  # Identity line (y = x)

#best ARIMA

forecast_value <- vector()
MSE_arima <- vector()

for (i in 1:50){
  data <- ts(y[1:(m+i-1)])
  best_model <- auto.arima(data)
  forecast_value[i] <- forecast(best_model, h=1)$mean[1]
  MSE_arima[i] <- MSE(forecast_value[i],y[m+i])
}
sqrt(mean(MSE_arima))


plot(y[indices],forecast_value, xlab = 'Actual', ylab = 'Predicted',pch = 16,
      cex.main=2)
points(x2[indices],x2_pred, col = "red",pch = 16)
max_val <- max(c(y[indices], forecast_value[indices]))
min_val <- min(c(y[indices], forecast_value[indices]))
abline(a = 0, b = 1, col = "blue", lty = 2)  # Identity line (y = x)
