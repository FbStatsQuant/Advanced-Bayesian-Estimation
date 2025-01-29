library(readxl)
library(splines)
library(horseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)

set.seed(2024)
setwd("C:/Users/felip/OneDrive - Rice University/Overcomplete random series priors/Coding")
df <- read.csv("SP500.csv", header = TRUE, sep = ",")
y <- rev(as.vector(unlist(df['Close'])))
         
n <- length(y)
x <- c(1:length(y))
b <- 0.6
J <- floor(n^b)


t <- c(0.98*min(x), 0.99*min(x), seq(1,n, length.out = J-4),max(x)*1.01, max(x)*1.02)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = FALSE)
Hs <- horseshoe(y, B, method.tau = c("truncatedCauchy"),  
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat <- as.vector(unlist(Hs[1]))

sqrt(MSE(y,B%*%beta_hat))

plot(x,y, xlab = "Observarions",
     ylab = "Closing Value", col = "lightblue", cex.main = 2)
points(x,B%*%beta_hat, type = "l",
       col = "red")
legend("bottomright",
       legend = c("Estimation", "Actual values"),
       col = c("red", "black"),
       lty = c(1, 1))
