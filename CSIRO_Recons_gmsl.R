library(splines)
library(MLmetrics)
library(Hmisc)
library(sparsehorseshoe)


set.seed(2024)
setwd("C:/Users/felip/OneDrive - Rice University/Overcomplete random series priors/Coding")

df <- read.table("CSIRO_Recons_gmsl_yr_2015.txt"
                 , header = FALSE, sep = "", dec = ".")
x <- unlist(df["V1"])
y <- unlist(df["V2"])
n <- length(x)
b <- 0.9
J <- floor(n^b)

t <- c(1880.48, 1880.49, seq(1880.5,2013.5, length.out = J-4),2013.51, 2013.52)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs5 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), s=0.3,  
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat5 <- as.vector(unlist(Hs5[1]))

MSE(y,as.numeric(B%*%beta_hat5))
plot(x,y, xlab = "Years",
     ylab = " Millimeters", col = "black", cex.main = 2)
points(x,B%*%beta_hat5, type = "l",
      col = "red")
# Add a legend
legend("topleft", 
       legend = c("Estimation", "Function f(t)"), 
       col = c("red", "black"), 
       lty = c(1, 1))

