library(splines)
library(MLmetrics)
library(Hmisc)
library(sparsehorseshoe)


set.seed(2024)
n <- 1000
J <- 20


x <- seq(0,1,length.out = n)
f <- function(x){
  exp(-2*x^2)}
y <- f(x)+rnorm(n,mean=0, sd=0.4)
plot(x, y, xlab = "x", ylab = "y", col = "lightblue", type = "l",cex.main = 2)
lines(x, f(x), type = "l", col = "black")

legend("topright",                       
       legend = c("Random Points", "Actual Function"),  
       col = c("lightblue", "black"),   
       lty = 1,                          
       bty = "n")
t <- c(-0.0002, -0.0001, seq(0,1, length.out = J-4),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs1 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = (0.4)^2, tau=0.4,
                  method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat1 <- as.vector(unlist(Hs1[1]))

MSE(y,as.numeric(B%*%beta_hat1))

plot(x,f(x), xlab = "z",
     ylab = "x", col = "lightblue", cex.main = 2)
lines(x,B%*%beta_hat1, type = "l",
      col = "black")
# Adding a legend
legend("topright",                       # Position of the legend
       legend = c("Estimation", "Actual Function"),  # Labels for the legend
       col = c("lightblue", "black"),   # Colors corresponding to the lines
       lty = 1,                          # Line type
       bty = "n")                        # No box around the legend