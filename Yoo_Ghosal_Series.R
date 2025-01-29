library(splines)
library(MLmetrics)
library(Hmisc)
library(sparsehorseshoe)


n <- 1000
b <- 0.7
set.seed(2024)
J <- floor(n^b)
x <- seq(0,1,length.out = n)
s <- function(m,x) {
  sqrt(2)*sum(seq(m)^(-1.5)*sin(seq(m))*cos(((seq(m))-0.5)*pi*x))
}
z <- vector()
for (j in 1:n) {
  z[j] = s(1000,j/n)
}

plot(x,z, xlab = "x",
     ylab = "y", col = "blue", type = "l")
y <- z+rnorm(n,mean=0, sd=0.4)
plot(x, y, xlab = "x", ylab = "y", col = "lightblue", type = "l",
     main = "Function vs Random Points", cex.main = 2)
lines(x, z, type = "l", col = "black")

legend("topright",                       
       legend = c("Random Points", "Actual Function"),  
       col = c("lightblue", "black"),   
       lty = 1,                          
       bty = "n")                        

t <- c(-0.0002, -0.0001, seq(0,1, length.out = J-4),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs2 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = 0.16, tau=0.4, s=0.5,
                  method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)

beta_hat2 <- as.vector(unlist(Hs2[1]))

MSE(y,as.numeric(B%*%beta_hat2))


plot(x,B%*%beta_hat2, xlab = "z",
     ylab = "x", col = "lightblue", cex.main = 2)
lines(x,z, type = "l",
      col = "black")

legend("topright",                       
       legend = c("Estimation", "Actual Function"),
       col = c("lightblue", "black"),
       lty = 1,                      
       bty = "n")                    


##Grid and MSE
b_grid2 <- seq(0.7, 1.5, by = 0.05)
mse_values2 <- numeric(length(b_grid2))
J_values <- numeric(length(b_grid2))

for (i in seq_along(b_grid2)) {
  J_values[i] <- floor(n^(b_grid2[i]))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J_values[i]-4),1.0001, 1.0002)
  Sp_order <- 3
  B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                    sparse = TRUE)
  Hs2 <- horseshoesp(y, B, 
                     method.tau = "truncatedCauchy",
                     Sigma2 = 0.16,
                     tau = 0.4,
                     s = 0.5,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  beta_hat2 <- as.vector(unlist(Hs2[1]))
  mse_values2[i] <- MSE(y, as.numeric(B %*% beta_hat2))
  cat(sprintf("b=%.2f | J=%d | MSE=%.3f\n", b_grid2[i], J_values[i], mse_values2[i]))
}
results2 <- data.frame(
  J = J_values,
  MSE = mse_values2
)
print(results2)


