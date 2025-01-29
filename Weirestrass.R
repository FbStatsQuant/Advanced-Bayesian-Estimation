library(splines)
library(MLmetrics)
library(Hmisc)
library(sparsehorseshoe)

set.seed(2024)
n <- 1000
b <- 0.7
J <- floor(n^b)

x <- seq(0,1,length.out = n)
a1 <- 0.5
a2 <- 7
s <- function(m,x) {
  sum(a1^(seq(m))*cos(a2^seq(10)*pi*x))
}
z <- vector()
for (j in 1:n) {
  z[j] = s(1000,j/n)
}
plot(x,z, xlab = "x",
     ylab = "y", col = "blue", type = "l")
y <- z+rnorm(n,mean=0, sd=0.4)
t <- c(-0.0002, -0.0001, seq(0,1, length.out = J-4),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs3 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = (0.16, tau=0.4, s=0.5,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)

beta_hat3 <- as.vector(unlist(Hs3[1]))

MSE(y,as.numeric(B%*%beta_hat3))

indices <- seq(1, length(x), by = 4)

plot(x[indices], (B%*%beta_hat3)[indices], cex.main = 1,
     xlab = "z",
     ylab = "x",
     col = "lightblue", pch = 21, bg = "lightblue")
lines(x[indices], z[indices], type = "l", col = "black")

legend("topleft",
       legend = c("Estimation", "Function f(z)"),
       col = c("lightblue", "black"),
       lty = c(1, 1))


##Grid and MSE
b_grid3 <- seq(0.7, 1.5, by = 0.05)
mse_values3 <- numeric(length(b_grid3))
J_values <- numeric(length(b_grid3))

for (i in seq_along(b_grid3)) {
  J_values[i] <- floor(n^(b_grid3[i]))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J_values[i]-4),1.0001, 1.0002)
  Sp_order <- 3
  B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                    sparse = TRUE)
  Hs3 <- horseshoesp(y, B, 
                     method.tau = "truncatedCauchy",
                     Sigma2 = 0.16,
                     tau = 0.4,
                     s = 0.5,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  beta_hat3 <- as.vector(unlist(Hs3[1]))
  mse_values3[i] <- MSE(y, as.numeric(B %*% beta_hat3))
  cat(sprintf("b=%.2f | J=%d | MSE=%.3f\n", b_grid3[i], J_values[i], mse_values3[i]))
}
results3 <- data.frame(
  J = J_values,
  MSE = mse_values3
)
print(results3)
