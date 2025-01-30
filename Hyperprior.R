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
y <- z+rnorm(n,mean=0, sd=0.4)

lambda <- 125
J_values <- 97:155


poisson_probs <- dpois(J_values, lambda)
poisson_probs <- poisson_probs / sum(poisson_probs)


MSE_values <- numeric(length(J_values))

for (j in seq_along(J_values)) {
  J <- J_values[j]
  
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J-4), 1.0001, 1.0002)

  B <- splineDesign(t, x, ord = 3, outer.ok = TRUE, sparse = TRUE)

  Hs2 <- horseshoesp(y, B, method.tau = "truncatedCauchy", Sigma2 = 0.16, 
                     tau = 0.4, s = 0.5, method.sigma = "fixed", 
                     burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
  
  beta_hat2 <- as.vector(unlist(Hs2[1]))

  y_hat <- as.numeric(B %*% beta_hat2)

  MSE_values[j] <- mean((y - y_hat)^2)
}

weighted_MSE <- sum(MSE_values * poisson_probs)

weighted_MSE
