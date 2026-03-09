library(horseshoe)

set.seed(1991)

sd  <- 0.15
b   <- 0.8
f   <- function(x) abs(x - 0.5)

fourier_basis <- function(x, K) {
  basis <- matrix(1, nrow = length(x), ncol = K + 1)
  for (k in 1:K) basis[, k + 1] <- cos(pi * k * x)
  basis
}

n_seq <- seq(1000, 10000, by = 1000)
a_n   <- numeric(length(n_seq))

for (i in seq_along(n_seq)) {

  n   <- n_seq[i]
  J   <- floor(n^b)
  x   <- seq(0, 1, length.out = n)
  z   <- f(x)
  y   <- z + rnorm(n, mean = 0, sd = sd)

  B   <- fourier_basis(x, J)

  Hs  <- horseshoe(y, B,
                   method.tau   = "fixed", tau = 1 / J,
                   method.sigma = "fixed", Sigma2 = sd^2,
                   nmc = 4000, burn = 1000)

  theta      <- apply(Hs$BetaSamples, 1, median)
  a_n[i]     <- sum(abs(theta) > 1 / (J*n^(0.2)))

  cat(sprintf("n = %5d | J = %4d | nonzeros = %d\n", n, J, a_n[i]))
}

# Log-log regression: log(a_n) = log(A) + a * log(n)
fit   <- lm(log(a_n) ~ log(n_seq))
a_hat <- coef(fit)[["log(n_seq)"]]

cat(sprintf("\nEstimated polynomial exponent a_hat = %.4f\n", a_hat))

# Plot
plot(log(n_seq), log(a_n),
     xlab = "log(n)", ylab = "log(a_n)",
     main = "Log-log plot of nonzero coefficients vs n",
     pch = 16)
abline(fit, col = "steelblue", lwd = 2)
legend("topleft",
       legend = sprintf("slope = %.4f", a_hat),
       col = "steelblue", lwd = 2, bty = "n")
