library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)
library(writexl)


set.seed(1991)
n <- 1000
x <- seq(0,1, length.out = n)
sd <- 0.3

a1 <- 0.5
a2 <- 10
a1*a2 > 1 + 1.5*pi
alpha <- -log(a1)/log(a2)
s <- function(m,x) {
  sum(a1^(seq(m))*cos(a2^seq(10)*pi*x))
}

z <- sapply(x, function(xi) s(1000, xi))
y <- z + rnorm(n, mean = 0, sd = sd)

J_values <- seq(50, 500, by = 50)
mse_results <- data.frame(J = J_values, mse_s_hs = NA, mse_s_ridge = NA, mse_f_hs = NA, 
                          mse_f_ridge = NA, mse_l_hs = NA, mse_l_ridge = NA)

for (J in J_values) {
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J),1.0001, 1.0002)
  Sp_order <- 3
  B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = FALSE)
  Hs_s <- horseshoe(y, B_s, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
  theta_s_hs <- unlist(Hs_s$BetaHat)
  f_hat_s_hs <- B_s%*%theta_s_hs
  lambda_s <- cv.glmnet(B_s,z, alpha=0, nfolds = 20)$lambda.min
  precision_matrix <- (lambda_s/((sd^2)) * diag(ncol(B_s)) + (1/(sd^2)) * crossprod(B_s))
  theta_s_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_s, y))
  f_hat_s_n <- B_s%*%theta_s_n
  
  fourier_basis <- function(x, K) {
    
    basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
    
    for (m in 1:K) {
      basis[, 2 * m]   <- sin(2 * pi * m * x)
      basis[, 2 * m + 1] <- cos(2 * pi * m * x)
    }
    
    return(basis)
  }
  
  B_f <- fourier_basis(x, floor(J/2))
  
  Hs_f <- horseshoe(y, B_f, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
  theta_f_hs <- unlist(Hs_f$BetaHat)
  f_hat_f_hs <- B_f%*%theta_f_hs
  lambda_f <- cv.glmnet(B_f,z, alpha=0, nfolds = 20)$lambda.min
  precision_matrix <- (lambda_f/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
  theta_f_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
  f_hat_f_n <- B_f%*%theta_f_n
  
  x_scaled <- 2 * x - 1
  B_l <- matrix(0, nrow = length(x), ncol = J+1)
  B_l[, 1] <- 1
  if(J >= 2) B_l[, 2] <- x_scaled
  for (k in 3:(J+1)) {
    m <- k - 2  
    B_l[, k] <- ((2*m - 1) * x_scaled * B_l[, k-1] - (m - 1) * B_l[, k-2]) / m
  }
  Hs_l <- horseshoe(y, B_l, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
  theta_l_hs <- unlist(Hs_l$BetaHat)
  f_hat_l_hs <- B_l%*%theta_l_hs
  lambda_l <- cv.glmnet(B_l,z, alpha=0, nfolds = 20)$lambda.min
  precision_matrix <- (lambda_l/((sd^2)) * diag(ncol(B_l)) + (1/(sd^2)) * crossprod(B_l))
  theta_l_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_l, y))
  f_hat_l_n <- B_l%*%theta_l_n
  
  mse_results[mse_results$J == J, ] <- c(J, mean((z - f_hat_s_hs)^2), mean((z - f_hat_s_n)^2), mean((z - f_hat_f_hs)^2), mean((z - f_hat_f_n)^2), mean((z - f_hat_l_hs)^2), mean((z - f_hat_l_n)^2))
}

write_xlsx(mse_results, "MSE_W.xlsx")

