library(horseshoe)
library(glmnet)
library(splines)
library(ggplot2)

set.seed(1991)
n <- 1000
x <- seq(0,1,length.out=n)
f <- function(x) 0.5*(abs(x-0.2))^5 + 2*(abs(x-0.5))^5 -0.5*(abs(x-0.9))^5
z <- f(x)
sd <- 0.05
y <- z + rnorm(n,0,sd)

fourier_basis <- function(x,K){
  basis <- matrix(1,nrow=length(x),ncol=2*K+1)
  for(k in 1:K){
    basis[,2*k]   <- sin(2*pi*k*x)
    basis[,2*k+1] <- cos(2*pi*k*x)
  }
  basis
}

J_seq <- seq(140,300,20)
mse_s_hs <- mse_s_ridge <- mse_f_hs <- mse_f_ridge <- numeric(length(J_seq))

for(i in seq_along(J_seq)){
  J <- J_seq[i]
  t <- seq(0.0001,0.9999,length.out=J)
  B_s <- bSpline(x,knots=t,degree=3,intercept=FALSE)
  Hs_s <- horseshoe(y,B_s,method.tau="halfCauchy",method.sigma="fixed",Sigma2=sd^2)
  theta_s_hs <- unlist(Hs_s$BetaHat)
  f_hat_s_hs <- B_s%*%theta_s_hs
  lambda_s <- cv.glmnet(B_s,y,alpha=0,nfolds=20)$lambda.min
  theta_s_n <- solve(crossprod(B_s)+lambda_s*diag(ncol(B_s)),crossprod(B_s,y))
  f_hat_s_n <- B_s%*%theta_s_n
  B_f <- fourier_basis(x,floor(J/2))
  Hs_f <- horseshoe(y,B_f,method.tau="halfCauchy",method.sigma="fixed",Sigma2=sd^2)
  theta_f_hs <- unlist(Hs_f$BetaHat)
  f_hat_f_hs <- B_f%*%theta_f_hs
  lambda_f <- cv.glmnet(B_f,y,alpha=0,nfolds=20)$lambda.min
  theta_f_n <- solve(crossprod(B_f)+lambda_f*diag(ncol(B_f)),crossprod(B_f,y))
  f_hat_f_n <- B_f%*%theta_f_n
  mse_s_hs[i]    <- mean((z-f_hat_s_hs)^2)
  mse_s_ridge[i] <- mean((z-f_hat_s_n)^2)
  mse_f_hs[i]    <- mean((z-f_hat_f_hs)^2)
  mse_f_ridge[i] <- mean((z-f_hat_f_n)^2)
}

mse_df <- data.frame(
  J             = J_seq,
  spline_hs     = mse_s_hs,
  spline_ridge  = mse_s_ridge,
  fourier_hs    = mse_f_hs,
  fourier_ridge = mse_f_ridge
)

library(writexl)
write_xlsx(mse_df, "mse_results_f3.xlsx")

