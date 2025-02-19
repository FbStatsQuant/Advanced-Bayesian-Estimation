library(readxl)
library(splines)
library(horseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)
library(Matrix)
library(mvtnorm)
library(microbenchmark)

set.seed(2024)
n <- 10000
b <- 0.5
J <- floor(n^b)
set.seed(2024)
#J <- 10

x <- seq(0,1,length.out = n)

s <- function(m,x) {
  sqrt(2)*sum(seq(m)^(-1.5)*sin(seq(m))*cos(((seq(m))-0.5)*pi*x))
}
z <- vector()
for (j in 1:n) {
  z[j] = s(1000,j/n)
}

sd <- 0.4
y <- z+rnorm(n,mean=0, sd=sd)

t <- c(-0.0002, -0.0001, seq(0,1, length.out = J-4),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)
B1 <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                   sparse = FALSE)

X <- B
Q_star <- crossprod(X)

tau <- 1
sigma_sq <- 1
p   <- ncol(B)           
lambda <- rep(1, p)
lambda_star <- tau * lambda
r <- stats::rnorm(p)

C <- Cholesky((1 / sigma_sq) * (Q_star + Diagonal(x = 1 / (lambda_star^2))), 
              perm=FALSE, LDL=FALSE)

L_sparse <- sparseMatrix(i = as.integer(C@i) + 1, 
                         p = as.integer(C@p), 
                         x = as.numeric(C@x), 
                         dims = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])), 
                         index1 = TRUE)

L_sparse%*%t(L_sparse)-(1 / sigma_sq) * (Q_star + Diagonal(x = 1 / (lambda_star^2)))
mu_sparse <- solve(t(L_sparse),
                   solve(L_sparse,t(B)%*%y/sigma_sq))
u_sparse <- solve(t(L_sparse),r)
u_sparse+mu_sparse


Q_1 <- t(B1)%*%B1
L<- chol((1 / sigma_sq) * (Q_1 + Diagonal(x = 1 / (lambda_star^2))))
t(L)%*%L-(1 / sigma_sq) * (Q_1 + Diagonal(x = 1 / (lambda_star^2)))
v <- solve(t(L),t(t(y)%*%B)/sigma_sq)
mu <- solve(L,v)
u <- solve(L,r)
mu+u

microbenchmark(
  C <- Cholesky((1 / sigma_sq) * (Q_star + Diagonal(x = 1 / (lambda_star^2))), perm=FALSE, LDL=FALSE),
  L_sparse <- sparseMatrix(i = as.integer(C@i) + 1, 
                           p = as.integer(C@p), 
                           x = as.numeric(C@x), 
                           dims = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])), 
                           index1 = TRUE),
  mu_sparse <- solve(t(L_sparse),
                     solve(L_sparse,t(B)%*%y/sigma_sq)),
  u_sparse <- solve(t(L_sparse),r),
  as.numeric(mu_sparse+u_sparse)
)


microbenchmark(
  L<- chol((1 / sigma_sq) * (Q_1 + Diagonal(x = 1 / (lambda_star^2)))),
  v <- solve(t(L),t(t(y)%*%B)/sigma_sq),
  mu <- solve(L,v),
  u <- solve(L,r),
  as.numeric(mu + u)
)

I_n=diag(n)
l0=rep(0,p)
l1=rep(1,n)
l2=rep(1,p)


lambda_star=tau*lambda
U=(lambda_star^2)*t(X)
## step 1 ##
u=stats::rnorm(l2,l0,lambda_star)
v=as.numeric(X%*%u + stats::rnorm(n))
## step 2 ##
v_star=solve((X%*%U+I_n),((y/sqrt(sigma_sq))-v),type = "cg")
sqrt(sigma_sq)*(u+U%*%v_star)

microbenchmark(
  U=(lambda_star^2)*t(X),
  u=stats::rnorm(l2,l0,lambda_star),
  v=as.numeric(X%*%u + stats::rnorm(n)),
  v_star=solve((X%*%U+I_n),((y/sqrt(sigma_sq))-v),type = "cg"),
  sqrt(sigma_sq)*(u+U%*%v_star)
)



