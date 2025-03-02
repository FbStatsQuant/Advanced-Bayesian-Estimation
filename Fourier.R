library(ggplot2)
library(horseshoe)
library(glmnet)
library(MLmetrics)


set.seed(1991)
n <- 1000
x <- seq(0,1, length.out=n)

J_true <- 10
theta_true <- rnorm(J_true, mean = 0, sd = 4)

fourier_basis <- function(x, K) {
  
  basis_matrix <- matrix(0, nrow = length(x), ncol = 2 * K)
  for (k in 1:K) {
    basis_matrix[, 2 * k - 1] <- sin(2 * pi * k * x)
    basis_matrix[, 2 * k] <- cos(2 * pi * k * x)
  }
  
  return(basis_matrix)
}

B_true <- fourier_basis(x, floor(J_true/2))
f_true <- B_true%*%theta_true
sd <- 7
y <- f_true + rnorm(n, mean=0, sd=sd)
mean(f_true^2)/mean((f_true-y)^2)


data <- data.frame(x, f_true, y)

ggplot(data, aes(x = x)) +
  geom_line(aes(y = f_true), color = "black", linewidth = 1.0) +
  geom_point(aes(y = y), color = "blue", size = 0.5, alpha = 0.5) +
  labs(title = "True Function and random points", x = "x", y = "y") +
  theme(plot.title = element_text(size = 16)) +
  scale_color_manual(values = c("lightblue", "black")) +
  theme_minimal()



B_50 <- fourier_basis(x, 25)

Hs_50 <- horseshoe(y, B_50, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
theta_hs_50 <- as.vector(unlist(Hs_50[1]))

k_50 <- cv.glmnet(B_50,f_true, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (1/((sd^2)/k_50) * diag(ncol(B_50)) + (1/(sd^2)) * crossprod(B_50))
theta_nr_50 <- as.vector(solve(precision_matrix, (1/(sd^2)) * crossprod(B_50, y)))

theta_true
theta_hs_50
theta_nr_50

MSE(f_true, B_50%*%theta_hs_50)
MSE(f_true, B_50%*%theta_nr_50)


B_100 <- fourier_basis(x, 50)

Hs_100 <- horseshoe(y, B_100, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
theta_hs_100 <- as.vector(unlist(Hs_100[1]))

k_100 <- cv.glmnet(B_100,f_true, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (1/((sd^2)/k_100) * diag(ncol(B_100)) + (1/(sd^2)) * crossprod(B_100))
theta_nr_100 <- as.vector(solve(precision_matrix, (1/(sd^2)) * crossprod(B_100, y)))

theta_true
theta_hs_100
theta_nr_100

MSE(f_true, B_100%*%theta_hs_100)
MSE(f_true, B_100%*%theta_nr_100)


B_200 <- fourier_basis(x, 100)

Hs_200 <- horseshoe(y, B_200, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
theta_hs_200 <- as.vector(unlist(Hs_200[1]))

k_200 <- cv.glmnet(B_200,f_true, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (1/((sd^2)/k_200) * diag(ncol(B_200)) + (1/(sd^2)) * crossprod(B_200))
theta_nr_200 <- as.vector(solve(precision_matrix, (1/(sd^2)) * crossprod(B_200, y)))

theta_true
theta_hs_200
theta_nr_200

MSE(f_true, B_200%*%theta_hs_200)
MSE(f_true, B_200%*%theta_nr_200)



B_400 <- fourier_basis(x, 200)

Hs_400 <- horseshoe(y, B_400, method.tau = c("halfCauchy"), method.sigma = c("fixed"), Sigma2 = sd^2)
theta_hs_400 <- as.vector(unlist(Hs_400[1]))

k_400 <- cv.glmnet(B_400,f_true, alpha=0, nfolds = 20)$lambda.min
precision_matrix <- (1/((sd^2)/k_400) * diag(ncol(B_400)) + (1/(sd^2)) * crossprod(B_400))
theta_nr_400 <- as.vector(solve(precision_matrix, (1/(sd^2)) * crossprod(B_400, y)))

theta_true
theta_hs_400
theta_nr_400

MSE(f_true, B_400%*%theta_hs_400)
MSE(f_true, B_400%*%theta_nr_400)

