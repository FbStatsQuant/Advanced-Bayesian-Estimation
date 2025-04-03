library(horseshoe)
library(ggplot2)
library(splines)
library(glmnet)
library(writexl)
library(readxl)

n_values <- seq(1000, 20000, by = 1000)
MSE_hs_b <- rep(0, length(n_values))
MSE_hs_f <- rep(0, length(n_values))
MSE_hs_l <- rep(0, length(n_values))
MSE_no_b <- rep(0, length(n_values))
MSE_no_f <- rep(0, length(n_values))
MSE_no_l <- rep(0, length(n_values))
b <- 0.7
sd <- 0.3

s <- function(m, x) {
  sqrt(2) * sum(seq(m)^(-1.5) * sin(seq(m)) * cos(((seq(m)) - 0.5) * pi * x))
}


#Hs Splines


for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J+1),1.0001, 1.0002)
  Sp_order <- 3
  B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  Hs <- horseshoesp(y, B_s, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                    Sigma2 = sd^2)
  
  f_hat_hs <- as.numeric(B_s %*% as.vector(unlist(Hs$BetaHat)))
  MSE_hs_b[i] <- mean((z - f_hat_hs)^2)
}

my_data_hs_b <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_hs_b)
)

fit_hs_bs <- lm(log_MSE ~ log_n, data = my_data_hs_b)
summary_fit <- summary(fit_hs_bs)

slope <- round(coef(fit_hs_bs)[2], 3)

ggplot(my_data_hs_b, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): B-splines / Horseshoe",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_hs_b$log_n)  +2, y = min(my_data_hs_b$log_MSE) + 0.1,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()



#Hs Fourier



fourier_basis <- function(x, K) {
  
  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
  
  for (m in 1:K) {
    basis[, 2 * m]   <- sin(2 * pi * m * x)
    basis[, 2 * m + 1] <- cos(2 * pi * m * x)
  }
  
  return(basis)
}

for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = floor(n^b)+1),1.0001, 1.0002)
  Sp_order <- 3
  B_f <- fourier_basis(x, floor(J/2))
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  Hs <- horseshoesp(y, B_f, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                    Sigma2 = sd^2)
  
  f_hat_hs <- as.numeric(B_f %*% as.vector(unlist(Hs$BetaHat)))
  MSE_hs_f[i] <- mean((z - f_hat_hs)^2)
}

my_data_hs_f <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_hs_f)
)

fit_hs_f <- lm(log_MSE ~ log_n, data = my_data_hs_f)
summary_fit <- summary(fit_hs_f)

slope <- round(coef(fit_hs_f)[2], 3)

ggplot(my_data_hs_f, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): Fourier / Horseshoe",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_hs_f$log_n)  +2, y = min(my_data_hs_f$log_MSE) + 0.3,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()



#Hs Legendre



for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  x_scaled <- 2 * x - 1
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = floor(n^b)+1),1.0001, 1.0002)
  Sp_order <- 3
  B_l <- matrix(0, nrow = length(x), ncol = J+1)
  B_l[, 1] <- 1
  if(J >= 2) B_l[, 2] <- x_scaled
  for (k in 3:(J+1)) {
    m <- k - 2  
    B_l[, k] <- ((2*m - 1) * x_scaled * B_l[, k-1] - (m - 1) * B_l[, k-2]) / m
  }
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  Hs <- horseshoesp(y, B_l, method.tau = c("halfCauchy"), method.sigma = c("fixed"),
                    Sigma2 = sd^2)
  
  f_hat_hs <- as.numeric(B_l %*% as.vector(unlist(Hs$BetaHat)))
  MSE_hs_l[i] <- mean((z - f_hat_hs)^2)
}

my_data_hs_l <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_hs_l)
)

fit_hs_l <- lm(log_MSE ~ log_n, data = my_data_hs_l)
summary_fit <- summary(fit_hs_l)

slope <- round(coef(fit_hs_l)[2], 3)

ggplot(my_data_hs_l, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): Legendre / Horseshoe",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_hs_l$log_n)  +2, y = min(my_data_hs_l$log_MSE) + 0.5,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()



#Normal Splines



for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J+1),1.0001, 1.0002)
  Sp_order <- 3
  B_s <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = FALSE)
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  lambda <- cv.glmnet(B_s,z, alpha=0, nfolds = 10)$lambda.min
  precision_matrix <- (lambda/((sd^2)) * diag(ncol(B_s)) + (1/(sd^2)) * crossprod(B_s))
  theta_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_s, y))
  f_hat_n <- B_s%*%theta_n
  MSE_no_b[i] <- mean((z - f_hat_n)^2)
}

my_data_no_b <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_no_b)
)

fit_no_bs <- lm(log_MSE ~ log_n, data = my_data_no_b)
summary_fit <- summary(fit_no_bs)

slope <- round(coef(fit_no_bs)[2], 3)

ggplot(my_data_no_b, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): B-splines / Normal",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_no_b$log_n)  +2, y = min(my_data_no_b$log_MSE) + 0.1,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()


#Normal Fourier


fourier_basis <- function(x, K) {
  
  basis <- matrix(1, nrow = n, ncol = 2 * K + 1)
  
  for (m in 1:K) {
    basis[, 2 * m]   <- sin(2 * pi * m * x)
    basis[, 2 * m + 1] <- cos(2 * pi * m * x)
  }
  
  return(basis)
}

for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = floor(n^b)+1),1.0001, 1.0002)
  Sp_order <- 3
  B_f <- fourier_basis(x, floor(J/2))
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  
  lambda <- cv.glmnet(B_f,z, alpha=0, nfolds = 10)$lambda.min
  precision_matrix <- (lambda/((sd^2)) * diag(ncol(B_f)) + (1/(sd^2)) * crossprod(B_f))
  theta_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_f, y))
  f_hat_n <- B_f%*%theta_n
  MSE_no_f[i] <- mean((z - f_hat_n)^2)
}

my_data_no_f <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_no_f)
)

fit_no_f <- lm(log_MSE ~ log_n, data = my_data_no_f)
summary_fit <- summary(fit_no_f)

slope <- round(coef(fit_no_f)[2], 3)

ggplot(my_data_no_f, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): Fourier / Normal",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_no_f$log_n)  +2, y = min(my_data_no_f$log_MSE) + 0.3,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()



#Normal Legendre



for (i in seq_along(n_values)) {
  n <- n_values[i]
  J <- floor(n^b)
  set.seed(n)
  x <- seq(0,1, length.out = n)
  x_scaled <- 2 * x - 1
  z <- sapply(x, function(xi) s(1000, xi))
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = floor(n^b)+1),1.0001, 1.0002)
  Sp_order <- 3
  B_l <- matrix(0, nrow = length(x), ncol = J+1)
  B_l[, 1] <- 1
  if(J >= 2) B_l[, 2] <- x_scaled
  for (k in 3:(J+1)) {
    m <- k - 2  
    B_l[, k] <- ((2*m - 1) * x_scaled * B_l[, k-1] - (m - 1) * B_l[, k-2]) / m
  }
  y <- z + rnorm(n, mean = 0, sd = sd)
  
  lambda <- cv.glmnet(B_l,z, alpha=0, nfolds = 10)$lambda.min
  precision_matrix <- (lambda/((sd^2)) * diag(ncol(B_l)) + (1/(sd^2)) * crossprod(B_l))
  theta_n <- solve(precision_matrix, (1/(sd^2)) * crossprod(B_l, y))
  f_hat_n <- B_l%*%theta_n
  MSE_no_l[i] <- mean((z - f_hat_n)^2)
}

my_data_no_l <- data.frame(
  log_n = log(n_values),
  log_MSE = log(MSE_no_l)
)

fit_no_l <- lm(log_MSE ~ log_n, data = my_data_no_l)
summary_fit <- summary(fit_no_l)

slope <- round(coef(fit_no_l)[2], 3)

ggplot(my_data_no_l, aes(x = log_n, y = log_MSE)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Regression of log(MSE) on log(n): Legendre / Normal",
       x = "log(n)", y = "log(MSE)") +
  annotate("text", x = min(my_data_no_l$log_n)  +2, y = min(my_data_no_l$log_MSE) + 0.5,
           label = paste("Slope =", slope), hjust = 0, size = 5, color = "black") +
  theme_minimal()



df_mse <- data.frame(
  n = n_values,
  MSE_hs_bspline = MSE_hs_b,
  MSE_hs_fourier = MSE_hs_f,
  MSE_hs_legendre = MSE_hs_l,
  MSE_lasso_bspline = MSE_no_b,
  MSE_lasso_fourier = MSE_no_f,
  MSE_lasso_legendre = MSE_no_l
)

write_xlsx(df_mse, path = "MSE_results.xlsx")


df_mse <- read_xlsx("MSE_results.xlsx")
df_mse$log_n <- log(df_mse$n)


df_mse$log_MSE_hs_bspline     <- log(df_mse$MSE_hs_bspline)
df_mse$log_MSE_hs_fourier     <- log(df_mse$MSE_hs_fourier)
df_mse$log_MSE_hs_legendre    <- log(df_mse$MSE_hs_legendre)
df_mse$log_MSE_lasso_bspline  <- log(df_mse$MSE_lasso_bspline)
df_mse$log_MSE_lasso_fourier  <- log(df_mse$MSE_lasso_fourier)
df_mse$log_MSE_lasso_legendre <- log(df_mse$MSE_lasso_legendre)
df_mse$log_log_n <- log(df_mse$log_n)

cat("\n--- Regressions: log(MSE) ~ log(n) ---\n")
summary(lm(log_MSE_hs_bspline     ~ log_n, data = df_mse))
summary(lm(log_MSE_hs_fourier     ~ log_n, data = df_mse))
summary(lm(log_MSE_hs_legendre    ~ log_n, data = df_mse))
summary(lm(log_MSE_lasso_bspline  ~ log_n, data = df_mse))
summary(lm(log_MSE_lasso_fourier  ~ log_n, data = df_mse))
summary(lm(log_MSE_lasso_legendre ~ log_n, data = df_mse))


## Last 10

df_last10 <- tail(df_mse, 10)

df_last10$log_MSE_hs_bspline     <- log(df_last10$MSE_hs_bspline)
df_last10$log_MSE_hs_fourier     <- log(df_last10$MSE_hs_fourier)
df_last10$log_MSE_hs_legendre    <- log(df_last10$MSE_hs_legendre)
df_last10$log_MSE_lasso_bspline  <- log(df_last10$MSE_lasso_bspline)
df_last10$log_MSE_lasso_fourier  <- log(df_last10$MSE_lasso_fourier)
df_last10$log_MSE_lasso_legendre <- log(df_last10$MSE_lasso_legendre)
df_last10$log_log_n <- log(df_last10$log_n)

cat("\n--- Regressions (Last 10): log(MSE) ~ log(n) ---\n")
summary(lm(log_MSE_hs_bspline     ~ log_n, data = df_last10))
summary(lm(log_MSE_hs_fourier     ~ log_n, data = df_last10))
summary(lm(log_MSE_hs_legendre    ~ log_n, data = df_last10))
summary(lm(log_MSE_lasso_bspline  ~ log_n, data = df_last10))
summary(lm(log_MSE_lasso_fourier  ~ log_n, data = df_last10))
summary(lm(log_MSE_lasso_legendre ~ log_n, data = df_last10))


