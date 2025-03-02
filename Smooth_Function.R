library(ggplot2)
library(splines)
library(MLmetrics)
library(Hmisc)
library(Matrix)
library(tidyr)
library(dplyr)

n <- 1000
b <- (1/3)
set.seed(2024)
J <- floor(2*n^b)


x <- seq(0, 1, length.out = n)
f <- function(x) exp(-2 * x^2)
z <- f(x)
sd <- 0.15
sigma_sq <- sd^2
(SNR <- var(z)/sd^2)
y <- z + rnorm(n, mean = 0, sd = sd)

lambda <- c(0.1, 1, 10)

(sigma_theta <- sd / sqrt(lambda)) 

data <- data.frame(x = x, y = y, z = z)

ggplot(data, aes(x = x)) +
  geom_line(aes(y = z), color = "blue") +
  labs(x = "x", y = "y")

ggplot(data, aes(x = x)) +
  geom_point(aes(y = y), color = "lightblue") +
  geom_line(aes(y = z), color = "black") +
  labs(title = "Function vs Random Points", x = "x", y = "y") +
  theme(plot.title = element_text(size = 16)) +
  scale_color_manual(values = c("lightblue", "black")) +
  theme_minimal()




t <- c(-0.0002, -0.0001, seq(0,1, length.out = J),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs1 <- horseshoesp(y, B, method.tau = c("halfCauchy"), Sigma2 = (sd)^2, tau=1, s=2,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
beta_hat1 <- as.vector(unlist(Hs1[1]))


# Create a data frame
plot_data <- data.frame(
  x = x,
  True_Function = z,
  Estimated = as.numeric(B%*%beta_hat1)
)

# Plot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = True_Function), color = "red", linewidth = 0.8) +  # True function (line)
  geom_point(aes(y = Estimated), color = "blue", size = 0.5, alpha = 0.5) +  # Estimates (points)
  labs(
    title = "True Function vs. Posterior Estimates",
    x = "x",
    y = "Value"
  ) +
  theme_minimal()

MSE(z,as.numeric(B%*%beta_hat1))




num_seeds <- 1
random_seeds <- sample(0:2000, num_seeds)
J_values <- seq(10, 100, by = 10)

B_list <- lapply(J_values, function(j) {
  knots <- c(-0.0002, -0.0001, seq(0, 1, length.out = j - 4), 1.0001, 1.0002)
  splineDesign(knots, x, ord = Sp_order, outer.ok = TRUE, sparse = FALSE)
})

results <- list()


for (tau in sigma_theta) {
  tau_sq <- tau  # Using τ directly for clarity
  MSE_matrix <- matrix(0, nrow = num_seeds, ncol = length(J_values))
  
  for (i in 1:num_seeds) {
    set.seed(random_seeds[i])
    y <- z + rnorm(n, 0, sd)
    
    for (j in seq_along(J_values)) {
      B <- B_list[[j]]
      p <- ncol(B)
      
      # Bayesian linear regression computation
      precision_matrix <- (1/tau_sq) * diag(p) + (1/sigma_sq) * crossprod(B)
      theta <- solve(precision_matrix, (1/sigma_sq) * crossprod(B, y))
      
      MSE_matrix[i, j] <- MSE(z, B %*% theta)
    }
  }
  results[[paste0("tau_", tau)]] <- data.frame(
    J = J_values,
    Avg_MSE = colMeans(MSE_matrix),
    Tau = factor(paste("τ =", tau), 
                 levels = paste("τ =", sigma_theta))
  )
}

plot_data <- do.call(rbind, results)

# Create visualization
ggplot(plot_data, aes(x = J, y = Avg_MSE, color = Tau)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Average MSE vs Number of Knots for Different Prior Variances",
       x = "Number of Knots (J)",
       y = "Average MSE",
       color = "Prior Variance (τ)") +
  theme_minimal() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))


sigma_theta <- sigma_theta[1]


t <- c(-0.0002, -0.0001, seq(0,1, length.out = J),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs1 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = (sd)^2, tau=0.5, s=0.5,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
beta_hat1 <- as.vector(unlist(Hs1[1]))




df_estimation <- data.frame(x=x, estimated_y=as.numeric(B%*%beta_hat1), actual_y=z)

ggplot(df_estimation, aes(x=x))+
  geom_point(aes(y = estimated_y, color = "Estimation"), size = 0.8) +
  geom_line(aes(y = actual_y, color = "Actual Function"), size = 1) +
  labs(x = "x", y = "y", color=NULL) +
  scale_color_manual(values = c("lightblue", "black")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

precision_matrix <- (1/(sigma_theta)^2) * diag(ncol(B)) + (1/sigma_sq) * crossprod(B)
theta <- solve(precision_matrix, (1/sigma_sq) * crossprod(B, y))


Hs1$TauHat
MSE(z,as.numeric(B%*%beta_hat1))
MSE(z, as.numeric(B%*%theta))





J_values <- c(seq(200,1000, by =200))
mse_values_3_1 <- numeric(length(J_values))
mse_values_3_2 <- numeric(length(J_values))
mse_values_3_3 <- numeric(length(J_values))

for (i in seq_along(J_values)) {
  t <- c(-0.0002, -0.0001, seq(0,1, length.out = J_values[i]-4),1.0001, 1.0002)
  Sp_order <- 3
  B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                    sparse = TRUE)
  Hs1 <- horseshoesp(y, B,
                     method.tau = "halfCauchy",
                     Sigma2 = (sd)^2,
                     tau = 1,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  Hs2 <- horseshoesp(y, B,
                     method.tau = "truncatedCauchy",
                     Sigma2 = (sd)^2,
                     tau = 1,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  Hs3 <- horseshoesp(y, B,
                     method.tau = "truncatedCauchy",
                     Sigma2 = (sd)^2,
                     tau = 1,
                     s = 1,
                     method.sigma = "fixed",
                     burn = 1000,
                     nmc = 5000)
  beta_hat1 <- as.vector(unlist(Hs1[1]))
  beta_hat2 <- as.vector(unlist(Hs2[1]))
  beta_hat3 <- as.vector(unlist(Hs3[1]))
  mse_values_3_1[i] <- MSE(z, as.numeric(B %*% beta_hat1))
  mse_values_3_2[i] <- MSE(z, as.numeric(B %*% beta_hat2))
  mse_values_3_3[i] <- MSE(z, as.numeric(B %*% beta_hat3))
  cat(sprintf("J=%d | MSE_1=%.4f, | MSE_2=%.4f, | MSE_3=%.4f\n", J_values[i], 
              mse_values_3_1[i], mse_values_3_2[i], mse_values_3_3[i]))
}

df <- data.frame(J_values, mse_values_3_1, mse_values_3_2, mse_values_3_3) %>%
  pivot_longer(cols = -J_values, names_to = "MSE_Type", values_to = "MSE")

# Rename MSE_Type for better legend display
df$MSE_Type <- factor(df$MSE_Type, labels = c("MSE 3_1", "MSE 3_2", "MSE 3_3"))

# Plot using ggplot
ggplot(df, aes(x = J_values, y = MSE, color = MSE_Type, shape = MSE_Type)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "MSE vs Number of Knots",
       x = "Number of Knots (J)", 
       y = "MSE") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green")) +
  theme(legend.title = element_blank())