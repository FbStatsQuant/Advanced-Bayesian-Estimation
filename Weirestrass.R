library(splines)
library(MLmetrics)
library(Hmisc)
library(horseshoe)
library(Matrix)
library(tidyr)
library(dplyr)

set.seed(2024)
n <- 1000
x <- seq(0,1,length.out = n)
a1 <- 0.8
a2 <- 8
a1*a2 > 1 + 1.5*pi
alpha <- -log(a1)/log(a2)
b <- 1/((2*alpha)+1)
J <- floor(1*(n^b))
c(alpha,b,J)
s <- function(m,x) {
  sum(a1^(seq(m))*cos(a2^seq(10)*pi*x))
}
z <- vector()
for (j in 1:n) {
  z[j] = s(1000,j/n)
}
sd <- 1
y <- z+rnorm(n,mean=0, sd=sd)

data <- data.frame(x = x, y = y, z = z)

# Subsample data
data_subsampled <- data[seq(1, nrow(data), by = 3), ]

# Scatter plot of y-values with subsampled points
ggplot(data_subsampled, aes(x = x)) +
  geom_point(aes(y = y), color = "black", size = 0.8) +  # Random points
  geom_line(aes(y = z), color = "lightblue", size = 0.6) +  # Function line
  labs(title = "Function vs Random Points", x = "x", y = "y") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16))

t <- c(-0.0002, -0.0001, seq(0,1, length.out = (J+4)-4),1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs1 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = (sd)^2, tau=1, s=1,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat1 <- as.vector(unlist(Hs1[1]))

MSE(z,as.numeric(B%*%beta_hat1))

# Create a data frame for estimation vs. actual values
df_estimation <- data.frame(
  x = x,
  estimated_y = as.numeric(B %*% beta_hat1),
  actual_y = z
)

# Subsample data: Keep every 2nd point
df_estimation_subsampled <- df_estimation[seq(1, nrow(df_estimation), by = 3), ]

# Plot estimated vs. actual function with legend in the top right
ggplot(df_estimation_subsampled, aes(x = x)) +
  geom_point(aes(y = estimated_y, color = "Estimation"), size = 1.2) +  # Estimated points
  geom_line(aes(y = actual_y, color = "Actual Function"), size = 0.5) +  # Actual function line
  labs(x = "x", y = "y", color = NULL) +
  scale_color_manual(values = c("Estimation" = "black", "Actual Function" = "lightblue")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.85, 0.85),  # Move legend to top right
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)  # Use linewidth instead of size
  )

J_values <- c(seq(200,2000, by =200))
mse_values_1_1 <- numeric(length(J_values))
mse_values_1_2 <- numeric(length(J_values))
mse_values_1_3 <- numeric(length(J_values))

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
  mse_values_1_1[i] <- MSE(z, as.numeric(B %*% beta_hat1))
  mse_values_1_2[i] <- MSE(z, as.numeric(B %*% beta_hat2))
  mse_values_1_3[i] <- MSE(z, as.numeric(B %*% beta_hat3))
  cat(sprintf("J=%d | MSE_1=%.4f, | MSE_2=%.4f, | MSE_3=%.4f\n", J_values[i], 
              mse_values_1_1[i], mse_values_1_2[i], mse_values_1_3[i]))
}

df <- data.frame(J_values, mse_values_1_1, mse_values_1_2, mse_values_1_3) %>%
  pivot_longer(cols = -J_values, names_to = "MSE_Type", values_to = "MSE")

# Rename MSE_Type for better legend display
df$MSE_Type <- factor(df$MSE_Type, labels = c("MSE 1_1", "MSE 1_2", "MSE 1_3"))

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