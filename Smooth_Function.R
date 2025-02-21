library(splines)
library(MLmetrics)
library(Hmisc)
library(ggplot2)
library(sparsehorseshoe)
library(Matrix)

set.seed(2024)
n <- 1000
b <- 0.5
J <-floor(n^b)
x <- seq(0,1,length.out = n)
f <- function(x){
  exp(-x^2)}
z <- f(x)
sd <- 0.3
y <- z+rnorm(n,mean=0, sd=sd)

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

t <- c(-0.002, -0.001, seq(0,1, length.out = J),1.001, 1.002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = TRUE)

Hs1 <- horseshoesp(y, B, method.tau = c("halfCauchy"), Sigma2 = (sd)^2, tau=1,
                   method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat1 <- as.vector(unlist(Hs1[1]))

MSE(z,as.numeric(B%*%beta_hat1))

df_estimation <- data.frame(x=x, estimated_y=as.numeric(B%*%beta_hat1), actual_y=f(x))

ggplot(df_estimation, aes(x=x))+
  geom_point(aes(y = estimated_y, color = "Estimation"), size = 0.8) +
  geom_line(aes(y = actual_y, color = "Actual Function"), size = 1) +
  labs(x = "x", y = "y", color=NULL) +
  scale_color_manual(values = c("lightblue", "black")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))



J_values <- c(seq(200,2000, by =200))
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