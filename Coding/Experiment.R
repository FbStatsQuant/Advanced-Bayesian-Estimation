library(ggplot2)
library(splines)

# Data generation
n <- 1000
b <- 1
set.seed(2024)
J <- floor(n^b)
x <- seq(0, 1, length.out = n)
t <- c(-0.0002, -0.0001, seq(0, 1, length.out = J - 4), 1.0001, 1.0002)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)

# Sparse random coefficients
theta <- rep(0, J - Sp_order)
non_zero_indices <- sample(1:(J - Sp_order), size = floor(0.1 * (J - Sp_order)))
theta[non_zero_indices] <- rnorm(length(non_zero_indices), mean = 1, sd = 3)

# Generate z and y
z <- B %*% theta
y <- as.vector(z) + rnorm(n)  # Adding random noise

# Data frame for plotting
df <- data.frame(x = x, y = y, z = as.vector(z))

# Plotting: Scatter plot for (x, y) and line plot for (x, z)
ggplot(df, aes(x = x)) +
  geom_point(aes(y = y), color = "blue", alpha = 0.5, size = 1.5) +  # Scatter plot for (x, y)
  geom_line(aes(y = z), color = "red", size = 1) +                  # Line plot for (x, z)
  labs(title = "Scatter Plot of (x, y) with Line Plot of (x, z)",
       x = "x",
       y = "Values") +
  theme_minimal()


Hs <- horseshoesp(y, B, method.tau = "truncatedCauchy", Sigma2 = 1, tau = 1, k = max(0.5,b/2),
                  method.sigma = "Jeffreys", burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)

theta_hat <- as.vector(unlist(Hs[1]))[non_zero_indices]
theta_true <- theta[c(non_zero_indices)]                      # True theta values
indices <- 1:length(theta_true)          # Index for plotting

# Create data frame for plotting
df_theta <- data.frame(Index = indices, 
                       True_Theta = theta_true, 
                       Estimated_Theta = theta_hat)

ggplot(df_theta, aes(x = Index)) +
  geom_point(aes(y = True_Theta, color = "True θ"), size = 2, alpha = 0.8) +          # True θ scatter points
  geom_point(aes(y = Estimated_Theta, color = "Estimated θ̂"), size = 2, shape = 17, alpha = 0.8) +  # Estimated θ̂ scatter points (different shape)
  labs(title = expression("True " ~ theta ~ " vs Estimated " ~ hat(theta)),
       x = "Index",
       y = expression(theta)) +
  scale_color_manual(values = c("True θ" = "blue", "Estimated θ̂" = "red")) +  # Custom colors
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "top")
