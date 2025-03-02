set.seed(123)

library(splines)
library(MLmetrics)
library(ggplot2)

# Simulation parameters
n <- 1000
Sp_order <- 3
x <- seq(0, 1, length.out = n)
tau_values <- c(0.2, 0.5, 1.0)
J_values <- seq(200, 1200, by = 200)
num_seeds <- 1

# True function and noise parameters
s <- function(m, x) {
  sqrt(2) * sum(seq(m)^(-1.5) * sin(seq(m)) * cos(((seq(m)) - 0.5) * pi * x))
}
z <- vector()
for (j in 1:n) {
  z[j] <- s(1000, j / n)
}

sd <- 0.2
sigma_sq <- sd^2

# Precompute basis matrices for all J values
B_list <- lapply(J_values, function(j) {
  knots <- c(-0.0002, -0.0001, seq(0, 1, length.out = j - 4), 1.0001, 1.0002)
  splineDesign(knots, x, ord = Sp_order, outer.ok = TRUE, sparse = FALSE)
})

# Initialize storage
results <- list()
random_seeds <- sample(0:2000, num_seeds)

# Main simulation loop
for (tau in tau_values) {
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
  
  # Store results with clear labeling
  results[[paste0("tau_", tau)]] <- data.frame(
    J = J_values,
    Avg_MSE = colMeans(MSE_matrix),
    Tau = factor(paste("τ =", tau), 
                 levels = paste("τ =", tau_values))
  )  # Added closing ) for data.frame and factor
}  # Closing bracket for tau loop

# Combine all results
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

num_Hs_seeds <- 1
Hs_seeds <- sample(0:2000, num_Hs_seeds)  # Generate 10 different seeds
Hs_seeds <- 2024
# Initialize storage for Horseshoe results
Hs_MSE_matrix <- matrix(0, nrow = num_Hs_seeds, ncol = length(J_values))

# Loop over seeds and J values
for (s in 1:num_Hs_seeds) {
  set.seed(Hs_seeds[s])  # Set seed for this iteration
  y <- z + rnorm(n, 0, sd)  # Generate new noisy observations
  
  for (j in seq_along(J_values)) {
    J <- J_values[j]  # Get current J
    
    # Construct B-spline basis matrix
    knots <- c(-0.0002, -0.0001, seq(0, 1, length.out = J - 4), 1.0001, 1.0002)
    B <- splineDesign(knots, x, ord = Sp_order, outer.ok = TRUE, sparse = TRUE)
    
    # Run Horseshoe prior regression
    Hs1 <- horseshoesp(y, B, method.tau = c("truncatedCauchy"), Sigma2 = (sd)^2, tau=1, s=0.5,
                       method.sigma = c("fixed"), burn = 1000, nmc = 5000, thin = 1)
    
    
    # Extract estimated coefficients
    beta_hat1 <- as.vector(unlist(Hs1[1])) 
    
    # Compute MSE for Horseshoe model
    Hs_MSE_matrix[s, j] <- MSE(z, as.numeric(B %*% beta_hat1))
  }
}

# Compute average MSE across 10 seeds
Hs_MSE_avg <- colMeans(Hs_MSE_matrix)

# Create data frame for Horseshoe results
hs_results <- data.frame(
  J = J_values,
  Avg_MSE = Hs_MSE_avg,
  Tau = factor("Horseshoe Prior")
)

# Combine all results
plot_data <- rbind(plot_data, hs_results)

# Plot the results
ggplot(plot_data, aes(x = J, y = Avg_MSE, color = Tau)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "MSE vs Number of Knots (J) for Different Priors",
       x = "Number of Knots (J)",
       y = "MSE",
       color = "Prior Type") +
  theme_minimal() +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "black")) +  # Black for HS prior
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))


tauout <- Hs1$tauSamples
hist(tauout)

