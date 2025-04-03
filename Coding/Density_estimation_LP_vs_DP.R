library(dirichletprocess)
library(ggplot2)

data(faithful)  # Bimodal density example
x <- faithful$waiting

# Create Dirichlet Process object with Normal mixture
dp <- DirichletProcessGaussian(x)
dp <- Fit(dp, 1000)  # Fit with 1000 MCMC iterations

# Generate posterior predictive samples
posterior_dens <- data.frame(x = seq(min(x)-5, max(x)+5, length.out=200))
posterior_dens$y <- PosteriorFunction(dp)(posterior_dens$x)

# Plot true data vs. DPM estimate
ggplot() +
  geom_histogram(aes(x = x, y = after_stat(density)), data = data.frame(x), bins = 30, alpha = 0.5) +
  geom_line(aes(x = x, y = y), data = posterior_dens, color = "red", linewidth = 1) +
  labs(title = "Density Estimation with Dirichlet Process Mixture",
       x = "Waiting Time (min)", y = "Density") +
  theme_minimal()

# KDE for comparison
kde_dens <- density(x, bw = "SJ")  # Sheather-Jones bandwidth selector

# Overlay plots
ggplot() +
  geom_histogram(aes(x = x, y = after_stat(density)), data = data.frame(x), bins = 30, alpha = 0.5) +
  geom_line(aes(x = x, y = y), data = posterior_dens, color = "red", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = kde_dens$x, y = kde_dens$y), color = "blue", linewidth = 1, linetype = "dashed") +
  labs(title = "DPM (Red) vs. KDE (Blue) Density Estimates",
       x = "Waiting Time (min)", y = "Density") +
  theme_minimal()