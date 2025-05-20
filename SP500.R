library(readxl)
library(horseshoe)
library(MLmetrics)
library(Hmisc)
library(forecast)
library(ggplot2)

set.seed(2024)
df <- read.csv("SP500.csv", header = TRUE, sep = ",")
y <- rev(as.vector(unlist(df['Close'])))

n <- length(y)
x_orig <- 1:n                       # original time index
x_scaled <- x_orig / n             # rescaled to (0,1)
b <- 0.8
J <- floor(n^b)

# Fourier basis using scaled x
fourier_basis <- function(x, K) {
  basis <- matrix(1, nrow = length(x), ncol = 2 * K + 1)
  for (k in 1:K) {
    basis[, 2 * k]     <- sin(2 * pi * k * x)
    basis[, 2 * k + 1] <- cos(2 * pi * k * x)
  }
  basis
}

B_f <- fourier_basis(x_scaled, floor(J / 2))

# Horseshoe estimation
Hs_f <- horseshoe(y, B_f, method.tau = "halfCauchy", method.sigma = "Jeffreys")
theta_f_hs <- unlist(Hs_f$BetaHat)
f_hat_f_hs <- B_f %*% theta_f_hs

# Exclude first and last 10 points
cut <- 10
data <- data.frame(
  Time = x_orig[(cut+1):(n - cut)],
  True = y[(cut+1):(n - cut)],
  Estimated = as.vector(f_hat_f_hs)[(cut+1):(n - cut)]
)

# Subset every 2nd point for estimated values
scatter_data <- data[seq(1, nrow(data), by = 10), ]

ggplot() +
  geom_line(data = data, aes(x = Time, y = True, color = "True (Observed)"), size = 0.8) +
  geom_point(data = scatter_data, aes(x = Time, y = Estimated, color = "Estimated (Horseshoe + Fourier)"),
             size = 1.1, alpha = 0.7) +
  scale_color_manual(values = c("True (Observed)" = "lightblue",
                                "Estimated (Horseshoe + Fourier)" = "blue")) +
  labs(
    title = "True vs Estimated S&P 500 (Horseshoe + Fourier)",
    x = "Time Index",
    y = "Closing Price",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = c(1, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(hjust = 0.5)
  )
