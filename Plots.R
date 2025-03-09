library(readxl)
library(ggplot2)
library(tidyr)

# Load the Excel file
data <- read_excel("MSE_ghosal_3.xlsx")

long_data <- pivot_longer(data, -J, names_to = "Metric", values_to = "MSE")

long_data <- pivot_longer(data, -J, names_to = "Metric", values_to = "MSE")

colors <- c(
  "mse_s_hs" = "#1f77b4",        # Blue
  "mse_s_ridge" = "#ff7f0e",     # Orange
  "mse_f_hs" = "#2ca02c",       # Green (Fourier HS)
  "mse_f_ridge" = "#98df8a",     # Light Green
  "mse_l_hs" = "#d62728",        # Red
  "mse_l_ridge" = "#9467bd"      # Purple
)

labels <- c(
  "mse_s_hs" = "MSE (B-Splines, HS)",
  "mse_s_ridge" = "MSE (B-Splines, Ridge)",
  "mse_f_hs" = "MSE (Fourier, HS)",
  "mse_f_ridge" = "MSE (Fourier, Ridge)",
  "mse_l_hs" = "MSE (Legendre, HS)",
  "mse_l_ridge" = "MSE (Legendre, Ridge)"
)

ggplot(long_data, aes(x = J, y = MSE, color = Metric)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors, labels = labels) +
  labs(x = "J", y = "Mean Squared Error", color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")


data_fourier <- data[, c("J", "mse_f_hs", "mse_f_ridge")]

long_data_fourier <- pivot_longer(data_fourier, -J, names_to = "Metric", values_to = "MSE")

colors_fourier <- c("mse_f_hs" = "#2ca02c", "mse_f_ridge" = "#98df8a")
labels_fourier <- c("mse_f_hs" = "MSE (Fourier, HS)", "mse_f_ridge" = "MSE (Fourier, Ridge)")

ggplot(long_data_fourier, aes(x = J, y = MSE, color = Metric)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors, labels = labels) +
  labs(x = "J", y = "Mean Squared Error", color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")

data_legendre <- data[, c("J", "mse_l_hs", "mse_l_ridge")]

long_data_legendre <- pivot_longer(data_legendre, -J, names_to = "Metric", values_to = "MSE")

colors_legendre <- c("mse_l_hs" = "#1f77b4", "mse_l_ridge" = "#aec7e8")
labels_legendre <- c("mse_l_hs" = "MSE (Legendre, HS)", "mse_l_ridge" = "MSE (Legendre, Ridge)")

ggplot(long_data_legendre, aes(x = J, y = MSE, color = Metric)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_legendre, labels = labels_legendre) +
  labs(x = "J", y = "Mean Squared Error", color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")

data_splines <- data[, c("J", "mse_s_hs", "mse_s_ridge")]

long_data_splines <- pivot_longer(data_splines, -J, names_to = "Metric", values_to = "MSE")

colors_splines <- c("mse_s_hs" = "black", "mse_s_ridge" = "red")
labels_splines <- c("mse_s_hs" = "MSE (B-Splines, HS)", "mse_s_ridge" = "MSE (B-Splines, Ridge)")

ggplot(long_data_splines, aes(x = J, y = MSE, color = Metric)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_splines, labels = labels_splines) +
  labs(x = "J", y = "Mean Squared Error", color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")




