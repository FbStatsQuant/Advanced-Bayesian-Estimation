library(splines)
library(MLmetrics)
library(Hmisc)
library(horseshoe)


set.seed(2024)

df <- read.csv("CSIRO_Recons_gmsl_yr_2019.csv", header = TRUE, sep = ",")
x <- unlist(df["Time"])
y <- unlist(df["GMSL"])
n <- length(x)
b <- 0.8
J <- floor(n^b)

t <- c(1879.98, 1879.99, seq(1880,2019, length.out = J-4),2019.01, 2019.02)
Sp_order <- 3
B <- splineDesign(t, x, ord = Sp_order, outer.ok = TRUE,
                  sparse = FALSE)

Hs5 <- horseshoe(y, B, method.tau = c("halfCauchy"),
                method.sigma = c("Jeffreys"), burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
beta_hat5 <- as.vector(unlist(Hs5[1]))

MSE(y,as.numeric(B%*%beta_hat5))
plot(x,y, xlab = "Years",
     ylab = " Millimeters", col = "black", cex.main = 2)
points(x,B%*%beta_hat5, type = "l",
      col = "red")
# Add a legend
legend("topleft",
       legend = c("Estimation", "Function f(t)"),
       col = c("red", "black"),
       lty = c(1, 1))



# Compute fitted values
fitted_values <- as.vector(B %*% beta_hat5)

# Create data frame
df_plot <- data.frame(
  x = x,
  y = y,
  fit = fitted_values
)

# Plot using ggplot
ggplot(df_plot, aes(x = x)) +
  geom_point(aes(y = y, color = "Function f(t)"), size = 1, alpha = 0.7) +
  geom_line(aes(y = fit, color = "Estimation"), linewidth = 0.8) +
  scale_color_manual(values = c("Estimation" = "red", "Function f(t)" = "black")) +
  labs(
    title = "CSIRO Global Mean Sea Level Estimation",
    x = "Years",
    y = "Millimeters",
    color = "Legend"
  ) +
  theme_minimal()
