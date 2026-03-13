# ============================================================
#  F1_chains.R
#  Trial run: n = 100, 200, 300, 400, 500
#  Stores per n: x, y, BetaSamples (post burn-in), J, n
#  Output: single .rds saved to /content/
# ============================================================

library(horseshoe)

# ============================================================
#  Settings
# ============================================================
set.seed(1991)

sd_val <- 0.07
b      <- 0.8
f <- function(x) 1.5*(abs(x-0.1))^3 - 5*(abs(x-0.4))^3

fourier_basis <- function(x, K) {
  basis <- matrix(1, nrow = length(x), ncol = K + 1)
  for (k in 1:K) basis[, k + 1] <- cos(pi * k * x)
  basis
}

n_seq   <- seq(1000, 10000, by = 400)
results <- vector("list", length(n_seq))
names(results) <- paste0("n_", n_seq)

# ============================================================
#  Main loop
# ============================================================
for (i in seq_along(n_seq)) {
  
  n <- n_seq[i]
  J <- floor(n^b)
  x <- seq(0, 1, length.out = n)
  y <- f(x) + rnorm(n, mean = 0, sd = sd_val)
  B <- fourier_basis(x, J)
  
  cat(sprintf("\n[%d/%d] n = %d | J = %d\n", i, length(n_seq), n, J))
  
  Hs <- horseshoe(y, B,
                  method.tau   = "fixed", tau = 1 / J,
                  method.sigma = "fixed", Sigma2 = sd_val^2,
                  nmc = 4000, burn = 1000)
  
  # BetaSamples is already post burn-in (nmc rows x J cols, transposed)
  results[[i]] <- list(
    n           = n,
    J           = J,
    x           = x,
    y           = y,
    BetaSamples = Hs$BetaSamples   # dimensions: J x nmc
  )
  
  cat(sprintf("  Done. BetaSamples dim: %d x %d\n",
              nrow(Hs$BetaSamples), ncol(Hs$BetaSamples)))
}

# ============================================================
#  Save to Google Drive
# ============================================================
out_path <- "/content/F2_chains_trial.rds"
saveRDS(results, file = out_path)
cat(sprintf("\nSaved to: %s\n", out_path))
