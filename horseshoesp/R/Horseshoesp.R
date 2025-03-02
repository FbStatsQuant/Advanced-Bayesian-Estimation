horseshoesp <- function(y, X,
                        method.tau = c("fixed", "truncatedCauchy", "halfCauchy"),
                        tau = 1,
                        method.sigma = c("fixed", "Jeffreys"),
                        Sigma2 = 1,
                        s = 2,
                        burn = 1000,
                        nmc = 5000,
                        thin = 1,
                        # extra numerical safeguards
                        ridge = 1e-12,     # small ridge added to Cholesky
                        eps_gamma = 1e-12, # floors for scale in gamma calls
                        eps_F = 1e-12      # floors for Fub, etc.
)
{
  ## Match arguments
  method.tau   <- match.arg(method.tau, choices = c("fixed", "truncatedCauchy", "halfCauchy"))
  method.sigma <- match.arg(method.sigma, choices = c("fixed", "Jeffreys"))

  ## Ensure X is sparse
  X <- Matrix::Matrix(X, sparse = TRUE)

  N       <- burn + nmc
  effsamp <- (N - burn) / thin

  n <- nrow(X)
  p <- ncol(X)

  ## Initialize parameters
  Beta      <- rep(0, p)
  lambda    <- rep(1, p)
  sigma_sq  <- Sigma2

  ## Storage for posterior samples
  betaout    <- matrix(0, nrow = p, ncol = effsamp)
  tauout     <- numeric(effsamp)
  sigmaSqout <- numeric(effsamp)

  ## Decide algorithm based on dimensions
  algo <- if (p > n) 1 else 2

  ## Precompute Q_star if p <= n
  if (algo == 2) {
    Q_star <- Matrix::crossprod(X)  # p x p
  }

  ## Identity
  I_n <- Matrix::Diagonal(n)

  for (i in seq_len(N)) {

    ##-----------------------------------------------------------
    ## 1) Sample Beta
    ##-----------------------------------------------------------
    lambda_star <- tau * lambda  # vector length p

    if (algo == 1) {
      # p > n approach
      U <- Matrix::Diagonal(x = as.numeric(lambda_star^2)) %*% Matrix::t(X)
      u <- stats::rnorm(p, mean = 0, sd = lambda_star)
      v <- as.numeric(X %*% u) + stats::rnorm(n)

      CU <- X %*% U + I_n + ridge * I_n
      C  <- Matrix::Cholesky(CU, perm = FALSE, LDL = FALSE)

      L <- Matrix::sparseMatrix(i     = as.integer(C@i) + 1,
                                p     = as.integer(C@p),
                                x     = as.numeric(C@x),
                                dims  = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                                index1= TRUE)

      rhs <- (y / sqrt(sigma_sq)) - v
      z   <- Matrix::solve(L, rhs)
      vv  <- Matrix::solve(Matrix::t(L), z)

      Beta <- as.numeric(sqrt(sigma_sq) * (u + U %*% vv))

    } else {
      # p <= n approach
      M <- (1/sigma_sq) * (Q_star + Matrix::Diagonal(x = 1/(lambda_star^2))) +
        ridge * Matrix::Diagonal(p)
      C <- Matrix::Cholesky(M, perm = FALSE, LDL = FALSE)

      L <- Matrix::sparseMatrix(i     = as.integer(C@i) + 1,
                                p     = as.integer(C@p),
                                x     = as.numeric(C@x),
                                dims  = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                                index1= TRUE)

      rhs <- Matrix::crossprod(X, y) / sigma_sq
      z   <- Matrix::solve(L, as.numeric(rhs))
      mu  <- Matrix::solve(Matrix::t(L), z)

      u    <- Matrix::solve(Matrix::t(L), stats::rnorm(p))
      Beta <- as.numeric(mu + u)
    }

    ##-----------------------------------------------------------
    ## 2) Sample local scales lambda
    ##-----------------------------------------------------------
    eta  <- 1 / (lambda^2)
    upsi <- stats::runif(p, 0, 1/(1 + eta))

    tempps <- Beta^2 / (2 * sigma_sq * tau^2)
    ub   <- (1 - upsi) / upsi
    Fub  <- 1 - exp(-tempps * ub)

    Fub[Fub < eps_F] <- eps_F
    Fub[Fub > 1 - eps_F] <- 1 - eps_F

    up <- stats::runif(p, 0, Fub)
    one_minus_up <- 1 - up
    one_minus_up[one_minus_up < eps_F] <- eps_F

    eta    <- -log(one_minus_up) / tempps
    lambda <- 1 / sqrt(eta)

    ##-----------------------------------------------------------
    ## 3) Sample global scale tau (truncated Cauchy update FIXED)
    ##-----------------------------------------------------------
    if (method.tau == "truncatedCauchy") {
      tempt <- sum((Beta / lambda)^2) / (2 * sigma_sq)
      tempt <- max(tempt, eps_gamma)

      et    <- 1 / tau^2
      utau  <- stats::runif(1, 0, 1/(1 + et))

      # **FIXED**: Ensure tau ≤ 1
      ubt_1 <- 1
      ubt_2 <- min((1 - utau) / utau, 1)  # **Upper bound set to 1**

      shape_par <- (p + 1)/2
      Fubt_1 <- stats::pgamma(ubt_1, shape = shape_par, scale = 1/tempt)
      Fubt_2 <- stats::pgamma(ubt_2, shape = shape_par, scale = 1/tempt)

      # Clamp to avoid NA values
      if (is.na(Fubt_1)) Fubt_1 <- eps_F
      if (is.na(Fubt_2)) Fubt_2 <- 1 - eps_F

      Fubt_1 <- min(max(Fubt_1, eps_F), 1 - eps_F)
      Fubt_2 <- min(max(Fubt_2, eps_F), 1 - eps_F)

      if (Fubt_1 > Fubt_2) {
        tmp <- Fubt_1
        Fubt_1 <- Fubt_2
        Fubt_2 <- tmp
      }

      ut <- stats::runif(1, Fubt_1, Fubt_2)
      et <- stats::qgamma(ut, shape = shape_par, scale = 1/tempt)

      tau <- min(max(1 / sqrt(et), 1e-8), 1)  # **Clamp tau to ≤ 1**
    }

    ##-----------------------------------------------------------
    ## 4) Sample sigma^2
    ##-----------------------------------------------------------
    if (method.sigma == "Jeffreys") {
      resid <- y - as.numeric(X %*% Beta)
      E_1   <- max(sum(resid^2), 1e-8)
      E_2   <- max(sum(Beta^2 / (tau * lambda)^2), 1e-8)
      sigma_sq <- 1 / stats::rgamma(1, shape = (n + p)/2, rate = (E_1 + E_2)/2)
    } else {
      sigma_sq <- Sigma2
    }

    ##-----------------------------------------------------------
    ## Store samples after burn-in
    ##-----------------------------------------------------------
    if (i > burn && (i - burn) %% thin == 0) {
      idx <- (i - burn) / thin
      betaout[, idx]    <- Beta
      tauout[idx]       <- tau
      sigmaSqout[idx]   <- sigma_sq
    }
  }

  ##-----------------------------------------------------------
  ## Posterior summaries
  ##-----------------------------------------------------------
  pMean   <- rowMeans(betaout)
  pMedian <- apply(betaout, 1, stats::median)
  pSigma  <- mean(sigmaSqout)
  pTau    <- mean(tauout)

  result <- list(
    BetaHat     = pMean,
    BetaMedian  = pMedian,
    Sigma2Hat   = pSigma,
    TauHat      = pTau,
    betaSamples = betaout,
    tauSamples  = tauout,
    sigma2Samps = sigmaSqout
  )

  return(result)
}
