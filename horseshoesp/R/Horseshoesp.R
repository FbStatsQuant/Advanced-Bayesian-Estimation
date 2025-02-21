#' Horseshoe Sparse Regression
#'
#' This function performs Horseshoe regression optimized for sparse matrices.
#'
#' @param y Response vector.
#' @param X Sparse design matrix (class "dgCMatrix" or similar).
#' @param method.tau Method for sampling global scale tau: "fixed", "truncatedCauchy", or "halfCauchy".
#' @param tau Fixed value of tau (used if method.tau = "fixed").
#' @param method.sigma Method for sigma^2: "fixed" or "Jeffreys".
#' @param Sigma2 Fixed value of sigma^2 (used if method.sigma = "fixed").
#' @param s Exponent used in truncated Cauchy sampling of tau.
#' @param burn Number of burn-in iterations.
#' @param nmc Number of post-burn-in MCMC draws.
#' @param thin Thinning interval.
#'
#' @return A list with posterior means and medians of coefficients, plus mean sigma^2 and tau.
#' @importFrom Matrix Cholesky solve Diagonal crossprod
#' @export
horseshoesp <- function(y, X,
                        method.tau = c("fixed", "truncatedCauchy","halfCauchy"),
                        tau = 1,
                        method.sigma = c("fixed", "Jeffreys"),
                        Sigma2 = 1,
                        s = 2,
                        burn = 1000,
                        nmc = 5000,
                        thin = 1)
{
  ## Match arguments
  method.tau   <- match.arg(method.tau,   choices = c("fixed", "truncatedCauchy", "halfCauchy"))
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
  #  algo=1 if p>n, else algo=2
  algo <- if (p > n) 1 else 2

  ## Precompute Q_star if p <= n
  if (algo == 2) {
    Q_star <- Matrix::crossprod(X)  # p x p
  }

  ## Identity and placeholders
  I_n <- Matrix::Diagonal(n)

  for (i in seq_len(N)) {

    ## 1) Sample Beta
    lambda_star <- tau * lambda  # vector length p

    if (algo == 1) {
      # p > n approach
      # U is p x n
      U <- Matrix::Diagonal(x = as.numeric(lambda_star^2)) %*% Matrix::t(X)

      # Draw u ~ N(0, diag(lambda_star^2))
      # rnorm(p, mean=0, sd=lambda_star) is vectorized
      u <- stats::rnorm(p, mean = 0, sd = lambda_star)

      # v = X u + e, where e ~ N(0, I_n)
      v <- as.numeric(X %*% u) + stats::rnorm(n)

      # Construct and factor (X U + I_n)
      CU <- X %*% U + I_n  # n x n
      C  <- Matrix::Cholesky(CU, perm = FALSE, LDL = FALSE)

      # Convert the Cholesky factor to a sparse matrix "L"
      L <- Matrix::sparseMatrix(i     = as.integer(C@i) + 1,
                                p     = as.integer(C@p),
                                x     = as.numeric(C@x),
                                dims  = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                                index1= TRUE)

      # Solve for v in the system L t(L) v = rhs
      rhs <- (y / sqrt(sigma_sq)) - v
      z   <- Matrix::solve(L, rhs)
      vv  <- Matrix::solve(Matrix::t(L), z)

      # Beta = sqrt(sigma_sq)*(u + U v)
      Beta <- as.numeric( sqrt(sigma_sq) * (u + U %*% vv) )

    } else {
      # p <= n approach
      # Factor (1/sigma_sq)*(Q_star + diag(1/lambda_star^2))
      C <- Matrix::Cholesky((1/sigma_sq) * (Q_star +
                                              Matrix::Diagonal(x = 1/(lambda_star^2))),
                            perm = FALSE, LDL = FALSE)

      L <- Matrix::sparseMatrix(i     = as.integer(C@i) + 1,
                                p     = as.integer(C@p),
                                x     = as.numeric(C@x),
                                dims  = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                                index1= TRUE)

      rhs <- Matrix::crossprod(X, y) / sigma_sq  # p-vector
      z   <- Matrix::solve(L, as.numeric(rhs))
      mu  <- Matrix::solve(Matrix::t(L), z)

      # Add normal draw
      u    <- Matrix::solve(Matrix::t(L), stats::rnorm(p))
      Beta <- as.numeric(mu + u)
    }

    ## 2) Sample local scales lambda
    #    Polya-Gamma type update in HS can be done via slice sampler or this approach
    eta   <- 1 / (lambda^2)
    upsi  <- stats::runif(p, 0, 1/(1 + eta))  # uniform(0, 1/(1+eta))

    tempps <- Beta^2 / (2 * sigma_sq * tau^2)
    ub     <- (1 - upsi) / upsi
    Fub    <- 1 - exp(-tempps * ub)

    # Floor Fub away from zero
    Fub[Fub < 1e-8] <- 1e-8

    up   <- stats::runif(p, 0, Fub)
    eta  <- -log(1 - up) / tempps
    lambda <- 1 / sqrt(eta)

    ## 3) Sample global scale tau
    if (method.tau == "halfCauchy") {
      tempt <- sum((Beta / lambda)^2) / (2 * sigma_sq)
      et    <- 1 / tau^2

      utau  <- stats::runif(1, 0, 1/(1 + et))
      ubt   <- (1 - utau) / utau
      Fubt  <- stats::pgamma(ubt, shape = (p + 1)/2, scale = 1/tempt)
      Fubt  <- max(Fubt, 1e-8)

      ut    <- stats::runif(1, 0, Fubt)
      et    <- stats::qgamma(ut, shape = (p + 1)/2, scale = 1/tempt)
      tau   <- max(1 / sqrt(et), 1e-8)

    } else if (method.tau == "truncatedCauchy") {
      tempt <- sum((Beta / lambda)^2) / (2 * sigma_sq)
      et    <- 1 / tau^2

      utau  <- stats::runif(1, 0, 1/(1 + et))
      ubt_1 <- 1
      ubt_2 <- min( (1 - utau)/utau, p^s )

      Fubt_1 <- stats::pgamma(ubt_1, shape = (p + 1)/2, scale = 1/tempt)
      Fubt_2 <- stats::pgamma(ubt_2, shape = (p + 1)/2, scale = 1/tempt)

      # Draw from [Fubt_1, Fubt_2]
      ut     <- stats::runif(1, Fubt_1, Fubt_2)
      et     <- stats::qgamma(ut, shape = (p + 1)/2, scale = 1/tempt)
      tau    <- max(1 / sqrt(et), 1e-8)

    }
    # if method.tau == "fixed", do nothing with tau

    ## 4) Sample sigma^2
    if (method.sigma == "Jeffreys") {
      # Residual sum of squares
      resid <- y - as.numeric(X %*% Beta)
      E_1   <- max(sum(resid^2), 1e-8)

      # Penalty from Beta part
      E_2   <- max(sum(Beta^2 / (tau * lambda)^2), 1e-8)

      sigma_sq <- 1 / stats::rgamma(1, shape = (n + p)/2, rate = (E_1 + E_2)/2)
    } else {
      # method.sigma == "fixed"
      sigma_sq <- Sigma2
    }

    ## Store samples after burn-in
    if (i > burn && (i - burn) %% thin == 0) {
      idx <- (i - burn) / thin
      betaout[, idx]    <- Beta
      tauout[idx]       <- tau
      sigmaSqout[idx]   <- sigma_sq
    }
  }

  ## Posterior summaries
  pMean   <- rowMeans(betaout)
  pMedian <- apply(betaout, 1, stats::median)
  pSigma  <- mean(sigmaSqout)
  pTau    <- mean(tauout)

  result <- list(
    BetaHat     = pMean,
    BetaMedian  = pMedian,
    Sigma2Hat   = pSigma,
    TauHat      = pTau,
    betaSamples = betaout,    # if you want to inspect full posterior draws
    tauSamples  = tauout,
    sigma2Samps = sigmaSqout
  )

  return(result)
}
