#' Horseshoe Sparse Regression
#'
#' This function performs Horseshoe regression optimized for sparse matrices.
#'
#' @param y Response vector.
#' @param X Sparse design matrix
#' @return A fitted model object.
#' @importFrom Matrix Cholesky solve Diagonal
#' @export

horseshoesp = function(y,X, method.tau = c("fixed", "truncatedCauchy","halfCauchy"), tau = 1,
                       method.sigma = c("fixed", "Jeffreys"), Sigma2 = 1, s=2,
                       burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
{
  method.tau <- match.arg(method.tau, choices = c("fixed", "truncatedCauchy", "halfCauchy"))
  method.sigma <- match.arg(method.sigma, choices = c("fixed", "Jeffreys"))


  X= Matrix::Matrix(X, sparse = TRUE)
  N=burn+nmc
  effsamp=(N-burn)/thin
  n=nrow(X)
  p=ncol(X)

  Beta=rep(0,p); lambda=rep(1,p)
  sigma_sq = Sigma2;

  ## output ##
  betaout=matrix(0,p,effsamp)
  lambdaout=matrix(0,p,effsamp)
  tauout=rep(0,effsamp)
  sigmaSqout=rep(0,effsamp)

  if (p > n) {
    algo <- 1
  } else {
    algo <- 2
  }

  I_n=diag(n)
  l0=rep(0,p)
  l1=rep(1,n)
  l2=rep(1,p)
  if(algo==2)
  {
    Q_star = Matrix::crossprod(X)
  }
  for(i in 1:N){

    #Beta
    lambda_star <- tau * lambda
    if(algo==1){
      U = Matrix::Diagonal(x = as.numeric(lambda_star^2)) %*% Matrix::t(X)
      u=stats::rnorm(l2,l0,lambda_star)
      v=X%*%u + stats::rnorm(n)
      C=Matrix::Cholesky((X%*%U+I_n),
                         perm=FALSE, LDL=FALSE)
      L=Matrix::sparseMatrix(i = as.integer(C@i) + 1,
                             p = as.integer(C@p),
                             x = as.numeric(C@x),
                             dims = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                             index1 = TRUE)
      rhs = matrix((y/sqrt(sigma_sq)) - v, ncol = 1)
      z = Matrix::solve(L, as.numeric(rhs))
      v = Matrix::solve(Matrix::t(L), z)
      Beta=as.numeric(sqrt(sigma_sq)*(u+U%*%v))
    }
    if(algo==2){
      C=Matrix::Cholesky((1 / sigma_sq) * (Q_star + Matrix::Diagonal(x = 1 / (lambda_star^2))),
                         perm=FALSE, LDL=FALSE)
      L=Matrix::sparseMatrix(i = as.integer(C@i) + 1,
                             p = as.integer(C@p),
                             x = as.numeric(C@x),
                             dims = c(as.integer(C@Dim[1]), as.integer(C@Dim[1])),
                             index1 = TRUE)

      rhs=Matrix::crossprod(X, y) / sigma_sq
      z=Matrix::solve(L, as.numeric(rhs))
      mu=Matrix::solve(Matrix::t(L), as.numeric(z))
      u=Matrix::solve(Matrix::t(L),stats::rnorm(p))
      Beta=as.numeric(mu + u)
    }

    #Lambda
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = Beta^2/(2*sigma_sq*tau^2)
    ub = (1-upsi)/upsi
    Fub = 1 - exp(-tempps*ub)
    Fub[Fub < (1e-4)] = 1e-4;
    up = stats::runif(p,0,Fub)
    eta = -log(1-up)/tempps
    lambda = 1/sqrt(eta);

    ##Tau
    if(method.tau == "halfCauchy"){
      tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt = (1-utau)/utau
      Fubt = stats::pgamma(ubt,(p+1)/2,scale=1/tempt)
      Fubt = max(Fubt,1e-8) # for numerical stability
      ut = stats::runif(1,0,Fubt)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = max(1/sqrt(et), 1e-6)
    }

    if(method.tau == "truncatedCauchy"){
      tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt_1=1
      ubt_2 = min((1-utau)/utau,p^(s))
      Fubt_1 = stats::pgamma(ubt_1,(p+1)/2,scale=1/tempt)
      Fubt_2 = stats::pgamma(ubt_2,(p+1)/2,scale=1/tempt)
      ut = stats::runif(1,Fubt_1,Fubt_2)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = max(1/sqrt(et), 1e-6)
    }

    #Sigma
    if(method.sigma == "Jeffreys"){
      if(algo==1)
      {
        E_1=max(t(y-as.numeric(X%*%Beta))%*%(y-as.numeric(X%*%Beta)),(1e-8))
        E_2=max(sum(Beta^2/((tau*lambda))^2),(1e-8))
      }

      else
      {
        E_1=max(t(as.numeric(X%*%Beta))%*%(y-as.numeric(X%*%Beta)),1e-6)
        E_2=max(sum(Beta^2/((tau*lambda))^2),1e-6)
      }
      sigma_sq=1/stats::rgamma(1,(n+p)/2,scale=2/(E_1+E_2))
    } else {
      sigma_sq=Sigma2
    }

    if(i > burn && i%%thin== 0)
    {
      betaout[,(i-burn)/thin] = Beta
      lambdaout[,(i-burn)/thin] = lambda
      tauout[(i-burn)/thin]=tau
      sigmaSqout[(i-burn)/thin]=sigma_sq
    }
  }

  pMean=apply(betaout,1,mean)
  pMedian=apply(betaout,1,stats::median)
  pLambda=apply(lambdaout,1,mean)
  pSigma=mean(sigmaSqout)
  pTau=mean(tauout)

  result=list("BetaHat"=pMean, "BetaMedian"=pMedian,
              "Sigma2Hat"=pSigma,"TauHat"=pTau, "LambdaHat"=pLambda)
  return(result)

}

