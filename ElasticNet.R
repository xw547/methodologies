require(glmnet)
require(R.matlab)
fil = readMat('lasso.mat')
X = fil$X
y = fil$y

soft_threshhold <- function(innerprod, lambda) {
  if (lambda >= abs(innerprod)){
    result = 0
  } 
  if ((lambda < abs(innerprod))&&(innerprod > lambda)){
    result = innerprod - lambda
  } 
  if ((lambda < abs(innerprod))&&(innerprod <  -lambda)){
    result = innerprod + lambda
  }
  return(result)
}

get_inner_prod <- function(y, X, beta, j){
  n = length(y)
  inner_xj_rj = sum( X[,j] * (y - X[,-j] %*% beta[-j]) )
  return(inner_xj_rj)
}

update_beta <- function(y, X, beta, lambda, a){
  p = length(beta)
  beta_new =c()
  for(j in 1:p){
  beta_new[j] = soft_threshhold(get_inner_prod(y, X, beta, j)/a[j], length(y)*lambda/a[j])
  }
  return(beta_new)
}

LASSO_Coordinate <- function(y, X, beta0, lambda, niter = 1e+4, tol = 1e-5){
  error = 1000
  n = 1
  all.beta = beta0
  a = colSums(X^2)
  while (error > tol) {
    beta <- c(update_beta(y, X, beta0, lambda, a))
    all.beta <- rbind(all.beta, beta)
    if (any(abs(beta) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(beta - beta0)^0.5
    if (n == niter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    beta0 <- beta
  }
  return(all.beta)
}


Elastic <- function(y, X, beta0, alpha, lambda, niter = 1e+4, tol = 1e-5){
  error = 1000
  n = 1
  all.beta = beta0
  gamma = lambda*alpha
  a = colSums(X^2)
  while (error > tol) {
    beta <- c(update_beta(y, X, beta0, gamma, a))
    for (j in 1: length(beta0)){
      beta[j] = beta[j]/(1+lambda*(1-alpha)*length(y)/a[j])
    }
    all.beta <- rbind(all.beta, beta)
    if (any(abs(beta) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(beta - beta0)^0.5
    if (n == niter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    beta0 <- beta
  }
  return(all.beta)
}

bet = c(lm(y~X)$coeff)[2:6]
for (i in 1:10000){
  bet = update_beta(y, X, bet, 0, a)
  }
