require(glmnet)
require(R.matlab)

soft_threshhold <- function(innerprod, lambda) {
  return(sign(innerprod)*max(0, abs(innerprod)- lambda))
}


update_beta <- function(y, X, beta_new, lambda, alpha){
  beta_new_2 =c()
  n = length(y)
  a = colSums(X^2)
  p = length(beta_new)
  for(j in 1:p){
    beta_new[j] = soft_threshhold((crossprod(X[,j], (y-X[,-j]%*%beta_new[-j])))/n, lambda*alpha)
    beta_new[j] = beta_new[j]/(a[j]/n+lambda*(1-alpha))
  }
  return(beta_new)
}


Elastic <- function(y, X, beta0, alpha, lambda, niter = 1e+4, tol = 1e-5){
  error = 1000
  n = 1
  all.beta = beta0
  gamma = lambda*alpha
  a = colSums(X^2)
  while (error > tol) {
    beta_new <- c(update_beta(y, X, beta0, lambda, alpha))
    all.beta <- rbind(all.beta, beta_new)
    if (any(abs(beta_new) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(beta_new - beta0)^0.5
    if (n == niter) {
      warning("Max iteration reached")
      break
    }
    beta0 <- beta_new
    n <- n + 1
  }
  return(beta_new)
}
