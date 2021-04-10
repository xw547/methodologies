require(doParallel)
require(ggplot2)
require(Matrix)
cl <- makeCluster(8)
registerDoParallel(cl)

###Generate the upper triangle matrix without the diagonal by finding the lower
###half and then delete them from the input matrix.
getutri <- function(A){
  return(as.matrix(A-tril(A)))
}

###Same idea from the previous one.
getltri <- function(A){
  return(as.matrix(A-triu(A)))
}


###Jacobi Method

Jacobi <- function(A, b, x){
  a <- diag(A)
  diag(A) <- 0
  x_new <- c()
  for (i in 1:length(x)) {
    x_new[i] <- (b[i] - crossprod(A[i, ], x))/a[i]
  }
  return(x_new)
}

seq_Jacobi <- function(A, b, x0, maxiter = 10000, tol =1e-5){
  error = 1000
  n = 1
  all.x = x0
  while (error > tol) {
    x <- c(Jacobi(A, b, x0))
    all.x <- rbind(all.x, x)
    if (any(abs(x) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(x - x0)^0.5
    if (n == maxiter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(all.x)
}

###Parallel Jacobi
PJacobi <- function(A, b, x){
  a <- diag(A)
  diag(A) <- 0
  outlist = foreach (i = 1:length(x),  .multicombine = TRUE)  %dopar% {
    (b[i] - crossprod(A[i, ], x))/a[i]
  }
  x_new = as.vector(unlist(outlist))
  return(x_new)
}

para_Jacobi <- function(A, b, x0, maxiter = 10000, tol =1e-5){
  error = 1000
  n = 1
  all.x = x0
  while (error > tol) {
    x <- c(PJacobi(A, b, x0))
    all.x <- rbind(all.x, x)
    if (any(abs(x) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(x - x0)^0.5
    if (n == maxiter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(all.x)
}