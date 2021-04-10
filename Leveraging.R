###Load necessary packages
require(matrixStats)
require(ggplot2)

###Generate Data
n = 500

#x_i \sim t_6
X = rt(500, 6)

#epsilon \sim normal(0,1)
epsilon = rnorm(500, 0, 1)
Y = -X + epsilon

###Basic parameter setuup
draws = 500
Pi_unif=rep(1/n,n)
XTXinv=solve(t(X)%*%X)
H=X%*%XTXinv%*%t(X)
H_ii=diag(H)/sum(diag(H))
Pi_weighted=H_ii

###unif sampling strategy
unif_algo <- function(r,Pi) {
  Beta_Hat_unif_500 = c()
  for (i in 1:draws){
    index_1=sample(c(1:n), size = r, replace = T,prob=Pi_unif) 
    #Phi_1=diag(Pi_unif[index_1])
    X_star_1=X[index_1]
    Y_star_1=Y[index_1]
    Phi_1=1/ (Pi_unif[index_1]^{1/2})
    model_1 <- lm(Y_star_1 ~ 0 + X_star_1 , weights=Phi_1)
    
    beta_hat_1=model_1$coefficients
    Beta_Hat_unif_500[i]=beta_hat_1
  }
  return(Beta_Hat_unif_500)
}

###beta_unif with unif-algo

weighted_algo <- function(r,Pi) {
  Beta_Hat_weighted_500 = c()
  for (i in 1:draws){
    index_2=sample(c(1:n), size = r, replace = T,prob=Pi_weighted) 
    #Phi_2=diag(Pi_weighted[index_2])
    X_star_2=X[index_2]
    Y_star_2=Y[index_2]
    Phi_2=1/ (Pi_weighted[index_2]^{1/2})
    model_2 <- lm(Y_star_2 ~ 0 + X_star_2 , weights=Phi_2)
    
    beta_hat_2=model_2$coefficients
    Beta_Hat_weighted_500[i]=beta_hat_2
  }
  #mean(MSE_2_200)
  #sd(MSE_2_200)
  #rowMeans(Beta_Hat_2_200)
  #rowSds(Beta_Hat_2_200)
  #norm(beta-rowMeans(Beta_Hat_2_200), type = "2")
  
  return(Beta_Hat_weighted_500)
}
