xiaohan_knn <- function(input, point, k){
  m = dim(input)[1]
  n = dim(input)[2]
  distmat  = input - matrix(rep(c(point, 0), n), m, n)
  distance = colSums(distmat[1:(m-1),]^2)
  tosort = rbind(distance, distmat[m,])
  sorted = tosort[, order(tosort[1,])]
  return(round(sum(sorted[2,1:k])/k))
}

# ##Testing part
# library(ggplot2)
# samples = t(mvrnorm(2000, c(0,0), diag(rep(1,2))))
# samplest = t(mvrnorm(2000, c(2,2), diag(rep(1,2))))
# samples = cbind(samples, samplest)
# samplewtag = rbind(samples, c(rep(1,2000), rep(0,2000)))
# 
# 
# testing = xiaohan_knn(samplewtag, c(1,1), 5)
# plotdata = cbind(samplewtag, c(1,1,testing+3))
# 
# 
# plotdata = as.data.frame(t(plotdata))
# plotdata[,3] = as.factor(plotdata[,3])
# colnames(plotdata) <- c("xaxis", "yaxis", "catagories")
# tp <- ggplot(data = plotdata) +
#   geom_point(aes(x = xaxis, y = yaxis, color = catagories), size = .7)
# tp
# 
