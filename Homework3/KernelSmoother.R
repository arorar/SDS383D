rm(list = ls())

set.seed(100)

kernelFunc <- function(distMethod, bandwidth, x, x0) {
   
  d <- as.numeric(dist(c(x,x0), method = distMethod))
  dnorm(d, sd = bandwidth)
  
}

weightFunc <- function(x0, x, kfunc, bandwidth=1) {
  
  k <- 
    sapply(x, function(xi) {
      kfunc("euclidean", bandwidth, xi, x0)
    })
  
  w <- k/bandwidth
  w/sum(w)
}

x <- 1:100; y <- 2*sin(0.25*x) + sin(rnorm(length(x)))
x <- scale(x, scale = FALSE); y <- scale(y, scale = FALSE)



for(i in -20:20) {

  x0 <- i
  w <- weightFunc(x0, x, kernelFunc, 1.5)
  yhat <- crossprod(w, y)
  
  plot(x, y, type = "l")
  points(i, yhat, pch = 16, col = "red")
  Sys.sleep(1)
}
