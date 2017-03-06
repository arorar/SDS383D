rm(list=ls())

localLinearSmoother <- function(x, y, x_star, h=1, kernel = dnorm){
  
  knl <- match.fun(kernel); n_star <- length(x_star)
  weights <- matrix(NA , nrow = length(x), ncol = length(x_star))	
  
  yhat <- 
    sapply(1:n_star, function(b) {
    
    s <- sapply(1:2, function(j) {
      sum( knl((x_star[b]-x)/h) * (x-x_star[b])^j )
    }) 
    
    w = knl((x_star[b]-x)/h) * (s[2] - (x-x_star[b])*s[1])
    w = w / sum(w)
    crossprod(w,y)
  })
  
  yhat
}

utilities <- read.csv('utilities.csv', header=TRUE)
x <- utilities$temp; y <- utilities$gasbill/utilities$billingdays
x <- scale(x, scale = FALSE); y <- scale(y, scale = FALSE)

bGrid <- seq(0.01, 5, by = 0.05)

pred_err <- sapply(bGrid, function(h) {
  yhat <- localLinearSmoother(x = x, y = y, x_star = x, h=h)
  sum((y - yhat)^2)/length(x)
})

plot(bGrid, pred_err, type = "l", main = "Bandwidth vs prediction Error")
opt_bwidth <- bGrid[which.min(pred_err)]
yhat <- localLinearSmoother(x, y, x, h=opt_bwidth)

plot(x, y, 
     main=paste('Local Linear Smoothing with Gaussian Kernel, Bandwidth = ', 
                sep="", opt_bwidth))	
lines(sort(x), yhat[order(x)],col='blue',lwd=2, lty = 1)

yhat <- localLinearSmoother(x, y, x, h = opt_bwidth); resids <- (y - yhat)
plot(x, resids, main = 
       paste('Local Linear Smoothing Residuals. Bandwidth = ',sep="", opt_bwidth))	
lfit <- loess(resids~x, col='blue',lwd=2, lty = 1)
lines(sort(x), lfit$fitted[order(x)], col='blue',lwd=2, lty = 1)