rm(list=ls())

library(dplyr)

mathTestGibbsSampler <- function(data, Niter = 22000, burnin = 2000, thin = 2) {
  
  querry <- data %>% 
    group_by(school) %>% 
    summarise(ni = n(), mean=mean(mathscore)) 
  
  ybar <- querry$mean; ni <- querry$ni; n <- sum(ni); p <- length(ni)
  
  # Initialize the chain
  thetas.chain <- array(NA, dim = c(Niter, p))
  mu.chain <- sigma2.chain <- tau2.chain <- array(NA, dim = Niter)
  mu.chain[1] <- 0; sigma2.chain[1] <- 1; tau2.chain[1] <- 1
  
  for (i in 2:Niter){
    
    # Posterior thetas
    var.post <- tau2.chain[i-1] * sigma2.chain[i-1] / (ni * tau2.chain[i-1] + 1)
    mean.post <- (mu.chain[i-1] + tau2.chain[i-1] * ni * ybar) / (ni * tau2.chain[i-1] + 1)
    thetas.chain[i,] <- rnorm(p, mean.post, sqrt(var.post))
    
    # Posterior mu
    theta.bar <- mean(thetas.chain[i,])
    mu.chain[i] <- rnorm(1, theta.bar, sqrt(sigma2.chain[i-1] * tau2.chain[i-1] / p))
    
    # Posterior sigma2
    S.theta <- sum((thetas.chain[i,] - mu.chain[i])^2)
    S.mathscore <- sum((data$mathscore - rep(thetas.chain[i,], times = ni))^2)
    rate.new <- (1/2) * (S.mathscore + S.theta / tau2.chain[i-1])
    sigma2.chain[i] <- 1/rgamma(1, (n + p)/2, rate.new)
    
    # Posterior tau2
    rate.new <- S.theta / (2 * sigma2.chain[i])
    tau2.chain[i] <- 1/rgamma(1, p/2 - 1, rate.new)
  }
  
  # Thin the chains
  acutalSeq <- seq(burnin + 1, Niter, by = thin)
  thetas.chain <- thetas.chain[acutalSeq,]; mu.chain   <- mu.chain[acutalSeq]
  sigma2.chain <- sigma2.chain[acutalSeq];  tau2.chain <- tau2.chain[acutalSeq]
  
  list(theta.hat  = colMeans(thetas.chain), mu.hat = mean(sigma2.chain),
       sigma2.hat = mean(sigma2.chain), tau2.hat = mean(tau2.chain), 
       ybar = ybar, ssize = ni)
}

math <- read.csv(file = 'mathtest.csv')
mathTestGibbs <- mathTestGibbsSampler(data = math)
kappa <- (mathTestGibbs$ybar - mathTestGibbs$theta.hat) / mathTestGibbs$ybar
scatter.smooth(mathTestGibbs$ssize, abs(kappa), xlab="Sample size", 
     ylab="Shrinkage towards Mean", 
     main = "Shrinkage of Group Mean Towards Global Mean",
     ylim = c(-0.02,0.15), 
     col="red", degree = 2)
abline( h=0, lty  =2)
