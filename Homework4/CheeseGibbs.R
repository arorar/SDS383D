rm(list = ls())

library(dplyr)
library(mvtnorm)
library(MCMCpack)

# No fixed effects. Missing beta
# Sampling Model:		y_{it} ~ N(Z_{it} gamma_i, sig.sq)
#	Priors:			
#			gamma_i ~ N(0, Sigma)
#			Sigma   ~ Inv-Wishart(d,C)
#			sig.sq  ~ 1/sig.sq (Jeffrey's Prior)

gibbsCheese <- function(y, Z, allStores, d, C, iter= 100, burn=20, thin=2) {
  
  n <- length(y); p <- ncol(Z); nstores <- nlevels(allStores); sig.sq <- rep(0,iter)
  gamma <- array(0,dim=c(nstores,p,iter)); Sigma <- array(0, dim=c(p, p, iter))
  
  gamma[,,1] <- rep(0, nstores*p);  Sigma[,,1] <- diag(p); sig.sq[1] <- 1
  
  for (iter in 2:iter) {
    
    #Update gamma
    Sig.inv <- solve(Sigma[,,iter-1])
    
    gamma[,,iter] <- 
      t(sapply(1:nstores, function(ns) {	
        
        bStore <- which(ns==as.integer(allStores))
        Zi <-  Z[bStore,]; yi <- y[bStore]; nt <- length(yi)
        gamma.post.var <- solve(Sigma[,,iter-1] + (sig.sq[iter-1])^-1 * crossprod(Zi))
        gamma.post.mean <- gamma.post.var %*% (sig.sq[iter-1]^-1 * t(Zi) %*% yi)
        rmvnorm(1, gamma.post.mean, gamma.post.var)
      }))
    
    #Update sig.sq
    gamma.store <- gamma[as.integer(allStores),,iter]; SS <- sum((y - Z %*% t(gamma.store))^2)
    sig.sq[iter] <- 1/rgamma(1, n/2, SS/2)
    
    #Update Sigma values for IW
    S <- 0; for (s in 1:nstores) S <- S + gamma[s,,iter] %*% t(gamma[s,,iter])
    Cn = C + S
    dn = d + nstores
    
    Sigma[,,iter] <- riwish(dn, Cn)
  }
  
  thinseq <- seq(1,iter - 1, by=thin)
  gamma = gamma[,,-burn]; Sigma = Sigma[,,-burn];sig.sq = sig.sq[-burn]
  gamma = gamma[,,thinseq]; Sigma = Sigma[,,thinseq]; sig.sq = sig.sq[thinseq]
  
  #posterior medians
  gamma.post.median    <- apply(gamma, 2, median)
  sig.sq.post.median  <- median(sig.sq)
  Sigma.post.median   <- apply(Sigma, c(1,2), mean)
  
  list(gamma.post.median = gamma.post.median, 
              sig.sq.post.median = sig.sq.post.median,
              Sigma.post.mean = Sigma.post.median)
}

data <- read.csv('cheese.csv')
data <- data %>% 
  mutate(logP = log(price), logQ = log(vol), disp = disp) %>%
  as.data.frame() %>%
  dplyr::select(-one_of("vol","price"))

y  <- data$logQ; Z <- model.matrix(logQ ~ logP + disp, data=data)
p <- d <- 3; C = diag(p)
mcoutput <- gibbsCheese(y, Z, data$store, d, C)


