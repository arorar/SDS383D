rm(list = ls())

library(dplyr)
library(mvtnorm)
library(MCMCpack)

# Sampling Model:		y_{it} ~ N(X_{it} beta + Z_{it} gamma_i, sig.sq)
#	Priors:			
#			gamma_i ~ N(0, Sigma)
#			Sigma   ~ Inv-Wishart(d,C)
#			sig.sq  ~ 1/sig.sq (Jeffrey's Prior)
#     beta ~ 1 (flat prior)

gibbsCheese <- function(y, X, Z, allStores, d, C, iter= 10, burn=2, thin=2) {
  
  n <- length(y); p <- ncol(Z); nstores <- nlevels(allStores); 
  
  #setup the structure for the chain
  beta <- array(0, dim=c(ncol(X),iter)); gamma <- array(0,dim=c(nstores, ncol(Z), iter)); 
  Sigma <- array(0, dim=c(p, p, iter)); sig.sq <- rep(0,iter)
  
  #initialize the chain
  gamma[,,1] <- rep(0, nstores*p);  Sigma[,,1] <- diag(p)
  beta[,1] <- sig.sq[1] <- 1
  
  for (iter in 2:iter) {
    
    # Update fixed effect for beta
    beta.post.var <- solve(
      sig.sq[iter-1]^-1 * sum(
        sapply(1:nstores, function(ns) {
          bStore <- which(ns==as.integer(allStores)); 
          Xi <- X[bStore,,drop=FALSE]; crossprod(Xi)
        })
      ))
    
    beta.post.mean <- beta.post.var/sig.sq[iter-1] %*% sum(
      sapply(1:nstores, function(ns) {
        bStore <- which(ns==as.integer(allStores)); 
        Zi <-  Z[bStore,]; yi <- y[bStore]; Xi <- X[bStore,,drop=FALSE]
        t(yi - Zi %*% gamma[ns,,iter-1]) %*% Xi
      })
    )
    
    beta[,iter] <- rmvnorm(1, beta.post.mean, beta.post.var)
    
    # Update random effect gamma for each store
    Sig.inv <- solve(Sigma[,,iter-1])
    
    gamma[,,iter] <- 
      t(sapply(1:nstores, function(ns) {	
        
        bStore <- which(ns==as.integer(allStores))
        Xi <- X[bStore,,drop=FALSE]; Zi <-  Z[bStore,]; yi <- y[bStore]
        
        gamma.post.var <- solve(Sigma[,,iter-1] + (sig.sq[iter-1])^-1 * crossprod(Zi))
        gamma.post.mean <- 
          gamma.post.var %*% 
          (sig.sq[iter-1]^-1 * t(Zi) %*% (yi - Xi %*% beta[,iter,drop=FALSE]))
        
        rmvnorm(1, gamma.post.mean, gamma.post.var)
      }))
    
    
    # Update sig.sq
    
    SS <- sum(
      sapply(1:nstores, function(ns) {
        bStore <- which(ns==as.integer(allStores)); 
        Zi <-  Z[bStore,]; yi <- y[bStore]; Xi <- X[bStore,,drop=FALSE]
        crossprod(yi - Zi %*% gamma[ns,,iter] - Xi %*% beta[,iter,drop=FALSE])
      })
    )
    
    sig.sq[iter] <- 1/rgamma(1, n/2, SS/2)
    
    #Update Sigma values for IW
    S <- 0; for (s in 1:nstores) S <- S + gamma[s,,iter] %*% t(gamma[s,,iter])
    Cn = C + S; dn = d + nstores
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
mcoutput <- gibbsCheese(y, Z[,1, drop = FALSE], Z, data$store, d, C)


