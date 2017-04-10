rm(list = ls())

library(dplyr)
library(mvtnorm)
library(MCMCpack)
library(robustbase)

# Sampling Model:		y_{it} ~ N(X_{it} beta + Z_{it} gamma_i, sig.sq)
#	Priors:			
#			gamma_i ~ N(0, Sigma)
#			Sigma   ~ Inv-Wishart(d,C)
#			sig.sq  ~ 1/sig.sq (Jeffrey's Prior)
#     beta ~ 1 (flat prior)

gibbsCheese <- function(y, X, Z, allStores, d, C, iter= 1000, burn=2, thin=2) {
  
  n <- length(y); p <- ncol(Z); nstores <- nlevels(allStores); 
  
  #setup the structure for the chain
  beta <- array(0, dim=c(ncol(X),iter)); gamma <- array(0,dim=c(nstores, ncol(Z), iter)); 
  Sigma <- array(0, dim=c(p, p, iter)); sig.sq <- rep(0,iter)
  
  #initialize the chain
  beta[,1] <- rep(0, p) ;gamma[,,1] <- rep(0, nstores*p);  Sigma[,,1] <- diag(p)
  sig.sq[1] <- 1
  
  pb <- txtProgressBar(min = 2, max = iter, style = 3)
  
  for (iter in 2:iter) {
    
    # Update fixed effect for beta
    beta.post.var <- solve(
      sig.sq[iter-1]^-1 * Reduce('+',
        lapply(1:nstores, function(ns) {
          bStore <- which(ns==as.integer(allStores)); 
          Xi <- X[bStore,,drop=FALSE]; crossprod(Xi)
        })
      ))
    
    beta.post.mean <- (beta.post.var/sig.sq[iter-1]) %*% t(Reduce('+',
      lapply(1:nstores, function(ns) {
        bStore <- which(ns==as.integer(allStores)); 
        Zi <-  Z[bStore,]; yi <- y[bStore]; Xi <- X[bStore,,drop=FALSE]
        t(yi - Zi %*% gamma[ns,,iter-1]) %*% Xi
      })
    ))
    
    beta[,iter] <- rmvnorm(1, beta.post.mean, beta.post.var)
    
    # Update random effect gamma for each store
    Sig.inv <- solve(Sigma[,,iter-1])
    
    gamma[,,iter] <- 
      t(sapply(1:nstores, function(ns) {	
        
        bStore <- which(ns==as.integer(allStores))
        Xi <- X[bStore,,drop=FALSE]; Zi <-  Z[bStore,]; yi <- y[bStore]
        
        gamma.post.var <- solve((sig.sq[iter-1])^-1 * crossprod(Zi) + Sig.inv)
        gamma.post.mean <- 
          gamma.post.var %*% 
          (sig.sq[iter-1]^-1 * crossprod(Zi, yi - Xi %*% beta[,iter,drop=FALSE]))
        
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
    S <- Reduce('+',lapply(1:nstores, function(ns) {
      tcrossprod(gamma[ns,,iter])
    }))
    
    Cn = C + S; dn = d + nstores
    Sigma[,,iter] <- riwish(dn, Cn)
    
    setTxtProgressBar(pb, iter)
  }
  
  thinseq <- seq(1,iter - 1, by=thin)
  
  beta <- beta[,-burn]; gamma <- gamma[,,-burn];   Sigma <- Sigma[,,-burn];   sig.sq <- sig.sq[-burn]
  beta <- beta[,thinseq]; gamma <- gamma[,,thinseq]; Sigma <- Sigma[,,thinseq]; sig.sq <- sig.sq[thinseq]
  
  #posterior medians
  beta.post.median     <- apply(beta, 1, median)
  gamma.post.median    <- apply(gamma, c(1,2), median)
  sig.sq.post.median   <- median(sig.sq)
  Sigma.post.median    <- apply(Sigma, c(1,2), mean)
  
  list(beta.post = beta.post.median, 
              gamma.post = gamma.post.median, 
              sig.sq.post = sig.sq.post.median,
              Sigma.post = Sigma.post.median)
}

data <- read.csv('cheese.csv')
data <- data %>% 
  mutate(logP = log(price), logQ = log(vol), disp = disp) %>%
  as.data.frame() %>%
  dplyr::select(-one_of("vol","price"))

y  <- data$logQ; X <- Z <- model.matrix(logQ ~ 1 + logP + disp + logP:disp, data=data)
p <- d <- ncol(Z); C = diag(p)
mcoutput <- gibbsCheese(y, X, Z, data$store, d, C)

cols <- c('blue','red')
xgrid <- seq(min(X[,2]), max(X[,2]), length.out = 50)
stores <- as.integer(data$store)

par(mfrow = c(2,2))

for (i in c(3, 7, 18, 25)){
  
  plot(X[stores==i,2], y[stores == i], main = as.character(stores[i]), 
       col = cols[X[stores==i, "disp"]+1], xlab = "logPrice", ylab="logQuantity")
  lines(xgrid, (mcoutput$beta.post[1] + mcoutput$gamma.post[i,1]) + 
          xgrid * (mcoutput$beta.post[2] + mcoutput$gamma.post[i,2]), col = cols[1], lwd = 2)
  lines(xgrid, (mcoutput$beta.post[1] + mcoutput$beta.post[3] + 
                  mcoutput$gamma.post[i,1]+ mcoutput$gamma.post[i,3]) + 
          xgrid * (mcoutput$beta.post[2] + mcoutput$beta.post[4] + 
                     mcoutput$gamma.post[i,2] + mcoutput$gamma.post[i,4]), 
        col = cols[2], lwd = 2)              
  
  fit = lmrob(y[stores == i] ~ -1 + X[stores==i,])
  betas.lm <- fit$coefficients
  lines(xgrid, (betas.lm[1]) + xgrid * (betas.lm[2]), col = cols[1], lwd = 1, lty = 2)
  lines(xgrid, (betas.lm[1]+ betas.lm[3]) + xgrid * (betas.lm[2] + betas.lm[4]), 
        col = cols[2], lwd = 1, lty = 2)              
}

par(mfrow = c(1,1))
