library(dplyr)
library(LaplacesDemon)

# Sampling Model:		y_{it} ~ N(X_{it} beta_i, sig.sq)
#	Priors:			
#			beta_i ~ N(mu, Sigma)
#			(mu,Sigma) ~ Normal-Inv-Wishart(m,v,C,d)
#			sig.sq ~ 1/sig.sq (Jeffrey's Prior)

gibbsCheese <- function(y, X, allStores, m, v, C, d,iter= 1000, burn=500, thin=2) {
  
  n <- length(y); p <- ncol(X); nstores <- nlevels(allStores)
  beta <- array(0,dim=c(nstores,p,iter))
  mu <- matrix(0,iter,p); Sigma <- array(0, dim=c(p, p, iter))
  sig.sq <- rep(0,iter)
  
  beta[,,1] <- rep(0, nstores*p);  mu[1,] <- rep(mean(y),p); Sigma[,,1] <- diag(p); sig.sq[1] <- 1
  
  for (iter in 2:iter) {
    
    #Update beta
    Sig.inv <- solve(Sigma[,,iter-1])
    for (ns in 1:nstores){	
      
      bStore <- which(ns==as.integer(allStores))
      Xi <-  X[bStore,]; yi <- y[bStore]; nt <- length(yi)
      beta.post.var <- solve(Sigma[,,iter-1] + (nt/sig.sq[iter-1]) * crossprod(Xi))
      beta.post.mean <- beta.post.var %*% (Sig.inv %*% mu[iter-1,] + (nt/sig.sq[iter-1]) * t(Xi) %*% yi )
      beta[ns,,iter] <- rmvnorm(1, beta.post.mean, beta.post.var)
    }
    
    #Update sig.sq
    beta.lookup <- beta[idx$store,,iter]; SS <- sum((y - X %*% t(beta.lookup))^2)
    sig.sq[iter] <- 1/rgamma(1, n/2, SS/2)
    
    #Update mu and Sigma values for NIW
    beta.bar <- colMeans(beta[,,iter])
    S <- t(beta[,,iter] - beta.bar)	%*% (beta[,,iter] - beta.bar)
    
    mn = (v*m + nstores*beta.bar)/(v + nstores)
    vn = v + nstores
    Cn = C + S + (v*nstores / (v+nstores)) * (beta.bar - m) %*% t(beta.bar - m)
    dn = d + nstores
    
    niv.sample <- rnorminvwishart(1, mn, vn, Cn, dn)
    mu[iter,] <- niv.sample$mu; Sigma[,,iter] <- niv.sample$Sigma
  }
  
  thinseq <- seq(1,iter - 1, by=thin)
  beta = beta[,,-burn]; Sigma = Sigma[,,-burn]; mu = mu[-burn,]; sig.sq = sig.sq[-burn]
  beta = beta[,,thinseq]; Sigma = Sigma[,,thinseq]; mu = mu[thinseq,]; sig.sq = sig.sq[thinseq]
  
  #posterior medians
  beta.post.median    <- apply(beta, 2, median)
  sig.sq.post.median  <- median(sig.sq)
  mu.post.median      <- apply(mu, 2, median)
  Sigma.post.median   <- apply(Sigma,2, median)
  
  list(beta = beta,Sigma = Sigma, mu = mu, sig.sq = sig.sq,
              beta.post.median = beta.post.median, 
              sig.sq.post.median = sig.sq.post.median,
              mu.post.median = mu.post.median,
              Sigma.post.median = Sigma.post.median)
}

data <- read.csv('cheese.csv')
data <- data %>% 
  mutate(logP = log(price), logQ = log(vol), disp = disp) %>%
  as.data.frame() %>%
  dplyr::select(-one_of("vol","price"))

y  <- data$logQ; X <- model.matrix(logQ ~ 1 + logP + disp, data=data)
p <- 3; m <- rep(0, 1); v = 1; C = diag(p); d = p + 1

mcoutput <- gibbsCheese(y, X, data$store, m, v, C, d)


