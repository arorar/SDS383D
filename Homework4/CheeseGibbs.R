library(dplyr)
library(LaplacesDemon)

# Sampling Model:		y_{it} ~ N(X_{it} beta_i, sig.sq)
#	Priors:			
#			beta_i ~ N(mu, Sigma)
#			(mu,Sigma) ~ Normal-Inv-Wishart(m,v,C,d)
#			sig.sq ~ 1/sig.sq (Jeffrey's Prior)

gibbsCheese <- function(y, X, idx, m, v, C, d,iter= 10000, burn=1000, thin=2) {
  
  n <- length(y); p <- ncol(X); stores <- length(unique(idx$store))
  beta <- array(0,dim=c(stores,p,iter))
  mu <- matrix(0,iter,p); Sigma <- array(0, dim=c(p, p, iter))
  sig.sq <- rep(0,iter)
  
  beta[,,1] <- rep(0, stores*p);  mu[1,] <- rep(mean(y),p); Sigma[,,1] <- diag(p); sig.sq[1] <- 1
  
  for (k in 2:iter) {
    
    Sig.inv <- solve(Sigma[,,k-1])
    for (j in 1:stores){	
      
      bStore <- which(j==idx$storenum)
      Xi <-  X[bStore,]; yi <- y[bStore]; nt <- length(yi)
      beta.var <- solve( Sigma[,,k-1] + (nt/sig.sq[k-1]) * crossprod(Xi)	)
      beta.mean <- beta.var %*% (Sig.inv %*% mu[k-1,] + (nt/sig.sq[k-1]) * t(Xi) %*% yi )
      beta[j,,k] <- rmvnorm(1, beta.mean, beta.var)
    }
    
    #Update sig.sq
    beta.lookup <- beta[idx$store,,k]; SS <- sum((y - X %*% t(beta.lookup))^2)
    sig.sq[k] <- 1/rgamma(1, n/2, SS/2)
    
    #Update mu and Sigma values for NIW
    beta.bar <- colMeans(beta[,,k])
    S <- t(beta[,,k] - beta.bar)	%*% (beta[,,k] - beta.bar)
    
    mn = (v*m + stores*beta.bar)/(v + stores)
    vn = v + stores
    Cn = C + S + (v*stores / (v+stores)) * (beta.bar - m) %*% t(beta.bar - m)
    dn = d + stores
    
    niv.sample <- rnorminvwishart(1, mn, vn, Cn, dn)
    mu[k,] <- niv.sample$mu; Sigma[,,k] <- niv.sample$Sigma
  }
  
  thinseq <- seq(1,iter - 1, by=thin)
  beta = beta[,,-burn]; Sigma = Sigma[,,-burn]; mu = mu[-burn,]; sig.sq = sig.sq[-burn]
  beta = beta[,,thinseq]; Sigma = Sigma[,,thinseq]; mu = mu[thinseq,]; sig.sq = sig.sq[thinseq]
  
  #posterior medians
  beta.post.mean = apply(beta,c(1,2),median)
  Sigma.post.mean = apply(Sigma,c(1,2),median)
  mu.post.mean = apply(mu,2,median)
  sig.sq.post.mean = median(sig.sq)
  
  #Return results.
  list(beta=beta,Sigma=Sigma,mu=mu,sig.sq=sig.sq,
              beta.post.mean=beta.post.mean, 
              mu.post.mean=mu.post.mean,
              Sigma.post.mean=Sigma.post.mean,
              sig.sq.post.mean=sig.sq.post.mean)
}

data <- read.csv('cheese.csv')
data <- data %>% 
  mutate(logP = log(price), logQ = log(vol), disp = disp, 
         week = ave(1:nrow(data), store, FUN=seq),
         storenum = ave(1:nrow(data), week, FUN=seq)) %>%
  as.data.frame() %>%
  dplyr::select(-one_of("vol","price"))

y  <- data$logQ; X <- model.matrix(logQ ~ 1 + logP + disp, data=data)
idx <- cbind.data.frame(store=data$store, storenum=data$storenum, week=data$week)	
p <- 3; m <- rep(0, 0); v = 1; C = diag(p); d = p + 1

output <- gibbsCheese(y, X, idx, m, v, C, d)
