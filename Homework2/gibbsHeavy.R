rm(list = ls())

library(mvtnorm)

gibbsSampler.heavytails <- function(X, y, m, K, d, eta, h, niter = 10000,
                             burnin = 1000, thin = 2){

    n <- nrow(X); p <- ncol(X) - 1

    # Setup
    betas   <- matrix(NA, nrow = niter, ncol = p+1)
    lambdas <- matrix(NA, nrow = niter, ncol = n)
    omegas  <- rep(NA, niter)

    # Set initial values
    betas[1,] <- rep(0, p+1); lambdas[1,] <- rep(2, n); omegas[1] <- 2

    for (iter in 2:niter){
      
      Lambda <- diag(lambdas[iter-1,]); 
      XtLambdaX <- t(X) %*% Lambda %*% X
      
      # Update the hyperparameters for the posterior
      temp <- (K %*% m + t(X) %*% Lambda %*% y)
      K.star <- K + XtLambdaX; K.star.inv <- solve(K.star)
      mu.star <- K.star.inv %*% temp
      d.star <- d + n + p
      eta.star <- t(y) %*% Lambda %*% y + t(m) %*% K %*% m + eta -
        t(temp) %*% K.star.inv %*% temp
      
      # Update beta step
      betas[iter,] <- rmvnorm(1, mean = mu.star, sigma = K.star.inv/omegas[iter-1])

      # Update omega step
      omegas[iter] <- rgamma(1, d.star/2, eta.star/2)

      # Update lambda step
      lambda.start.rate <- (1/2) * (omegas[iter] * (y - X %*% betas[iter,])^2 + h)
      lambdas[iter,] <- rgamma(n, rep((h+1)/2, n), as.numeric(lambda.start.rate))
    }

    # Thin the chains to avoid autocorrelation
    vals <- seq(burnin + 1, niter, by = thin)
    betas   <- betas[vals,]; lambdas <- lambdas[vals,]; omegas  <- omegas[vals]

    return(list(betas = betas, lambdas = lambdas, omegas = omegas))
}

# Initialize
gdp <- read.csv(file = 'gdpgrowth.csv')
X <- cbind(1, gdp$DEF60); y <- gdp$GR6096
n <- nrow(gdp); p <- ncol(X) - 1

# Set hyperparameters
m <- rep(0, p + 1); K <- diag(c(0.0005, 0.0005)); d <- 0.01; eta <- 0.01; h <- 0.1

# Gibbs. Will get 5000 samples
niter <- 12000; burnin <- 2000; thin <- 2
posterior <- gibbsSampler.heavytails(X, y, m, K, d, eta, h)
betas <- posterior$betas; lambdas <- posterior$lambdas; omegas <- posterior$omegas

# Select the first 10 low precision/(high variance) observations 
cml <- colMeans(lambdas); scml <- sort(cml)[1:10]; idx <- which(cml %in% scml)

# Plots
linMod <- lm(y ~ X-1)
plot(X[,2], y, xlab = 'Defense Spending', ylab = 'GDP Growth',
     main = "Comparison of Heavy Tailed Bayesian vs Frequentist LM")
abline(coef = linMod$coefficients[1:2], col = "red")
abline(coef = colMeans(betas)[1:2], col = "blue")

# Display those countries
points(X[idx,2], y[idx], pch = 16, cex = 0.9, col = 'darkgreen')
text(X[idx,2], y[idx], as.character(gdp$CODE[idx]), cex=.75, adj=1.25)

legend('bottomright', bty = "n", legend = c('LM','HeavyTailBayesLM'), lty = 1,
       col = c('red', 'blue'))





