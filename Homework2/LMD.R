rm(list=ls())

gdp <- read.csv(file = 'gdpgrowth.csv')

X <- cbind(1, gdp$DEF60); y <- gdp$GR6096
n <- nrow(gdp); p <- ncol(X) - 1

Lambda <- diag(rep(1, n))
m <- rep(0, p + 1); K <- diag(c(0.005, 0.005))
d <- 0.01; eta <- 0.01

XtLambdaX <- t(X) %*% Lambda %*% X

# Update the hyperparameters for the posterior
temp <- (K %*% m + t(X) %*% Lambda %*% y)
K.star <- K + XtLambdaX; temp.K.star <- solve(K.star)
mu.star <- temp.K.star %*% temp
d.star <- d + n + p
eta.star <- t(y) %*% Lambda %*% y + t(m) %*% K %*% m + eta - 
  t(temp) %*% temp.K.star %*% temp

linMod = lm(y ~ X-1)
plot(X[,2], y, pch = 1, cex = 0.8, asp = 1, xlab = 'Defense Spending', 
     ylab = 'GDP Growth',
     main = "Comparison of Bayesian vs Frequentist LM")
abline(coef = c(linMod$coefficients[1], linMod$coefficients[2]), col = "red")
abline(coef = c(mu.star[1], b = mu.star[2]), col = 'blue')
legend('topright', bty = "n", legend = c('LM','BayesLM'), lty = 1, 
       col = c('red', 'blue'))



