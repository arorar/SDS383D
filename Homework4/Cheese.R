rm(list=ls())

library(lme4)
library(dplyr)
library(mosaic)

data <- read.csv('cheese.csv')

modData <- data %>% 
  group_by(store) %>%
  mutate(logP = log(price), logQ = log(vol)) %>%
  summarise(ni = n(), meanp=mean(logP), meanq = mean(logQ)) %>%
  as.data.frame()

par(mfrow = c(2,2))
plot(modData$meanp, pch=19,main='Avg Price by Store Over All Weeks', 
     xlab='Store',ylab='Price')
plot(modData$meanq, pch=19,main='Avg Volume by Store Over All Weeks', 
     xlab='Store',ylab='Vol')
plot(modData$meanq, modData$meanp, pch=19, 
     main='Avg Volume vs Price by Store Over All Weeks', 
     xlab='Vol',ylab='Price')
par(mfrow = c(1,1))

modData <- data %>% 
  mutate(logP = log(price), logQ = log(vol)) %>%
  as.data.frame()

hlm <- lmer(logQ ~  (-1 + ( -1 + logP + disp + disp:logP | store)), data = modData)

resid <- ranef(hlm, condVar = T)
dotplot(resid, scales=list(cex=c(.5,.5)),layout=c(2,1), main=T, 
        main.title='Random Effects by by Store',)


