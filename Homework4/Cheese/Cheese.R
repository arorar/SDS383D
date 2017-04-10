rm(list=ls())

library(lme4)
library(dplyr)
library(mosaic)

data <- read.csv('cheese.csv')

modData <- data %>% 
  mutate(logP = log(price), logQ = log(vol)) %>%
  dplyr::select(-one_of(c("price","vol"))) %>%
  as.data.frame()

par(mfrow = c(1,2))
boxplot(logP ~ disp, data = modData, col = 'gray', pch = 16, main="log price vs display")
boxplot(logQ ~ disp, data = modData, col = 'gray', pch = 16, main="log quant vs display")
par(mfrow = c(1,1))

xyplot(logQ ~ logP | store, data = modData, group = disp, col = c('blue','red'), 
       col.line = c('blue','red'), cex = 0.5, type = c("r", "p"), 
       par.strip.text=list(cex=.5))

hlm <- lmer(logQ ~  (1 + logP + disp + disp:logP | store), data = modData)
summary(hlm)

resid <- ranef(hlm, condVar = T)
dotplot(resid, scales=list(cex=c(.5,.5)),layout=c(4,1), main=T, 
        main.title='Random Effects by by Store')


