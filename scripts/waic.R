#calculating WAIC

# compute wAIC for model with covariate
library(R2jags)
samples.m1 <- jags.samples(cjs.rain.site$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m1$p_waic <- samples.m1$WAIC
samples.m1$waic <- samples.m1$deviance + samples.m1$p_waic
tmp <- sapply(samples.m1, sum)
waic.m1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


################
samples.m2 <- jags.samples(cjs.springrain.site$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m2$p_waic <- samples.m2$WAIC
samples.m2$waic <- samples.m2$deviance + samples.m2$p_waic
tmp <- sapply(samples.m2, sum)
waic.m2 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
##############3333

samples.m3 <- jags.samples(cjs.summer.rain.site$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m3$p_waic <- samples.m3$WAIC
samples.m3$waic <- samples.m3$deviance + samples.m3$p_waic
tmp <- sapply(samples.m3, sum)
waic.m3 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
####################

# wAIC of m1 the model with covariate << wAIC of m0 intercept only
d<-data.frame(totalrain = waic.m1, springrain = waic.m2, summerrain = waic.m3)
write.csv(d,'./output/waic.rain.csv', row.names = F)
