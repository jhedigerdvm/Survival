#calculating WAIC

# compute wAIC for model with covariate
library(R2jags)

################survival as a function of birth year rain, birth site, age, and the interaction bt bs and birthyear rain
samples.m1 <- jags.samples(cjs.sitexrain$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m1$p_waic <- samples.m1$WAIC
samples.m1$waic <- samples.m1$deviance + samples.m1$p_waic
tmp <- sapply(samples.m1, sum)
waic.m1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


#######survival as a function of birth year rain, birth site, age, and NO interaction bt bs and birthyear rain
samples.m2 <- jags.samples(cjs.site.rain$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m2$p_waic <- samples.m2$WAIC
samples.m2$waic <- samples.m2$deviance + samples.m2$p_waic
tmp <- sapply(samples.m2, sum)
waic.m2 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

#survival as a function of birth year rain, birth site, age, and the interaction bt bs and birthyear rain and site and age
samples.m3 <- jags.samples(cjs.site.rain.age$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m3$p_waic <- samples.m3$WAIC
samples.m3$waic <- samples.m3$deviance + samples.m3$p_waic
tmp <- sapply(samples.m3, sum)
waic.m3 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


#survival as a function of birth year rain, birth site, age, and the interaction of bs and age
samples.m4 <- jags.samples(cjs.sitexage$model,
                           c("WAIC","deviance"),
                           type = "mean",
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m4$p_waic <- samples.m4$WAIC
samples.m4$waic <- samples.m4$deviance + samples.m4$p_waic
tmp <- sapply(samples.m4, sum)
waic.m4 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
# 
# # ############## summer rain model
# samples.m5 <- jags.samples(cjs.summer.rain.site$model,
#                            c("WAIC","deviance"),
#                            type = "mean",
#                            n.iter = 5000,
#                            n.burnin = 1000,
#                            n.thin = 1)
# 
# samples.m5$p_waic <- samples.m5$WAIC
# samples.m5$waic <- samples.m5$deviance + samples.m5$p_waic
# tmp <- sapply(samples.m5, sum)
# waic.m5 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
# 
# # ############## summer rain model
# samples.m6 <- jags.samples(cjs.rain.site.ran.age$model,
#                            c("WAIC","deviance"),
#                            type = "mean",
#                            n.iter = 5000,
#                            n.burnin = 1000,
#                            n.thin = 1)
# 
# samples.m6$p_waic <- samples.m6$WAIC
# samples.m6$waic <- samples.m6$deviance + samples.m6$p_waic
# tmp <- sapply(samples.m6, sum)
# waic.m6 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


# wAIC of m1 the model with covariate << wAIC of m0 intercept only
d<-data.frame(sitexrain=waic.m1, site.rain=waic.m2, sitexrainxage=waic.m3, sitexage=waic.m4)
write.csv(d,'./output/waic.phi.csv')
