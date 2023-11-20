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

# The difference in wAIC tells us that the covariate has some effect
# wAIC of m1 the model with covariate << wAIC of m0 intercept only
data.frame(m1 = waic.m1)