#survival one site in nimble

#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
#library(mcmcr) 
library(viridis)
library(nimble)
library(here)
library(coda)

data <- read.csv('ch.pmdi.csv', header = T)

data <- data %>%  mutate(bs = recode(bs, ey = "control")) #rename ey to control to serve as a reference class 
data <- data %>%  mutate(bs = recode(bs, dmp = "wy")) #rename dmp to WY to merge wy and dmp into one bs 

#how many capture histories do we have
count <- sum(data$status == '1')
#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'year', values_from = 'status', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates 
ch[indices] <- 1

dim(ch)
#ch_new <- matrix(NA, nrow = 489, ncol = 15)

# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

#create birthsite vector
id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = control, 2  = dmp + wy
# bs.2 <- bs
# bs.2[bs.2 %in% c(3)] <-2 #binary variable of dmp versus pasture
# sum(bs.2==1) #182 individuals in DMP
# sum(bs.2==2) #307 individuals in pasture

sum(bs==1) #control
sum(bs==2) #tgt + dmp
unique(bs)
# 
# #create ageclass matrix with continuous age cov
# data$age.sc <- scale(data$ageclass)
# age.sc<- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
# age.sc <- as.matrix(age.sc[,-1])


#create ageclass matrix treating age as categorical
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

#create birth year vector
birthyear <- as.numeric(as.factor(id.bs.by$birth_year))

#create capture year vector
capyear <- f

# 
# create vector with last occasion for each individual, marked by 2, 15 for end of study
# rework h to only include capture myopathy or harvest, do not censor natural mortality
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 15) #change to equal number of columns/years
h
f-h #check for zero


#create spring pmdi data
pmdi.spring.sc<-pivot_wider(data, names_from = 'year', values_from = 'pmdi_spring.sc', id_cols = 'animal_id' )
pmdi.spring.sc<-as.matrix(pmdi.spring.sc[,-1])

##create simulated birthyear rainfall vector
nvalues <- 100
pmdi.spring.sc.sim <- seq(from = min(pmdi.spring.sc, na.rm = T), 
                          to = max(pmdi.spring.sc, na.rm = T), 
                          length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

#create initial values for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
dim(z.init)
for(i in 1:nrow(ch)){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i, 1:f[i]] <- NA   # allow sampler to fill this one
}

#create initial values for simulated capture history
ch.new.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
dim(ch.new.init)
for(i in 1:nrow(ch)){
  ch.new.init[i, f[i]:h[i]] <- 1
  ch.new.init[i, 1:f[i]] <- NA   # allow sampler to fill this one
}

#how many individuals were captured in the capture history after first occasion?
rowsum <- vector()
for (i in 1:nrow(ch)){
rowsum[i] <- sum(ch[i,(f[i]+1):h[i]])
}
N.obs <- sum(rowsum)

##############################################################################
# 1. MODEL CODE
###############################################################################

sink("wolf_jags.txt")
cat("
    model {
  
  # prior for recapture p
  p ~ dbeta(1, 1)
  
  # priors
    alpha ~ dnorm(0, 0.001)
    
    beta1[1] <- 0        # age reference
    eps1[1]  <- 0        # year RE reference
    
    beta3 ~ dnorm(0, 0.001)
  
  # age-class effects
    for(u in 2:15){
      beta1[u] ~ dnorm(0, 0.01)
    }
    
  # year RE
    for(u in 2:14){
      eps1[u] ~ dnorm(0, tau)
    }
    
  tau <- 1/(sigma * sigma)
  sigma ~ dunif(0, 100)
  
  
  # STATE + OBSERVATION PROCESS
  for(i in 1:nind){
    
    # latent state initialized at first capture
    z[i, f[i]] <- 1

    
    for(t in (f[i]+1):h[i]){
      
      # survival state process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      
      # z_rep[i,t] ~ dbern(mu1_rep[i,t])
      # mu1_rep[i,t] <- phi[i,t-1] * z_rep[i,t-1]
      
      logit(phi[i,t-1]) <- alpha +
                            beta1[ ageclass[i, t-1] ] +
                            beta3 * pmdi[i, t-1] +
                            eps1[ year[i] ]
      
      # observation model
      ch[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p * z[i,t]
      ch_new[i,t] ~ dbern(mu2[i,t])
      
    } #t
      indiv_sum[i] <- sum(ch_new[i, (f[i]+1):h[i]])
      
  } #i
        N_rep <- sum(indiv_sum[1:489])
        Bayesian_p <- step(N_rep - N.obs)

  #for the code below, i assume I need to assign values for N and first, I tried using with the indexing from my i,t above
  # for(i in 1:N) {
  #   for(t in first[i]:T) {
  #     N_obs<- N_obs + ch[i,t] #this produces a running count of number of individuals in that capture history
  #     N_rep <-N_rep + ch_new[i,t]
  #   }
  # }
  
  # # derived age-specific phi
  # for(i in 1:15){
  #   phi_age[i] <- exp(alpha + beta1[i]) /
  #     (1 + exp(alpha + beta1[i]))
  # }
  # 
  
    }",fill=TRUE)
sink()

win.data <- list(    ch = ch,
                   nind = nrow(ch),
                   f = f,
                   h = h,
                   ageclass = ageclass,
                   pmdi = pmdi.spring.sc,
                   year = capyear, N.obs = N.obs)


#provide initial values
inits <- function(){list(
  alpha  = rnorm(1, 0, 1),
  p    = 0.5,
  z    = z.init,
  ch.new = ch.new.init,
  beta1 = c(NA, rnorm(14, 0, 1)),
  beta3 = rnorm(1, 0, 1),
  eps1  = c(NA, rnorm(13, 0, 1)),
  sigma = runif(1, 0.5, 5))}

nb <- 500 #100000 #This will take about 20-30 minutes to run on a faster machine
ni <- 2000 #300000 
nt <- 10  #3
nc <- 3

# Parameters to estimate
params <- c(  #'alpha', #intercept
              #'beta1', # age effect
              #'beta3', # pmdi effect
              #'eps1', # random effect of cap year
              #'sigma', 'ch_new', # hyper parameter
              #'N_obs',
              'Bayesian_p',
              'N_rep', 
              'indiv_sum' )

out <- jags(win.data, inits, params, "wolf_jags.txt", n.adapt=1, n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt, parallel=T)
as.mcmc(out$summary)
