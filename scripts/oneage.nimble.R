#survival one site in nimble

#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(nimble)
library(here)

library(nimble)

data <- read.csv('./cleaned/ch.pmdi.csv', header = T)

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



##############################################################################
# 1. MODEL CODE
###############################################################################

code <- nimbleCode({
  
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
      
      logit(phi[i,t-1]) <- alpha +
        beta1[ ageclass[i, t-1] ] +
        beta3 * pmdi[i, t-1] +
        eps1[ year[i] ]
      
      # observation model
      ch[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p * z[i,t]
      
    } #t
  } #i
  
  # derived age-specific phi
  for(i in 1:15){
    phi_age[i] <- exp(alpha + beta1[i]) /
      (1 + exp(alpha + beta1[i]))
  }
  
})


##############################################################################
# 2. CONSTANTS AND DATA
###############################################################################

constants <- list(
  nind = nrow(ch),
  f = f,
  h = h,
  ageclass = ageclass,
  pmdi = pmdi.spring.sc,
  year = capyear
)

data_list <- list(
  ch = ch
)

##############################################################################
# 3. INITIAL VALUES
###############################################################################

# latent z initialization
z.init <- matrix(0, nrow = nrow(ch), ncol = ncol(ch))
for(i in 1:nrow(ch)){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i, f[i]] <- NA   # allow sampler to fill this one
}

inits <- list(
  alpha  = rnorm(1, 0, 1),
  p    = 0.5,
  z    = z.init,
  beta1 = c(0, rnorm(14, 0, 1)),
  beta3 = rnorm(1, 0, 1),
  eps1  = c(0, rnorm(13, 0, 1)),
  sigma = runif(1, 0.5, 5)
)

###############################################################################
# 4. BUILD & COMPILE MODEL
###############################################################################

Rmodel <- nimbleModel(code,
                      constants = constants,
                      data = data_list,
                      inits = inits)

Cmodel <- compileNimble(Rmodel)


###############################################################################
# 5. CONFIGURE MCMC
###############################################################################

conf <- configureMCMC(Rmodel,
                      monitors = c('alpha', #intercept
                                   'p', # recapture probability
                                   'beta1', # age effect
                                   'beta3', # pmdi effect
                                   'eps1', # random effect of cap year
                                   'sigma', # hyper parameter
                                   'phi_age' #derived parameter for age probability
                                   ))


###############################################################################
# 6. RUN MCMC
###############################################################################

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


###############################################################################
# 7. RETURN RESULTS
###############################################################################

samples <- runMCMC(Cmcmc,
                   niter = 10000,
                   nburnin = 5000,
                   thin = 10,
                   nchains = 3,
                   samplesAsCodaMCMC = TRUE
)

  

###############################################################################
# 8.  DIAGNOSTIC PLOTS
###############################################################################

plot(samples)

# convert CODA to matrix
post_mat<- as.matrix(samples)
     
  


