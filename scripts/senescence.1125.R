#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

data<- read.csv('./cleaned/final.ch1.csv', header = T)

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
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass matrix with continuous age cov
data$age.sc <- scale(data$ageclass)
age.sc<- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc <- as.matrix(age.sc[,-1])


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
h <- replace(h, is.infinite(h), 15)
h
f-h #check for zero

#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'animal_id' )
weight<- as.matrix(weight[,-1])

antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcsin', id_cols = 'animal_id' )
antlers<- as.matrix(antlers[,-1])

#rainfall
cy.rain<-pivot_wider(data, names_from = 'year', values_from = 'annual.sc', id_cols = 'animal_id' )
cy.rain<-as.matrix(cy.rain[,-1])

by.rain<-pivot_wider(data, names_from = 'year', values_from = 'by.rain.sc', id_cols = 'animal_id' )
by.rain<-as.matrix(by.rain[,-1])

##create simulated birthyear rainfall vector
nvalues <- 100
by.rain.sim <- seq(from = min(by.rain, na.rm = T), to = max(by.rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data



#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 


#survival as a function of birth year rain, birth site, age, and NO interaction bt bs and birthyear rain
# Specify model in JAGS language, top performing model
set.seed(100)
sink("cjs.site.rain.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

int ~ dnorm(0,0.01)
age.beta[1] <- 0
bs.beta[1] <- 0
age.bs.beta[1] <- 0
eps.capyear[1] <- 0

for (u in 2:15){
  age.beta[u] ~ dnorm(0,0.01)
}

for (u in 2:3){
  bs.beta[u] ~ dnorm(0, 0.01)
}

for (u in 2:15){
  age.bs.beta[u] ~ dnorm(0, 0.01) 
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

for (u in 2:15){
  eps.capyear[u] ~ dnorm(0, tau.capyear)
}
  tau.capyear <- 1/(sigma.capyear*sigma.capyear)
  sigma.capyear ~ dunif(0,100)

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int  + age.beta[ageclass[i,t-1]] 
                                      + bs.beta[bs[i]] + age.bs.beta[ageclass[i,t-1]]*bs[i]
                                      + eps.capyear[capyear[i]]
                                      

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
    for (j in 1:10){ #ageclass
      for (k in 1:3){
      survival[j,k] <- exp(int + age.beta[j] + bs.beta[k] + age.bs.beta[j]*k)/
                            (1 + exp(int + age.beta[j] + bs.beta[k] + age.bs.beta[j]*k))

      }
    }


}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}


# Bundle data
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, bs = bs, capyear = capyear )
# by.rain = by.rain, capyear = capyear,  bs = bs, by.rain.sim = by.rain.sim)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, 
                         age.beta = c(NA, rnorm(14,0,1)), 
                          bs.beta = c(NA, rnorm(2,0,1)),
                          age.bs.beta = c(NA, rnorm(14,0,1)),
                            eps.capyear = c(NA, runif(14, 0, 100)))}

# site.beta = c(NA, rnorm(2,0,1)), by.rain.beta = rnorm(1,0,1),
# eps.capyear = c(NA, rnorm(13,0,1)))} #age.beta = c(NA, rnorm(13,0,1)), , 'survival_diff',

parameters <- c('int', 'age.beta', 'bs.beta', 'age.bs.beta', 'survival')#, , 'site.beta','by.rain.beta', 'p', 'eps.capyear', 'survival' )

# MCMC settings
ni <- 15000
nt <- 1
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.site.rain <- jagsUI(jags.data, inits, parameters, "cjs.site.rain.jags", n.chains = nc,
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(cjs.site.rain)
MCMCtrace(cjs.site.rain)
