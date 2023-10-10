library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#load data 
data<-read.csv('./cleaned/caphx2022_nofawns.csv', header = T)

#rename column names to year
data<-data %>% rename('2008' = X2008, '2009' = X2009, '2010' = X2010, '2011' = X2011, '2012' = X2012, 
                      '2013' = X2013, '2014' = X2014, '2015' = X2015, '2016' = X2016, '2017' = X2017, 
                      '2018' = X2018, '2019' = X2019, '2020' = X2020, '2021' = X2021, '2022' = X2022)


CH<-data
known.fate<-data  #known deaths marked with 2
known.fate <- known.fate[, -c(1,17,18)]
known.fate<-as.matrix(known.fate) # 

#dataframe with only capture histories
CH <- CH[, -c(1,17,18)]
CH<-as.matrix(CH) # this will become y (i.e. our observations)

#create capture history with just 1s and 0s
indices <- which(CH == 2, arr.ind = TRUE)
CH[indices] <- 1

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$bs))
unique(bs) 
bs  #dmp is 1, e yana 2, w yana 3

#make age matrix, age class 1 is a 1.5 year old, age class 2 is a 2.5 year old 
ageclass<-data
for (i in 1:dim(ageclass)[1]){
  ageclass[i,2] <- 2007 - data$birth_year[i]  +1
  ageclass[i,3] <- 2007 - data$birth_year[i]  +2
  ageclass[i,4] <- 2007 - data$birth_year[i]  +3
  ageclass[i,5] <- 2007 - data$birth_year[i]  +4
  ageclass[i,6] <- 2007 - data$birth_year[i] +5
  ageclass[i,7] <- 2007 - data$birth_year[i] +6
  ageclass[i,8] <- 2007 - data$birth_year[i] +7
  ageclass[i,9] <- 2007 - data$birth_year[i] +8
  ageclass[i,10] <- 2007 - data$birth_year[i] +9
  ageclass[i,11] <- 2007 - data$birth_year[i] +10
  ageclass[i,12] <- 2007 - data$birth_year[i] +11
  ageclass[i,13] <- 2007 - data$birth_year[i] +12
  ageclass[i,14] <- 2007 - data$birth_year[i] +13
  ageclass[i,15] <- 2007 - data$birth_year[i] +14
  ageclass[i,16] <- 2007 - data$birth_year[i] +15
}

ageclass<- ageclass[,-c(1,17, 18)]

#Make any zeros or negatives NA
ageclass[ageclass<=0] <- NA
ageclass<-as.matrix(ageclass) 

# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0))
f <- apply(ageclass, 1, get.first) 
f #could use this to evaluate cohort effects by year

#create vector with last occasion for each individual, marked by 2, 15 for end of study 
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 15)
h

#create capture year vector
#make age matrix with 1 as fawn, 2 as immature, 3 as mature
capyear<-f
capyear
ageclass


# Specify model in JAGS language
sink("cjs-age-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
age.beta[1] <- 0
site.beta[1] <- 0

int ~ dnorm(0, 0.01)

for (u in 2:14){     #15 age classes but 14 survival periods                        
   age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
}

for (u in 2:3){ #3 sites
  site.beta[u] ~ dnorm(0, 0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1 #kept this in because we do know that for them to be kept in this dataset, they were alive at cap occ 2, for example row 142,
                        #the ind was born in 2013 and tagged then but was not recaught until 2020. We know that they were alive at cap occ 2

      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + site.beta[bs[i]] #+ eps.capyear[capyear[t-1]]int + 

        # Observation process
            CH[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i


#derived parameter
  for (i in 1:3){ #site 
    for (j in 1:14){ #age class
      survival[i,j] <- exp(int+ age.beta[j] + site.beta[i])/ (1 + exp(int+ age.beta[j] + site.beta[i]))
    }                     #delta method to convert from logit back to probability Powell et al. 2007
}
}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(CH), ncol = ncol(CH))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}


# Bundle data
jags.data <- list(h = h, CH = CH, f = f, nind = nrow(CH),  ageclass = ageclass, bs = bs)# ,h = h, n.occasions = ncol(CH),capyear=capyear, bs = bs


# Initial values
inits <- function(){list(int = rnorm(1, 0, 1), z = z.init, age.beta = c(NA, rnorm(13,0,1)), site.beta = c(NA, rnorm(2, 0, 1)))}
#,int = rnorm(1,0,1), ,  site.beta = c(NA, rnorm(2,0,1)),eps.capyear = rnorm(15,0,1), sigma.capyear = runif(1,0,10)

parameters <- c('int', 'age.beta', 'site.beta', 'p', 'survival')
# 'survival', 'site_diff' 'survival',, 'site_diff','eps.capyear' 'int','site.beta',

# MCMC settings
ni <- 4000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age.site <- jagsUI(jags.data, inits, parameters, "cjs-age-site.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
MCMCtrace(cjs.age.site)
print(cjs.age.site)