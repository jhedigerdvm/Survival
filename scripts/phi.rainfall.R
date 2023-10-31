#Survival model for age and birthsite, capyear as a random effect
#treatment and control survival are not different, but TGT survive less 
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

data<- read.csv('./cleaned/caphx.rainfall.long.csv', header = T)


ch <- matrix(data$status,nrow = 506, ncol = 15) #may need to create a matrix of CH instead of a vector 

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates 
ch[indices] <- 1

#create birthsite vector
bs <- as.numeric(factor(data$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass vector
ageclass <- data$ageclass

#create animal id vector
id <- as.numeric(factor(data$animal_id))

#create birth year vector
birthyear <- as.numeric(as.factor(data$birth_year))

#create capture year vector
capyear <- as.numeric(as.factor(data$year))

#scale and center covariate annual rainfall
annual.rain <- scale(data$annual)
annual.rain<- annual.rain[,1]


#create vector with first capture occasion which is equal to their birthyear
f <- birthyear

#create vector with last occasion for each individual, marked by 2, 15 for end of study 
#rework h to only include capture myopathy or harvest, do not censor natural mortality 
# get.last<- function(x) min(which(x>1))
# h <- apply(known.fate,1,get.last)
# h <- replace(h, is.infinite(h), 15)
# h


# Specify model in JAGS language
sink("cjs.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
phi ~ dbeta(1,1)


# Likelihood 
for (i in 1:7590){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):7590){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2

        # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i


}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(data), ncol = 15)
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[f[i]:15] <- 1
#   z.init[f[i]] <- NA
# }


# Bundle data
jags.data <- list(ch = ch, f = f)#, 

# Initial values
inits <- function(){list(z = z.init, phi = rbeta(1,1,1), p = rbeta(1,1,1))} #

parameters <- c('phi','p')
# 'survival', 'site_diff' 'survival',, 'site_diff','eps.capyear' 'int','site.beta',

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age.site <- jagsUI(jags.data, inits, parameters, "cjs.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
