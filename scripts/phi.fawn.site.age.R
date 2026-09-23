#models including fawn data
# ---- Load Packages ----
#final run of survival analyses
library(jagsUI)
library(ggplot2)
library(ggdist)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

# ---- Load data ----
# data <- read.csv('./cleaned/ch.pmdi.csv', header = T)
data <- read.csv('./cleaned/fawncaphx.csv', header = T)

#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'cap_year', values_from = 'status_cam', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #83 individuals with known fates 
ch[indices] <- 1


# Create vector with the occasion each indiv is marked
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

#create birthsite vector
id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = control, 2  = dmp + wy

sum(bs==1) #control
sum(bs==2) #tgt + dmp
unique(bs)

#create ageclass matrix treating age as categorical
  # ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
  # ageclass<- ageclass[,-1]
  # ageclass<-as.matrix(ageclass)

#scale and center age
  age.mean <- mean(data$age, na.rm = TRUE)
  age.sd <- sd(data$age, na.rm = TRUE)
  data$age.sc <- (data$age - age.mean) / age.sd
  
  age.sc <- pivot_wider(data, names_from = 'cap_year', values_from = 'age.sc', id_cols = 'animal_id' )
  age.sc<- as.matrix(age.sc[,-1])
  
  age.pred <- seq(0.5, 15.5, by = 1)
  
  age.pred.sc <- (age.pred - age.mean) / age.sd


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

#identify 0 and +1 indiv
  problem <- which((h - f) <= 1)
  
  length(problem)
  
  data.frame(
    id = id.bs.by$animal_id[problem],
    f = f[problem],
    h = h[problem]
  )

  known.fate[problem, ]
  
#identify the individuals to keep and adjust lists and vectors
  keep <- which((h - f) > 0)
  
  known.fate <- known.fate[keep, ]
  ch <- ch[keep, ]
  f <- f[keep]
  h <- h[keep]
  bs <- bs[keep]
  birthyear <- birthyear[keep]
  capyear <- capyear[keep]
  age.sc <- age.sc[keep, ]

#verify the fix
  which((h - f) <= 0)


  
  #age simulation for continuous model
  age.sim <- age.pred.sc
  nvalues <- length(age.sim)

# 
#---- Model1: phi ~ int+ age + age2 + site  ----


# Specify model in JAGS language
set.seed(100)
sink("model1.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)
  beta3[1] <- 0
  

  beta1 ~ dnorm(0, 0.001) #age
  beta2 ~ dnorm(0,0.001) #age ^2
  
  for ( u in 2:3){
    beta3[u] ~ dnorm(0, 0.001)
  }
  
 
  
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]   #age
                                  + (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
                                  + beta3[bs[i]]
                                

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
       for (j in 1:15) { #simulated age
        for ( i in 1:3){
            phi.age[i, j] <- exp( int + beta1*age.sim[j]
                                      + beta2*age.sim[j]*age.sim[j]
                                      + beta3[i]
                                     
                                      ) /
                               (1 + exp( int+ beta1*age.sim[j]
                                              + beta2*age.sim[j]*age.sim[j]
                                              + beta3[i]
                                              

                                             ) )

          }}


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc,
                  bs = bs, year = capyear, age.sim = age.sim)
                  
                  #morpho = weight, , density = density, density.sim = density.sim) NA_indices = NA_indices_weight, occasions = occasions_weight,morpho.sim = weight.sim, pmdi = pmdi.spring.sc,pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(2,0,1))#birth site
  
  
  
  # 
  # eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'phi.age')

# MCMC settings
ni <- 2000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model1<- jagsUI(jags.data, inits, parameters, "model1.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model1)





#---- Model2: phi ~ int+ age + age2 + site + site*age + site*age2  ----


# Specify model in JAGS language
set.seed(100)
sink("model2.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)
  beta3[1] <- 0
  beta4[1] <- 0
  beta5[1] <- 0

  beta1 ~ dnorm(0, 0.001) #age
  beta2 ~ dnorm(0,0.001) #age ^2
  
  for ( u in 2:3){
    beta3[u] ~ dnorm(0, 0.001)
  }
  
  for ( u in 2:3){
      beta4[u] ~ dnorm(0, 0.001)
  }
   
   for ( u in 2:3){
      beta5[u] ~ dnorm(0, 0.001)
   } 
  
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]   #age
                                  + (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
                                  + beta3[bs[i]]
                                  + beta4[bs[i]]*ageclass[i,t-1]
                                  + beta5[bs[i]]*(ageclass[i,t-1]*ageclass[i,t-1])

                                

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
       for (j in 1:15) { #simulated age
        for ( i in 1:3){
            phi.age[i, j] <- exp( int + beta1*age.sim[j]
                                      + beta2*age.sim[j]*age.sim[j]
                                      + beta3[i]
                                      + beta4[i]*age.sim[j]
                                      + beta5[i]*age.sim[j]*age.sim[j]
                                      ) /
                               (1 + exp( int+ beta1*age.sim[j]
                                              + beta2*age.sim[j]*age.sim[j]
                                              + beta3[i]
                                              + beta4[i]*age.sim[j]
                                              + beta5[i]*age.sim[j]*age.sim[j]

                                             ) )

          }}


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc,
                  bs = bs, year = capyear, age.sim = age.sim)

#morpho = weight, , density = density, density.sim = density.sim) NA_indices = NA_indices_weight, occasions = occasions_weight,morpho.sim = weight.sim, pmdi = pmdi.spring.sc,pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(2,0,1)),#birth site
  beta4 = c(NA, rnorm(2,0,1)),#,#birth site
  beta5 = c(NA, rnorm(2,0,1))#,#birth site
  
  
  # 
  # eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'phi.age')

# MCMC settings
ni <- 2000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model2<- jagsUI(jags.data, inits, parameters, "model2.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model2)

