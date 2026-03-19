#model selection for quadratic relationship of age
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
data <- read.csv('./cleaned/ch.pmdi.dens.csv', header = T)

# data <- data %>% filter(!birth_year < 2011 & !year < 2011) #filter for years where we have all three sites
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

sum(bs==1) #control
sum(bs==2) #tgt + dmp
unique(bs)

#create ageclass matrix treating age as categorical
data$age.sc <- scale(data$ageclass)
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

age.sc <- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc<- as.matrix(age.sc[,-1])

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

#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'animal_id' )
weight<- as.matrix(weight[,-1])
weight <- scale(weight) # scale and center


#function for weight matrix
weight.init <- weight
weight.init[is.na(weight.init)]<-0 #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
weight.init[!is.na(weight)]<-NA
for (i in 1:dim(weight.init)){ #cant have mean weight for years before animal was first captured
  weight.init[i,1:f[i]]<-NA
}

#need to find where NA values are in weight matrix, will use this information to build priors
weight <- as.data.frame(weight)
indices <- as.data.frame(which(is.na(weight), arr.ind=T))
indices <- indices %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
NA_indices_weight <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
for(i in 1:nrow(indices)){
  NA_indices_weight[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
}
weight<-as.matrix(weight)

#how many occasions does each individual have of an NA weight
occasions_weight <- rowSums(is.na(weight))


#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 


##now do the same thing we did with weight with antlers
antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcsin', id_cols = 'animal_id' )
antlers<- as.matrix(antlers[,-1])

antlers <- scale(antlers) # scale and center


#function for weight matrix
antlers.init <- antlers
antlers.init[is.na(antlers.init)]<-0 #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
antlers.init[!is.na(antlers)]<-NA
for (i in 1:dim(antlers.init)){ #cant have mean weight for years before animal was first captured
  antlers.init[i,1:f[i]]<-NA
}

#need to find where NA values are in weight matrix, will use this information to build priors
antlers <- as.data.frame(antlers)
indices_antlers <- as.data.frame(which(is.na(antlers), arr.ind=T))
indices_antlers <- indices_antlers %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
NA_indices_antlers <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
for(i in 1:nrow(indices_antlers)){
  NA_indices_antlers[indices_antlers[[i,1]],indices_antlers[[i,3]]] <- indices_antlers[[i,2]]
}
antlers<-as.matrix(antlers)

#how many occasions does each individual have of an NA weight
occasions_antlers <- rowSums(is.na(antlers))
#


#add simulated antler values
nvalues <- 100
antler.sim <- seq(from = min(antlers, na.rm = T), to = max(antlers, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data


#create spring pmdi data
pmdi.spring.sc<-pivot_wider(data, names_from = 'year', values_from = 'pmdi_spring.sc', id_cols = 'animal_id' )
pmdi.spring.sc<-as.matrix(pmdi.spring.sc[,-1])

##create simulated birthyear rainfall vector
nvalues <- 100
pmdi.spring.sc.sim <- seq(from = min(pmdi.spring.sc, na.rm = T), to = max(pmdi.spring.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

unique(bs) # two sites

#add density data
data$densitykm2.sc <- scale(data$densitykm2)
density <-pivot_wider(data, names_from = 'year', values_from = 'densitykm2.sc', id_cols = 'animal_id' )
density <- density[,-1]
density <- as.matrix(density)

##create simulated densities
nvalues <- 100
density.sim <- seq(from = min(density, na.rm = T), to = max(density, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

#age simulation for continuous model
nvalues <- 15
age.sim <- seq(from = min(age.sc, na.rm = T), to = max(age.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 




#---- Model1: full model ----


# Specify model in JAGS language
set.seed(100)
sink("model1.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  
  beta3[1] <- 0 #ey
  beta3[2] ~ dnorm(0, 0.01)  #wy
  
  beta4 ~ dnorm( 0 , 0.001) #density
  
  beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
                                  + beta3[ bs[ i ]] #birthsite
                                  + beta4*density[ i , t-1 ]
                                  + beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(1,0,1)),#birth site
  beta4 = rnorm(1,0,1), #density
  beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model1<- jagsUI(jags.data, inits, parameters, "model1.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(model1)


#---- Model2: phi ~ int+ age + age2 + site + density ----


# Specify model in JAGS language
set.seed(100)
sink("model2.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  
  beta3[1] <- 0 #ey
  beta3[2] ~ dnorm(0, 0.01)  #wy
  
  beta4 ~ dnorm( 0 , 0.001) #density
  
  #beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
                                  + beta3[ bs[ i ]] #birthsite
                                  + beta4*density[ i , t-1 ]
                                  #+ beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(1,0,1)),#birth site
  beta4 = rnorm(1,0,1), #density
  #beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model2<- jagsUI(jags.data, inits, parameters, "model2.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model2)


#---- Model3: phi ~ int+ age + age2 + site  ----


# Specify model in JAGS language
set.seed(100)
sink("model3.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  
  beta3[1] <- 0 #ey
  beta3[2] ~ dnorm(0, 0.01)  #wy
  
  #beta4 ~ dnorm( 0 , 0.001) #density
  
  #beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
                                  + beta3[ bs[ i ]] #birthsite
                                  #+ beta4*density[ i , t-1 ]
                                  #+ beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(1,0,1)),#birth site
  beta4 = rnorm(1,0,1), #density
  #beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model3<- jagsUI(jags.data, inits, parameters, "model3.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model3)

#---- Model4: int + age + age2 ----

# Specify model in JAGS language
set.seed(100)
sink("model4.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  
  # beta3[1] <- 0 #ey
  # beta3[2] ~ dnorm(0, 0.01)  #wy
  
  #beta4 ~ dnorm( 0 , 0.001) #density
  
  #beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
                                  #+ beta3[ bs[ i ]] #birthsite
                                  #+ beta4*density[ i , t-1 ]
                                  #+ beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  # beta3 = c(NA, rnorm(1,0,1)),#birth site
  # beta4 = rnorm(1,0,1), #density
  #beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model4<- jagsUI(jags.data, inits, parameters, "model4.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model4)

#---- Model5: phi ~ int   ----


# Specify model in JAGS language
set.seed(100)
sink("model5.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  # beta1 ~ dnorm(0, 0.001)
  # beta2 ~ dnorm(0,0.001)
  
  # beta3[1] <- 0 #ey
  # beta3[2] ~ dnorm(0, 0.01)  #wy
  
  #beta4 ~ dnorm( 0 , 0.001) #density
  
  #beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
            logit(phi[i,t-1]) <- int #+ beta1*ageclass[ i , t-1 ]   #age 
                                  #+ (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
                                  #+ beta3[ bs[ i ]] #birthsite
                                  #+ beta4*density[ i , t-1 ]
                                  #+ beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]

          # Log-likelihood for WAIC      
            loglik[i,t] <- logdensity.bern(ch[i,t], mu2[i,t])


      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  # beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  # beta2=rnorm(1,0,1),
  # beta3 = c(NA, rnorm(1,0,1)),#birth site
  # beta4 = rnorm(1,0,1), #density
  #beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1', 'loglik')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model5<- jagsUI(jags.data, inits, parameters, "model5.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model5)

# loglik_array <- model5$sims.list$loglik
# 
# library(loo)
# 
# # flatten i,t into columns
# loglik_mat <- apply(loglik_array, 1, c)
# loglik_mat <- t(loglik_mat)
# 
# # remove columns with NA (unobserved entries)
# loglik_mat <- loglik_mat[, colSums(is.na(loglik_mat)) == 0]
# loglik_mat <- loglik_mat[, apply(loglik_mat, 2, var) > 0]
# loo_result <- loo(loglik_mat)
# 
# print(loo_result)
# 
# print(loo_result)
# loo_result$diagnostics$pareto_k
# loo_result <- loo(loglik_mat, moment_match = TRUE)

#---- Model6: phi ~ int+ age + age2 + density  ----


# Specify model in JAGS language
set.seed(100)
sink("model6.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  
  # beta3[1] <- 0 #ey
  # beta3[2] ~ dnorm(0, 0.01)  #wy
  
  beta4 ~ dnorm( 0 , 0.001) #density
  
  #beta5 ~ dnorm ( 0 , 0.001) #pmdi spring
  


  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
                                  #+ beta3[ bs[ i ]] #birthsite
                                  + beta4*density[ i , t-1 ]
                                  #+ beta5*pmdi[ i , t-1 ]
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
       
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(1,0,1)),#birth site
  beta4 = rnorm(1,0,1), #density
  #beta5 = rnorm(1,0,1), # spring pmdi
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3',"beta4",'beta5',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model6<- jagsUI(jags.data, inits, parameters, "model6.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model6)



model_list <- list(
  m1 = c(model1$pD, model1$DIC),
  m2 = c(model2$pD, model2$DIC),
  m3 =  c(model3$pD, model3$DIC),
  m4 =  c(model4$pD, model4$DIC),
  m5 =  c(model5$pD, model5$DIC),
  m6 = c(model6$pD, model6$DIC)
)

model_mat <- do.call(rbind, model_list)

colnames(model_mat) <- c("pD", "DIC")

models<- as.data.frame(model_mat)
