#model selection procedure

# ---- Load Data ----
#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

# data <- read.csv('./cleaned/ch.pmdi.csv', header = T)
data <- read.csv('./cleaned/ch.pmdi.dens.csv', header = T)

# data <- data %>% filter(!birth_year < 2011 & !year < 2011) #filter for years where we have all three sites
data <- data %>%  mutate(bs = recode(bs, ey = "control")) #rename ey to control to serve as a reference class 
data <- data %>%  mutate(bs = recode(bs, dmp = "wy")) #rename dmp to WY to merge wy and dmp into one bs 



# ---- JAGS Data Input ----
#how many capture histories do we have
count <- sum(data$status == '1')
#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'year', values_from = 'status', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

# known.fate<- known.fate[18:19,]
# ch <- ch[18:19,]
# ch[2,2] <- 1
# ch 
# known.fate <- ch

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



# ---- Model1: phi ~ site + age + weight + spring pmdi + density ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
  }

  beta2[1] <- 0 # birth site, control
  beta2[2] ~ dnorm(0,0.001) #birthsite, west yana

  beta3 ~ dnorm(0,0.001)  #capture year spring pmdi

  beta4 ~ dnorm(0, 0.001) #density

  beta5 ~ dlnorm(0, 0.1) #weight

  beta6[1] <- 0               #age x morpho
  for ( u in 2:15) {
    beta6[u] ~ dnorm(0, 0.1)
  }

  for (u in 1:nind){      #prior for missing morphometrics
    for (j in 1:occasions[u]){
    morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
       }
  }


  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + beta2[bs[i]]            #birth site, 2 sites
                                      + beta3*pmdi[i, t-1]   #capture year pmdi spring
                                      + beta4*density[i,t-1] #population density
                                      + beta5*morpho[i,t-1]   #morphometric measurement
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
  beta3 = rnorm(1,0,1),             #spring pmdi
  beta4 =  rnorm(1,0,1),            # density
  beta5 =  rlnorm(1,0,1),            #morpho
  beta6 = c(NA, rnorm(14,0,1)),
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'eps1')

# MCMC settings
ni <- 20000
nt <- 10
nb <- 15000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model1.csv') #DIC 1787



# # ---- Model1a: phi ~ age + spring pmdi + age*density ----
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("phi.age.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#   int ~ dnorm(0, 0.001)
# 
#   beta1[1] <- 0 #age
#   for ( u in 2:15) {
#     beta1[u] ~ dnorm(0, 0.01)  #age
#   }
# 
#   # beta2[1] <- 0 # birth site, control
#   # beta2[2] ~ dnorm(0,0.001) #birthsite, west yana
# 
#   beta3 ~ dnorm(0,0.001)  #capture year spring pmdi
# 
#   beta4 ~ dnorm(0, 0.001) #density
# 
#   # beta5 ~ dnorm(0, 0.001) #weight
# 
#   beta6[1] <- 0  #interaction age x density
#   for(u in 2:15){
#     beta6[u] ~ dnorm( 0 , 0.001 )
#   }
# 
#   for (u in 1:nind){      #prior for missing morphometrics
#     for (j in 1:occasions[u]){
#     morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
#        }
#   }
# 
# 
#   eps1[1] <- 0 #capture year RE
#    for (u in 2:14){  #prior for year effect
#     eps1[u] ~ dnorm(0,tau)
#   }
# 
# 
#   tau <- 1/(sigma*sigma)
#   sigma ~ dunif(0,100)
# 
# 
# # Likelihood
# for (i in 1:nind){
#    # Define latent state at first capture, we know for sure the animal is alive
#       z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):h[i]){
#         # State process
#             z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
#             logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
#                                       #+ beta2[bs[i]]            #birth site, 2 sites
#                                       + beta3*pmdi[i, t-1]   #capture year pmdi spring
#                                       + beta4*density[i,t-1] #population density
#                                       #+ beta5*morpho[i,t-1]   #morphometric measurement
#                                       + beta6[ageclass[i,t-1]]*density[i, t-1]
#                                       + eps1[year[i]]           #capture year random effect
#             
#            
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
# 
# 
# 
#       } #t
#    } #i
# 
#    #derived parameters
# 
#         # # for (i in 1:2){ #birthsite
#         #   for (j in 1:12) { #age
#         #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
#         #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
#         #
#         #     }
# 
#         # for (i in 1:100){ #density sim
#         #   for (j in 1:2) { #site
#         #     phi.dens[i, j] <- exp( int + beta3[j] 
#         # + beta4*density.sim[i] 
#         # + beta5[j]*density.sim[i] ) / #
#         #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
#         #
#         #     }}
#         # for (k in 1:14) { #cap year effect eps1
#         #
#         #   phi.year[k] <- exp( int + eps1[k]   )/
#         #                     (1 + exp( int + eps1[k] ))
#         #
#         #     }
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #Function for latent state
# z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }
# 
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
#                   bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
#                   NA_indices = NA_indices_weight, occasions = occasions_weight,
#                   morpho = weight, year = capyear, density = density, density.sim = density.sim)
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1),
#   z = z.init,
#   morpho = weight.init,
#   beta1 = c(NA, rnorm(14,0,1)),     #age beta
#   beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
#   beta3 = rnorm(1,0,1),             #spring pmdi
#   beta4 =  rnorm(1,0,1),            # density
#   beta5 =  rnorm(1,0,1),            #morpho
#   beta6 = c(NA, rnorm(14, 0, 1)),
#   eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
# )
# }
# 
# 
# parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'eps1')
# 
# # MCMC settings
# ni <- 10000
# nt <- 10
# nb <- 5000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
#                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.age)
# MCMCtrace(phi.age)
# 
# write.csv(phi.age$summary, './output/model1a.csv')
# 
# 
# # ---- Model1b: phi ~ age + springpmdi*density ----
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("phi.age.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#   int ~ dnorm(0, 0.001)
# 
#   beta1[1] <- 0 #age
#   for ( u in 2:15) {
#     beta1[u] ~ dnorm(0, 0.01)  #age
#   }
# 
#   # beta2[1] <- 0 # birth site, control
#   # beta2[2] ~ dnorm(0,0.001) #birthsite, west yana
# 
#   beta3 ~ dnorm(0,0.001)  #capture year spring pmdi
# 
#   beta4 ~ dnorm(0, 0.001) #density
# 
#   # beta5 ~ dnorm(0, 0.001) #weight
# 
#   beta6 ~ dnorm( 0 , 0.001 ) #interaction pmdi x density
#   
#   for (u in 1:nind){      #prior for missing morphometrics
#     for (j in 1:occasions[u]){
#     morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
#        }
#   }
# 
# 
#   eps1[1] <- 0 #capture year RE
#    for (u in 2:14){  #prior for year effect
#     eps1[u] ~ dnorm(0,tau)
#   }
# 
# 
#   tau <- 1/(sigma*sigma)
#   sigma ~ dunif(0,100)
# 
# 
# # Likelihood
# for (i in 1:nind){
#    # Define latent state at first capture, we know for sure the animal is alive
#       z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):h[i]){
#         # State process
#             z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
#             logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
#                                       #+ beta2[bs[i]]            #birth site, 2 sites
#                                       + beta3*pmdi[i, t-1]   #capture year pmdi spring
#                                       + beta4*density[i,t-1] #population density
#                                       #+ beta5*morpho[i,t-1]   #morphometric measurement
#                                       + beta6*pmdi[i,t-1]*density[i, t-1]
#                                       + eps1[year[i]]           #capture year random effect
# 
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
# 
# 
# 
#       } #t
#    } #i
# 
#    #derived parameters
# 
#         # # for (i in 1:2){ #birthsite
#         #   for (j in 1:12) { #age
#         #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
#         #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
#         #
#         #     }
# 
#         # for (i in 1:100){ #density sim
#         #   for (j in 1:2) { #site
#         #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
#         #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
#         #
#         #     }}
#         # for (k in 1:14) { #cap year effect eps1
#         #
#         #   phi.year[k] <- exp( int + eps1[k]   )/
#         #                     (1 + exp( int + eps1[k] ))
#         #
#         #     }
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #Function for latent state
# z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }
# 
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
#                   bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
#                   NA_indices = NA_indices_weight, occasions = occasions_weight,
#                   morpho = weight, year = capyear, density = density, density.sim = density.sim)
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1),
#   z = z.init,
#   morpho = weight.init,
#   beta1 = c(NA, rnorm(14,0,1)),     #age beta
#   beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
#   beta3 = rnorm(1,0,1),             #spring pmdi
#   beta4 =  rnorm(1,0,1),            # density
#   beta5 =  rnorm(1,0,1),            #morpho
#   beta6 = rnorm(1,0,1),            
#   eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
# )
# }
# 
# 
# parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'eps1')
# 
# # MCMC settings
# ni <- 10000
# nt <- 10
# nb <- 5000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
#                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.age)
# MCMCtrace(phi.age)
# 
# write.csv(phi.age$summary, './output/model1a.csv')
# 
# 
# # ---- Model2: phi ~ site + age  + spring pmdi + density ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
  }

  beta2[1] <- 0 # birth site, control
  beta2[2] ~ dnorm(0,0.001) #birthsite, west yana

  beta3 ~ dnorm(0,0.001)  #capture year spring pmdi

  beta4 ~ dnorm(0, 0.001) #density




  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + beta2[bs[i]]            #birth site, 2 sites
                                      + beta3*pmdi[i, t-1]   #capture year pmdi spring
                                      + beta4*density[i,t-1] #population density
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
  beta3 = rnorm(1,0,1),             #spring pmdi
  beta4 =  rnorm(1,0,1),            # density
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model2.csv')

# 
# ---- Model3: phi ~ site + age  + density ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
  }

  beta2[1] <- 0 # birth site, control
  beta2[2] ~ dnorm(0,0.001) #birthsite, west yana

  beta4 ~ dnorm(0, 0.001) #density




  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + beta2[bs[i]]            #birth site, 2 sites
                                      + beta4*density[i,t-1] #population density
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
  beta4 =  rnorm(1,0,1),            # density
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1','beta2', 'beta4',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model3.csv')
# 
# 
# ---- Model4: phi ~ site + age  ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
  }

  beta2[1] <- 0 # birth site, control
  beta2[2] ~ dnorm(0,0.001) #birthsite, west yana





  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + beta2[bs[i]]            #birth site, 2 sites
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1','beta2',  'eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model4.csv')


# ---- Model5: phi ~  age + capyear RE TOP MODEL ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
  }


  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        for (k in 1:14) { #cap year effect eps1

          phi.year[k] <- exp( int + eps1[k]   )/
                            (1 + exp( int + eps1[k] ))

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1',  'eps1', 'phi.year')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
# 
write.csv(phi.age$summary, './output/model5.csv')
# 
# 
# ---- Model6: phi ~  age no random effect of year ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.01)  #age
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        for (k in 1:14) { #cap year effect eps1

          phi.year[k] <- exp( int + eps1[k]   )/
                            (1 + exp( int + eps1[k] ))

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(14,0,1))     #age beta
)
}


parameters <- c('int', 'beta1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model6.csv')


# 
#
# ---- Model7: intercept only ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)




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
            logit(phi[i,t-1]) <-  int

            z.new[i,t] ~ dbern(mu1[i,t])
            res[i,t-1] <- z[i,t] - mu1[i,t]
            res.new[i,t-1] <- z.new[i,t] - mu1[i,t]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]


      } #t
   } #i

   #derived parameters


        # # for (i in 1:2){ #birthsite
        #   for (j in 1:12) { #age
        #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
        #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
        #
        #     }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init
)
}


parameters <- c('int', 'fit', 'fit.new')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)

write.csv(phi.age$summary, './output/model7.csv')


