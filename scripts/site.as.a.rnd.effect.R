#phi as a function of antlers

#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

data<- read.csv('./cleaned/final.ch1.csv', header = T)
data <- data %>% filter(!birth_year < 2011 & !year < 2011) #filter for years where we have all three sites

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
bs.2 <- bs
bs.2[bs.2 %in% c(3)] <-2 #binary variable of dmp versus pasture
sum(bs.2==1) #182 individuals in DMP
sum(bs.2==2) #307 individuals in pasture

#create ageclass matrix with continuous age cov
data$age.sc <- scale(data$ageclass)
age.sc<- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc <- as.matrix(age.sc[,-1])


#create ageclass matrix treating age as categorical
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

#create age class of juvenile, mature, geriatric
age.bin <- ageclass
age.bin[age.bin %in% c(2,3)] <- 1 #juveniles 1.5 - 3.5 year old
age.bin[age.bin %in% c(4,5,6,7,8)] <- 2 #mature 4.5 - 8.5 year old
age.bin[age.bin %in% c(9,10,11,12,13,14,15)] <- 3 #geriatric 9.5 and older

#create birth year vector
birthyear <- as.numeric(as.factor(id.bs.by$birth_year))

#create capture year vector
capyear <- f

# 
# create vector with last occasion for each individual, marked by 2, 15 for end of study
# rework h to only include capture myopathy or harvest, do not censor natural mortality
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 12) #change to equal number of columns/years
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



#rainfall
cy.rain<-pivot_wider(data, names_from = 'year', values_from = 'annual.sc', id_cols = 'animal_id' )
cy.rain<-as.matrix(cy.rain[,-1])

by.rain<-pivot_wider(data, names_from = 'year', values_from = 'by.rain.sc', id_cols = 'animal_id' )
by.rain<-as.matrix(by.rain[,-1])

##create simulated birthyear rainfall vector
nvalues <- 100
by.rain.sim <- seq(from = min(by.rain, na.rm = T), to = max(by.rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data


##create simulated capture year rainfall vector
nvalues <- 100
cy.rain.sim <- seq(from = min(cy.rain, na.rm = T), to = max(cy.rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data


#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 


#add simulated antler values
nvalues <- 100
antler.sim <- seq(from = min(antlers, na.rm = T), to = max(antlers, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 

#add simulated age values
age.sim <- seq(from = min(age.sc, na.rm = T), to = max(age.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of age


# ---- Model1: phi ~ weight x site interaction ----

#survival as a function of weight (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.morpho.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

# int ~ dnorm(0,0.001)
# 
  beta1 ~ dnorm(0,0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous
  
for ( u in 1:3) { #interaction between morpho and site
  beta4[u] ~ dlnorm(0, 0.1)
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
}

for (u in 1:3){         #prior for random site effect
  eps1[u] ~ dlnorm(0,tau)
  }

for (u in 1:12){  #prior for year effect
  eps2[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <- eps1[bs[i]] + beta1*ageclass[i,t-1]
                                      + beta2*morpho[i, t-1]
                                      + beta3*by.rain[i, t-1]
                                      + beta4[bs[i]]*morpho[i,t-1]
                                      + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
            
          #predictive checks
            res[i,t] <- z[i,t]- mu1[i,t]
            z.new[i,t] ~ dbern(mu1[i,t], tau)
            res.new[i,t] ~ z.new[i,t] - mu1[i,t]
            
      } #t
   } #i

   #derived parameters
      for (j in 1:100 ) { #age simulation, beta1
      for (l in 1:3){ #site, eps1

      survival[j,l] <- exp( beta2*morpho.sim[j]  + eps1[l] + beta4[l]*morpho.sim[j] )/ 
                            (1 + exp(beta2*morpho.sim[j]  + eps1[l] + beta4[l]*morpho.sim[j]  )) 
  
    } # for j
    } # for l
        
        fit <- sum(res[])
        fit.new <- sum(res.new[])


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, by.rain = by.rain,
                   bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(#int = rnorm(1,0,1), 
                         z = z.init,
                         eps1 = rlnorm(3,0,1), #birth site random effect
                         morpho = weight.init, #initial values for NA morphos
                         beta1 = rnorm(1,0,1), #age beta
                         beta2 = rlnorm(1, 0, 1),#morpho beta
                         beta3 = rnorm(1, 0, 1), # rain beta
                         beta4 = rlnorm(3,0,1),#morpho and site interaction
                         eps2 = rnorm(12, 0, 1) #capture year random effect
                          )
                        }


parameters <- c('eps1', 'beta1','beta2', 'beta3', 'beta4', 'eps2',
                'fit', 'fit.new', 'survival')  

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.morpho <- jagsUI(jags.data, inits, parameters, "phi.morpho.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.morpho)
MCMCtrace(phi.morpho)
write.csv(phi.morpho$summary, './output/phi.weighttxsite.csv', row.names = T)


#create a tibble of the posterior draws
gather<- phi.morpho %>% gather_draws(survival[morpho, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$morpho == 2)[1] # 4500 values of antler 1
first_idx
#unscale and uncenter weight
morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
morpho.sim.usc1 <- for (i in morpho.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)

}

gather$bodymass <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=bodymass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Body mass (lbs)", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.weightxsite.jpg', phi.plot, width = 15, height = 10)


# ---- Model1a: phi ~ antlers x site interaction ----

#survival as a function of antlers (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.morpho.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

# int ~ dnorm(0,0.001)
# 
  beta1 ~ dnorm(0,0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous
  
for ( u in 1:3) { #interaction between morpho and site
  beta4[u] ~ dlnorm(0, 0.1)
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
}

for (u in 1:3){         #prior for random site effect
  eps1[u] ~ dlnorm(0,tau)
  }

for (u in 1:12){  #prior for year effect
  eps2[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <- eps1[bs[i]] + beta1*ageclass[i,t-1]
                                      + beta2*morpho[i, t-1]
                                      + beta3*by.rain[i, t-1]
                                      + beta4[bs[i]]*morpho[i,t-1]
                                      + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

   #derived parameters
      for (j in 1:100 ) { #antler simulation, beta2
      for (l in 1:3){ #site, eps1

      survival[j,l] <- exp( beta2*morpho.sim[j]  + eps1[l] + beta4[l]*morpho.sim[j] )/ 
                            (1 + exp(beta2*morpho.sim[j]  + eps1[l] + beta4[l]*morpho.sim[j]  )) 
  
    } # for j
    } # for l


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = antler.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_antlers, occasions = occasions_antlers,morpho = antlers, 
                  year = capyear)

# Initial values
inits <- function(){list(#int = rnorm(1,0,1), 
  z = z.init,
  eps1 = rlnorm(3,0,1), #birth site random effect
  morpho = antlers.init, #initial values for NA morphos
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = rlnorm(3,0,1),#morpho and site interaction
  eps2 = rnorm(12, 0, 1)#capture year random effect
)
}


parameters <- c('eps1', 'beta1','beta2', 'beta3', 'beta4', 'eps2', 'survival')#'int', 'beta2',  

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.morpho <- jagsUI(jags.data, inits, parameters, "phi.morpho.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.morpho)
MCMCtrace(phi.morpho)
write.csv(phi.morpho$summary, './output/phi.antlersxsite.csv', row.names = T)


#create a tibble of the posterior draws
gather<- phi.morpho %>% gather_draws(survival[morpho, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$morpho == 2)[1] # 3601 values of antler 1

#unscale and uncenter weight
morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
morpho.sim.usc1 <- for (i in morpho.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$bcs <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=bcs, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Antlers (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.antlersxsite.jpg', phi.plot, width = 15, height = 10)


# ---- Model2: phi ~ age x site interaction ----


#survival as a function of weight (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

  beta1 ~ dnorm(0,0.1)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.01)  #birth year rain beta continuous
  
for ( u in 1:3) { #interaction between age and site
  beta4[u] ~ dnorm(0, 1)
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
}

for (u in 1:3){         #prior for random site effect
  eps1[u] ~ dlnorm(0,tau)
  }

for (u in 1:12){  #prior for year effect
  eps2[u] ~ dnorm(0,tau)
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
            logit(phi[i,t-1]) <- eps1[bs[i]] + beta1*ageclass[i,t-1]
                                      + beta2*morpho[i, t-1]
                                      + beta3*by.rain[i, t-1]
                                      + beta4[bs[i]]*ageclass[i,t-1]
                                      + eps2[year[i]]
                                      

        # Observation process
              ch[i,t] ~ dbern(mu2[i,t])
              mu2[i,t] <- p * z[i,t]

      
        
     } #t
   } #i
   
    

   # #derived parameters
   #    for (j in 1:100 ) { #age simulation, beta1
   #    for (l in 1:3){ #site, eps1
   # 
   #    survival[j,l] <- exp( beta1*age.sim[j]  + eps1[l] + beta4[l]*age.sim[j] )/ 
   #                          (1 + exp(beta1*age.sim[j]  + eps1[l] + beta4[l]*age.sim[j]  )) 
   #  #   } #for k
   #  # }   #for j
   #  } # for i
   #  } # for l


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(#int = rnorm(1,0,1), 
  z = z.init,
  eps1 = rlnorm(3,0,1), #birth site random effect
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = rnorm(3,0,1),#age and site interaction
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('eps1', 'beta1','beta2', 'beta3', 'beta4', 'eps2','p', 
                'sigma', "fit_obs_sum", "fit_rep_sum", "p")#'int', 'beta2',  

# MCMC settings
ni <- 5000
nt <- 10
nb <- 4000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age <- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
write.csv(phi.age$summary, './output/phi.agexsite.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.age %>% gather_draws(survival[age, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$age == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
age.sim.usc <- (age.sim * sd(data$ageclass, na.rm = T)) + mean(data$ageclass, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
age.sim.usc1 <- for (i in age.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$ageclass <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=ageclass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.1),          # x, y inside the plot area
        legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.agexsite.jpg', phi.plot, width = 15, height = 10)

