#survival analysis using on site separately 

#final run of survival analyses
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

# data<- read.csv('./cleaned/final.ch1.csv', header = T)
data <- read.csv('./cleaned/ch.carryoverrain.csv', header = T)
data <- data %>% filter(!birth_year < 2011 & !year < 2011) #filter for years where we have all three sites

data<- data %>% filter(bs == "ey") 
unique(data$bs)

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
# 
# #create birthsite vector
# id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
# bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy
# bs.2 <- bs
# bs.2[bs.2 %in% c(3)] <-2 #binary variable of dmp versus pasture
# sum(bs.2==1) #182 individuals in DMP
# sum(bs.2==2) #307 individuals in pasture

#create ageclass matrix with continuous age cov
data$age.sc <- scale(data$ageclass)
age.sc<- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc <- as.matrix(age.sc[,-1])


#create ageclass matrix treating age as categorical
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

# #create age class of juvenile, mature, geriatric
# age.bin <- ageclass
# age.bin[age.bin %in% c(2,3)] <- 1 #juveniles 1.5 - 3.5 year old
# age.bin[age.bin %in% c(4,5,6,7,8)] <- 2 #mature 4.5 - 8.5 year old
# age.bin[age.bin %in% c(9,10,11,12,13,14,15)] <- 3 #geriatric 9.5 and older

#create birth year vector
birthyear <- as.numeric(as.factor(data$birth_year))

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


#carryover rainfall
cy.prior.rain<-pivot_wider(data, names_from = 'year', values_from = 'cy.rain.one.prior.sc', id_cols = 'animal_id' )
cy.prior.rain<-as.matrix(cy.prior.rain[,-1])

by.prior.rain<-pivot_wider(data, names_from = 'year', values_from = 'by.rain.one.prior.sc', id_cols = 'animal_id' )
by.prior.rain<-as.matrix(by.prior.rain[,-1])

##create simulated birthyear prior rainfall vector
nvalues <- 100
by.rain.prior.sim <- seq(from = min(by.prior.rain, na.rm = T), to = max(by.prior.rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of carryover rainfall in data


##create simulated capture year prior rainfall vector
nvalues <- 100
cy.rain.prior.sim <- seq(from = min(cy.prior.rain, na.rm = T), to = max(cy.prior.rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of carryover rainfall in data

#lets make rainfall categorical with low intermediate and high rainfall based upon quantiles
data <- data %>%
  mutate(cy.rain.group = ntile(cy.rain, 3))

data <- data %>%
  mutate(by.rain.group = ntile(by.rain, 3))

by.rain.group <- data$by.rain.group
cy.rain.group <- data$cy.rain.group

by.rain.group <- unique(data[, c("animal_id", "by.rain.group")])
by.rain.group <- as.numeric(factor(by.rain.group$by.rain.group)) # 1 = dmp, 2 = ey, 3 = wy
# by.rain.group <- by.rain.group %>% #two duplicates in the data
#   distinct(animal_id, .keep_all = TRUE)
# by.rain.group <- by.rain.group$by.rain.group


# ---- Model1: phi ~ weight x site interaction, fixed effect ----

#survival as a function of weight (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

  beta1 ~ dnorm(0,0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous

  
for ( u in 1:3) { 
  beta4[u] ~ dlnorm(0, 0.1)#interaction between morpho and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*by.rain[i, t-1]
                                  + beta4[bs[i]]*morpho[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
            
         
            
      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #weight simulation, beta2
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( beta2*morpho.sim[i]  + beta5[j] + beta4[j]*morpho.sim[i] )/
                            (1 + exp( beta2*morpho.sim[i]  + beta5[j] + beta4[j]*morpho.sim[i]))
      
    } # for j
    } # for l
    
      for (i in c(1,50,100)){
      for (j in 1:3){
        
        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(#int = rnorm(1,0,1), 
  z = z.init,
  morpho = weight.init, #initial values for NA morphos
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = rlnorm(3,0,1),#morpho and site interaction
  beta5 = rnorm(3, 0, 1), #site effect
  eps2 = rnorm(12, 0, 1) #capture year random effect
)
}


parameters <- c('beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival'
)  

# MCMC settings
ni <- 15000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.weight <- jagsUI(jags.data, inits, parameters, "phi.weight.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.weight)
MCMCtrace(phi.weight)
write.csv(phi.weight$summary, './output/phi.weighttxsite.fe.csv', row.names = T)


#create a tibble of the posterior draws
gather<- phi.weight %>% gather_draws(survival[morpho, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$morpho == 2)[1] # 4500 values of antler 1

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
ggsave('./figures/phi.weightxsite.fe.jpg', phi.plot, width = 10, height = 10)


# ---- Model2: phi ~ antlers x site interaction , fixed effect----

#survival as a function of antlers (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)

sink("phi.antlers.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

  beta1 ~ dnorm(0,0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous

  
for ( u in 1:3) { 
  beta4[u] ~ dlnorm(0, 0.1)#interaction between morpho and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*by.rain[i, t-1]
                                  + beta4[bs[i]]*morpho[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
            
         
            
      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #antler simulation, beta2
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( beta2*morpho.sim[i]  + beta5[j] + beta4[j]*morpho.sim[i] )/
                            (1 + exp( beta2*morpho.sim[i]  + beta5[j] + beta4[j]*morpho.sim[i]))
      
    } # for j
    } # for l
    
      for (i in c(1,50,100)){
      for (j in 1:3){
        
        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = antler.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_antlers, occasions = occasions_antlers,morpho = antlers, 
                  year = capyear)

# Initial values
inits <- function(){list(#int = rnorm(1,0,1), 
  z = z.init,
  morpho = antlers.init, #initial values for NA morphos
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = rlnorm(3,0,1),#morpho and site interaction
  beta5 = rnorm(3,0,1), #site effect
  eps2 = rnorm(12, 0, 1)#capture year random effect
)
}


parameters <- c( 'beta1','beta2', 'beta3', 'beta4', 'eps2', 'surv_diff', 'survival')#'int', 'beta2',  

# MCMC settings
ni <- 15000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.antlers <- jagsUI(jags.data, inits, parameters, "phi.antlers.jags", n.chains = nc,
                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.antlers)
MCMCtrace(phi.antlers)
write.csv(phi.antlers$summary, './output/phi.antlersxsite.fe.csv', row.names = T)


#create a tibble of the posterior draws
gather<- phi.antlers %>% gather_draws(survival[morpho, site]) #this creates a dataframe in long format with indexing
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
ggsave('./figures/phi.antlersxsite.fe.jpg', phi.plot, width = 10, height = 10)


# ---- Model3: phi ~ age x site interaction, fixed effect ----


#survival as a function of weight (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0,0.001)
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous

  
for ( u in 2:3) { 
  beta4[u] ~ dnorm(0, 0.001)#interaction between age and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*by.rain[i, t-1]
                                  + beta4[bs[i]]*ageclass[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
            
         
            
      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #age simulation, beta1
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta1*age.sim[i]  + beta5[j] + beta4[j]*age.sim[i] )/
                            (1 + exp( int+ beta1*age.sim[i]  + beta5[j] + beta4[j]*age.sim[i]))

    } # for j
    } # for l

      for (i in c(1,50,100)){
      for (j in 1:3){

        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  
  int = rnorm(1,0,1), 
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#age and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1) #random effect for capture year
  
)
}

parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')  

# MCMC settings
ni <- 20000
nt <- 10
nb <- 12000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age <- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
write.csv(phi.age$summary, './output/phi.agexsite.fe.csv', row.names = T)

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
        legend.position.inside = c(0.9,0.9),          # x, y inside the plot area
        legend.justification = c("right", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.agexsite.fe.jpg', phi.plot, width = 10, height = 10)



# 
# ---- Model4: phi ~ birth year rain x site interaction, fixed effect ----

# Specify model in JAGS language
set.seed(100)
sink("phi.rain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous

  
for ( u in 2:3) { 
  beta4[u] ~ dnorm(0, 0.001)#interaction between rain and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*rain[i, t-1]
                                  + beta4[bs[i]]*rain[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #rain simulation, beta3
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i] )/
                            (1 + exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i]))

    } # for j
    } # for l

      for (i in c(1,50,100)){
      for (j in 1:3){

        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.rain <- jagsUI(jags.data, inits, parameters, "phi.rain.jags", n.chains = nc,
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.rain)
MCMCtrace(phi.rain)
write.csv(phi.rain$summary, './output/phi.rainxsite.fe.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.rain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (by.rain.sim * sd(data$by.rain, na.rm = T)) + mean(data$by.rain, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$byrainfall <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=byrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Birthyear Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.byrainxsite.fe.jpg', phi.plot, width = 10, height = 10)



# ---- Model5: phi ~ capture year rain, fixed effect ----

# Specify model in JAGS language
set.seed(100)
sink("phi.rain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)

  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #rain beta continuous


for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
}

# for (u in 1:10){  #prior for year effect
#   eps2[u] ~ dnorm(0,tau)
# }

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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*rain[i, t-1]
                                  #+ eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #rain simulation, beta3

      survival[i] <- exp( int + beta3*rain.sim[i] )/
                            (1 + exp( int + beta3*rain.sim[i]))

    } # for l

      # for (i in c(1,50,100)){
      # for (j in 1:3){
      # 
      #   surv_diff[i,j] <- survival[i,j] - survival[i,1]
      # }
      # }


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, rain = cy.rain,
                  age.sim = age.sim, morpho.sim = weight.sim, rain.sim = cy.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1) # rain beta
  # beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  # beta5 = c(NA, rnorm(2,0,1)), #site effect
 # eps2 = rnorm(10, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'survival')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.rain <- jagsUI(jags.data, inits, parameters, "phi.rain.jags", n.chains = nc,
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.rain)
MCMCtrace(phi.rain)
write.csv(phi.rain$summary, './output/phi.cyrain.ey.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.rain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (cy.rain.sim * sd(data$cy.rain, na.rm = T)) + mean(data$cy.rain, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$cyrainfall <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=cyrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Capture Year Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.cyrainxsite.fe.jpg', phi.plot, width = 10, height = 10)



# 
# ---- Model6: phi ~ carryover birth year rain x site interaction, fixed effect ----


# Specify model in JAGS language
set.seed(100)
sink("phi.prior.rain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #birth year prior rain beta continuous

  
for ( u in 2:3) { 
  beta4[u] ~ dnorm(0, 0.001)#interaction between rain and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*rain[i, t-1]
                                  + beta4[bs[i]]*rain[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #rain simulation, beta3
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i] )/
                            (1 + exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i]))

    } # for j
    } # for l

      for (i in c(1,50,100)){
      for (j in 1:3){

        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, 
                  rain = by.prior.rain,rain.sim = by.rain.prior.sim,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, 
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.priorbyrain <- jagsUI(jags.data, inits, parameters, "phi.prior.rain.jags", n.chains = nc,
                          n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.priorbyrain)
MCMCtrace(phi.priorbyrain)
write.csv(phi.priorbyrain$summary, './output/phi.byraincarryover.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.priorbyrain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (by.rain.prior.sim * sd(data$by.rain.one.prior, na.rm = T)) + mean(data$by.rain.one.prior, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$bypriorrainfall <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=bypriorrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Birthyear Prior Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.bypriorrainxsite.fe.jpg', phi.plot, width = 10, height = 10)



# ---- Model7: phi ~ carryover capture year rain x site interaction, fixed effect ----

# Specify model in JAGS language
set.seed(100)
sink("phi.prior.cyrain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  beta3 ~ dnorm(0,0.001)  #rain beta continuous

  
for ( u in 2:3) { 
  beta4[u] ~ dnorm(0, 0.001)#interaction between rain and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3*rain[i, t-1]
                                  + beta4[bs[i]]*rain[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:100 ) { #rain simulation, beta3
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i] )/
                            (1 + exp( int + beta3*rain.sim[i]  + beta5[j] + beta4[j]*rain.sim[i]))

    } # for j
    } # for l

      for (i in c(1,50,100)){
      for (j in 1:3){

        surv_diff[i,j] <- survival[i,j] - survival[i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, 
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, 
                  rain = cy.prior.rain, rain.sim = cy.rain.prior.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.priorcyrain <- jagsUI(jags.data, inits, parameters, "phi.prior.cyrain.jags", n.chains = nc,
                          n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.priorcyrain)
MCMCtrace(phi.priorcyrain)
write.csv(phi.priorcyrain$summary, './output/phi.priorcyrainxsite.fe.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.priorcyrain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (cy.rain.prior.sim * sd(data$cy.rain.one.prior, na.rm = T)) + mean(data$cy.rain.one.prior, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$cypriorrainfall <- vector

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=cypriorrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Capture Year Prior Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.priorcyrainxsite.fe.jpg', phi.plot, width = 10, height = 10)



# ---- Model8: phi ~ capture year rain(cat) x site interaction, fixed effect ----

# Specify model in JAGS language
set.seed(100)
sink("phi.rain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)
  beta3[1] <- 0
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous

  
for ( u in 2:3) { 
  beta3[u] ~ dnorm(0, 0.001) #rain beta categorical
  beta4[u] ~ dnorm(0, 0.001)#interaction between rain and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3[rain[i]]
                                  + beta4[bs[i]]*rain[i]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:3 ) { #rain  beta3
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta3[i]  + beta5[j] + beta4[j]*i )/
                            (1 + exp( int + beta3[i]  + beta5[j] + beta4[j]*i))

    } # for j
    } # for l

      # for (i in c(1,50,100)){
      # for (j in 1:3){
      # 
      #   surv_diff[i,j] <- survival[i,j] - survival[i,1]
      # }
      # }


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, rain = cy.rain.group,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = cy.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = c(NA, rnorm(2, 0, 1)), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.rain <- jagsUI(jags.data, inits, parameters, "phi.rain.jags", n.chains = nc,
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.rain)
MCMCtrace(phi.rain)
write.csv(phi.rain$summary, './output/phi.cyrainxsite.group.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.rain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)
gather$rain <- as.factor(gather$rain)

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=rain, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = position_dodge(width = .4))+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Capture Year Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.cyrainxsite.cat.jpg', phi.plot, width = 10, height = 10)


# ---- Model9: phi ~ birth year rain(cat) x site interaction, fixed effect ----

# Specify model in JAGS language
set.seed(100)
sink("phi.rain.jags")
cat("
model {

#prior for recapture prob
  p ~ dbeta(1, 1)
  
#priors
  int ~ dnorm(0,0.001)
  beta3[1] <- 0
  beta4[1] <- 0
  beta5[1] <- 0
  
  beta1 ~ dnorm(0, 0.001)  # age beta continuous
  beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous

  
for ( u in 2:3) { 
  beta3[u] ~ dnorm(0, 0.001) #rain beta categorical
  beta4[u] ~ dnorm(0, 0.1)#interaction between rain and site
  beta5[u] ~ dnorm(0, 0.001) #site effect
}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
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
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]
                                  + beta2*morpho[i, t-1]
                                  + beta3[rain[i]]
                                  + beta4[bs[i]]*rain[i]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (i in 1:3 ) { #rain  beta3
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta3[i]  + beta5[j] + beta4[j]*i )/
                            (1 + exp( int + beta3[i]  + beta5[j] + beta4[j]*i))

    } # for j
    } # for l

      # for (i in c(1,50,100)){
      # for (j in 1:3){
      # 
      #   surv_diff[i,j] <- survival[i,j] - survival[i,1]
      # }
      # }


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = age.sc, rain = by.rain.group,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = cy.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), #age beta
  beta2 = rlnorm(1, 0, 1),#morpho beta
  beta3 = c(NA, rnorm(2, 0, 1)), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#rain and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = rnorm(12, 0, 1)
)
}

parameters <- c('int','beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.rain <- jagsUI(jags.data, inits, parameters, "phi.rain.jags", n.chains = nc,
                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.rain)
MCMCtrace(phi.rain)
write.csv(phi.rain$summary, './output/phi.byrainxsite.group.csv', row.names = T)

#prepare to create ggplot with posteriors
#create a tibble of the posterior draws
gather<- phi.rain %>% gather_draws(survival[rain, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)
gather$rain <- as.factor(gather$rain)

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=rain, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = position_dodge(width = .4))+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Birth Year Rainfall", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
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
ggsave('./figures/phi.byrainxsite.cat.jpg', phi.plot, width = 10, height = 10)



# ---- Model10: phi ~ age (categorical) x site interaction, fixed effect ----


#survival as a function of weight (continuous), age (continuous), site random effect, age x site, random effect of capture year
# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0,0.001)
  beta1[1] <- 0
  beta4[1] <- 0
  beta5[1] <- 0
  eps2[1] <- 0

for (u in 2:11){
  beta1[u] ~ dnorm(0, 0.001) #ageclass beta

}

  # beta2 ~ dlnorm(0, 0.1)    # morpho beta continuous
  # beta3 ~ dnorm(0,0.001)  #birth year rain beta continuous

  
for ( u in 2:3) { 
  beta5[u] ~ dnorm(0, 0.001) #site effect
    beta4[u] ~ dnorm(0, 0.001)#interaction between age and site

}

for (u in 1:nind){      #prior for missing morphometrics
  for (j in 1:occasions[u]){
  morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
     }
}

for (u in 2:12){  #prior for year effect
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
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]
                                  # + beta2*morpho[i, t-1]
                                  # + beta3*by.rain[i, t-1]
                                  + beta4[bs[i]]*ageclass[i,t-1]
                                  + beta5[bs[i]]
                                  + eps2[year[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
            
         
            
      } #t
   } #i

   #derived parameters
      for (i in 1:11 ) { #age beta1
      for (j in 1:3){ #site, beta5

      survival[i,j] <- exp( int + beta1[i]  + beta5[j] + beta4[j]*i )/
                            (1 + exp( int+ beta1[i]  + beta5[j] + beta4[j]*i))

    } # for j
    } # for l

      # for (i in 1:11){ # age
      # for (j in 1:3){ #site
      # 
      #   surv_diff[i,j] <- survival[i,j] - survival[i,1]
      # }
      # }


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ncol = ncol(ch), ageclass = ageclass, by.rain = by.rain,
                  bs = bs, age.sim = age.sim, morpho.sim = weight.sim, rain.sim = by.rain.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,morpho = weight, year = capyear)

# Initial values
inits <- function(){list(
  
  int = rnorm(1,0,1), 
  z = z.init,
  morpho = weight.init,
  beta1 = c(NA, rnorm(10,0,1)), #age beta
  # beta2 = rlnorm(1, 0, 1),#morpho beta
  # beta3 = rnorm(1, 0, 1), # rain beta
  beta4 = c(NA, rnorm(2,0,1)),#age and site interaction
  beta5 = c(NA, rnorm(2,0,1)), #site effect
  eps2 = c(NA, rnorm(11, 0, 1)) #random effect for capture year
  
)
}

parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps2', 'surv_diff', 'survival')  

# MCMC settings
ni <- 25000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age <- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
write.csv(phi.age$summary, './output/phi.agexsite.fe.csv', row.names = T)

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
        legend.position.inside = c(0.9,0.9),          # x, y inside the plot area
        legend.justification = c("right", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.agexsite.fe.jpg', phi.plot, width = 10, height = 10)

