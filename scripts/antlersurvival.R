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
h <- replace(h, is.infinite(h), 15)
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
NA_indices <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
for(i in 1:nrow(indices)){
  NA_indices[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
}
weight<-as.matrix(weight)

#how many occasions does each individual have of an NA weight
occasions <- rowSums(is.na(weight))

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

# 
# 
# #survival as a function of weight (continuous), age (1-15), site (all 3), age x site, random effect of capture year
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
# 
# int ~ dnorm(0,0.01)
# age.beta[1] <- 0      #ageclass
# bs.beta[1] <- 0       #birth site
# age.bs.beta[1] <- 0   #age site interaction
# eps.capyear[1] <- 0   #capture year random effect
# 
# for (u in 2:15){
#   age.beta[u] ~ dnorm(0,0.01)
# }
# 
# for (u in 2:3){
#   bs.beta[u] ~ dnorm(0, 0.01)
# }
# 
# for (u in 2:15){
#   age.bs.beta[u] ~ dnorm(0, 0.01)
# }
# 
# weight.beta ~ dunif(0, 100)
# antlers.beta ~ dunif(0,100)
# by.rain.beta ~ dnorm(0,0.001)
# 
# for (u in 1:nind){
#   for (j in 1:occasions[u]){  #prior for missing weights
#   weight[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
#      }
# }
# 
# 
# for (u in 1:nind){
#   for (j in 1:occasions_antlers[u]){  #prior for missing weights
#   antlers[u,NA_indices_antlers[u,j]] ~ dnorm( 0, 0.01)
#      }
# }
# 
# tau <- 1/(sigma*sigma)
# sigma ~ dunif(0,100)
# 
# for (u in 2:15){
#   eps.capyear[u] ~ dnorm(0, tau.capyear)
# }
#   tau.capyear <- 1/(sigma.capyear*sigma.capyear)
#   sigma.capyear ~ dunif(0,100)
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
#             logit(phi[i,t-1]) <- int  + age.beta[ageclass[i,t-1]]
#                                       + bs.beta[bs[i]] + age.bs.beta[ageclass[i,t-1]]*bs[i]
#                                       + eps.capyear[capyear[i]]
#                                       + weight.beta*weight[i, t-1]
#                                       + antlers.beta*antlers[i, t-1]
#                                       + by.rain.beta*by.rain[i, t-1]
#                                       
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
#       } #t
#    } #i
# 
# #derived parameters
#     for (j in 1:10){ #ageclass
#       for (k in 1:3){ #site
#       survival[j,k] <- exp(int + age.beta[j] + bs.beta[k] + age.bs.beta[j]*k)/
#                             (1 + exp(int + age.beta[j] + bs.beta[k] + age.bs.beta[j]*k))
#       } #for k
#     }   #for j
# 
# 
#     for (j in 1:9){
#       for (k in 1:3){
#         surv_diff[j,k] <- survival[j+1, k] - survival[j,k]
#       }
#     }
# 
#     for (j in 1:10){
#       for (k in 1:3){
#         site_diff[j,k] <- survival[j, 1] - survival[j,2]
#       }
#     }
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
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, bs = bs, by.rain = by.rain,
#                   capyear = capyear, weight = weight, NA_indices = NA_indices, occasions = occasions,
#                   antlers = antlers, NA_indices_antlers = NA_indices_antlers, occasions_antlers = occasions_antlers)
# 
# # Initial values
# inits <- function(){list(int = rnorm(1,0,1), z = z.init,
#                          age.beta = c(NA, rnorm(14,0,1)),
#                          bs.beta = c(NA, rnorm(2,0,1)),
#                          age.bs.beta = c(NA, rnorm(14,0,1)),
#                          eps.capyear = c(NA, runif(14, 0, 100)),
#                          weight.beta = runif(1, 0, 100),
#                          weight = weight.init,
#                          antlers.beta = runif(1, 0, 100),
#                          antlers = antlers.init,
#                          by.rain.beta = rnorm(1,0,1))}
# 
# 
# parameters <- c('int', 'age.beta', 'bs.beta', 'age.bs.beta', 'by.rain.beta', 'weight.beta', 'antlers.beta',
#                 'survival','surv_diff', 'site_diff')#, 
# 
# # MCMC settings
# ni <- 40000
# nt <- 10
# nb <- 30000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.age <- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
#                   n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.age)
# MCMCtrace(phi.age)
# write.csv(phi.age$summary, './output/phi.age.site.weight.ant.byrain.csv', row.names = T)
# 


####

#survival as a function of antlers (continuous), age (1-15), site (all 3)
# Specify model in JAGS language
set.seed(100)
sink("phi.antlers.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

int ~ dnorm(0,0.01)
age.beta[1] <- 0      #ageclass
bs.beta[1] <- 0       #birth site
antlers.bs.beta[1] <- 0   #antler site interaction
# eps.capyear[1] <- 0   #capture year random effect

for (u in 2:15){
  age.beta[u] ~ dnorm(0,0.01)
}
# 
for (u in 2:3){
  bs.beta[u] ~ dnorm(0, 0.01)
}

for (u in 2:3){
  antlers.bs.beta[u] ~ dnorm(0, 0.01) #antler birthsite interaction
}

antlers.beta ~ dunif(0, 100)
# by.rain.beta ~ dnorm(0,0.001)

for (u in 1:nind){
  for (j in 1:occasions_antlers[u]){  #prior for missing weights
  antlers[u,NA_indices_antlers[u,j]] ~ dnorm( 0, 0.01)
     }
}


tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

# 
# for (u in 2:15){
#   eps.capyear[u] ~ dnorm(0, tau.capyear)
# }
#   tau.capyear <- 1/(sigma.capyear*sigma.capyear)
#   sigma.capyear ~ dunif(0,100)

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int  + age.beta[ageclass[i,t-1]]
                                      + bs.beta[bs[i]] 
                                      #+ age.bs.beta[ageclass[i,t-1]]*bs[i]
                                      # + eps.capyear[capyear[i]]
                                      # + weight.beta*weight[i, t-1]
                                      + antlers.beta*antlers[i, t-1]
                                      + antlers.bs.beta*antlers[i,t-1]*bs[i]
                                      # + by.rain.beta*by.rain[i, t-1]
                                      
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i
# 
#derived parameters
    # for (j in 1:10){ #ageclass
      for (k in 1:3){ #site
        for (i in 1:100 ) { #antler simulation
      survival[i,k] <- exp(int + bs.beta[k] + antlers.beta*antler.sim[i] )/
                            (1 + exp(int  + bs.beta[k] + antlers.beta*antler.sim[i] ))
      } #for k 
   # }   #for j+ age.beta[j]
    } # for i 
# + age.bs.beta[j]*k+ age.bs.beta[j]*k
# 
#     for (j in 1:9){
#       for (k in 1:3){
#         surv_diff[j,k] <- survival[j+1, k] - survival[j,k]
#       }
#     }
# 
#     for (j in 1:10){
#       for (k in 1:3){
#         site_diff[j,k] <- survival[j, 1] - survival[j,2]
#       }
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, bs = bs,# by.rain = by.rain,
                  #capyear = capyear, weight = weight, NA_indices = NA_indices, occasions = occasions,
                  antlers = antlers, NA_indices_antlers = NA_indices_antlers, occasions_antlers = occasions_antlers,
                  antler.sim = antler.sim)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init,
                         age.beta = c(NA, rnorm(14,0,1)),
                         bs.beta = c(NA, rnorm(2,0,1)),
                         # age.bs.beta = c(NA, rnorm(14,0,1)),
                         # eps.capyear = c(NA, runif(14, 0, 100)),
                         # weight.beta = runif(1, 0, 100),
                         # weight = weight.init,
                         antlers.beta = runif(1, 0, 100),
                         antlers = antlers.init)}
#  by.rain.beta = rnorm(1,0,1))}


parameters <- c('int', 'age.beta','antlers.beta','bs.beta', 'survival') #, 'age.bs.beta', 'by.rain.beta', 'weight.beta', 'antlers.beta',
# 'survival','surv_diff', 'site_diff')#, 

# MCMC settings
ni <- 10000
nt <- 10
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.antlers <- jagsUI(jags.data, inits, parameters, "phi.antlers.jags", n.chains = nc,
                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.antlers)
MCMCtrace(phi.antlers)
write.csv(phi.age$summary, './output/phi.age.site.weight.ant.byrain.csv', row.names = T)


#create a tibble of the posterior draws
gather<- phi.antlers %>% gather_draws(survival[antler, site]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$antler == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
antler.sim.usc <- (antler.sim * sd(data$bcsin, na.rm = T)) + mean(data$bcsin, na.rm = T)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector <- numeric(0)
antler.sim.usc1 <- for (i in antler.sim.usc) {
  rep_i <- rep(i, times = 4500)
  vector <- c(vector,rep_i)
  
}

gather$antlerscore <- vector

#plot for ageclass 7

phi.plot.7<- gather %>% 
  # subset(gather, age %in% '7') %>%
  ggplot(aes(x=antlerscore, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Antler Score (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Survival by antler score and site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position= c(.8,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
ggsave('./figures/phi.site.rain.jpg', phi.plot.7, width = 15, height = 10)
