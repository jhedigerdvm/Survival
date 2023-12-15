#Survival model with interaction between rainfall and birth site 
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(here)

data<- read.csv('./cleaned/caphx.rainfall.nov.oct1.csv', header = T)

ch<- pivot_wider(data, names_from = 'year', values_from = 'status', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates 
ch[indices] <- 1


# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0))
f <- apply(ch, 1, get.first) 

#create birthsite vector
id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass matrix
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

#create animal id vector
id <- as.numeric(factor(id.bs.by$animal_id))

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


annual.rainfall<-pivot_wider(data, names_from = 'year', values_from = 'annual.sc', id_cols = 'animal_id' )
annual.rainfall<-annual.rainfall[,-1]
annual.rainfall<-as.matrix(annual.rainfall)

spring.rainfall<- pivot_wider(data, names_from = 'year', values_from = 'sum.march.apr.may.sc', id_cols = 'animal_id' )
spring.rainfall<-spring.rainfall[,-1]
spring.rainfall<-as.matrix(spring.rainfall)


summer.rainfall<- pivot_wider(data, names_from = 'year', values_from = 'sum.jun.jul.aug.sc', id_cols = 'animal_id' )
summer.rainfall<-summer.rainfall[,-1]
summer.rainfall<-as.matrix(summer.rainfall)

nvalues <- 1000
summer.rain.sim <- seq(from = min(summer.rainfall, na.rm = T), to = max(summer.rainfall, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
spring.rain.sim <- seq(from = min(spring.rainfall, na.rm = T), to = max(spring.rainfall, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
rain.sim <- seq(from = min(annual.rainfall, na.rm = T), to = max(annual.rainfall, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
rain.sim

#model evaluating the effects of total annual rainfall, birthsite and the interaction between site and rainfall
# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
int ~ dnorm(0,0.01)

site.beta[1] <- 0

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      survival[j,i] <- exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta', 'rain.site.beta', 'p','survival' )

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site <- jagsUI(jags.data, inits, parameters, "cjs-rain-site.jags", n.chains = nc,
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site)
# 
# ########################################################3
# #effects of spring rainfall
# # Specify model in JAGS language
# set.seed(100)
# sink("cjs-springrain-site.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
# int ~ dnorm(0,0.01)
# 
# site.beta[1] <- 0
# 
# for (u in 2:3){
#   site.beta[u] ~ dnorm(0,0.01)
# }
# 
# rain.beta ~ dnorm(0,0.01)
# 
# 
# rain.site.beta[1] <-0
# 
# for (u in 2:3){
#   rain.site.beta[u] ~ dnorm(0,0.01)
# }
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
#             logit(phi[i,t-1]) <- int + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1]
# 
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
#       } #t
#    } #i
# 
# #derived parameters
#   for (i in 1:3){ #site and rain.site.beta
#     for (j in 1:1000){ #rain
#       survival[j,i] <- exp(int + site.beta[i] + rain.beta*spring.rain.sim[j] + rain.site.beta[i]*spring.rain.sim[j])/
#                             (1 + exp(int + site.beta[i] + rain.beta*spring.rain.sim[j] + rain.site.beta[i]*spring.rain.sim[j]))
#       }                     #delta method to convert from logit back to probability Powell et al. 2007
#     }
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
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = spring.rainfall, spring.rain.sim = spring.rain.sim)
# 
# # Initial values
# inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1),
#                          site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)))} #
# 
# parameters <- c('int', 'site.beta', 'rain.beta', 'rain.site.beta', 'p','survival' )
# 
# # MCMC settings
# ni <- 1000
# nt <- 1
# nb <- 500
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# cjs.springrain.site <- jagsUI(jags.data, inits, parameters, "cjs-springrain-site.jags", n.chains = nc,
#                         n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
# 
# print(cjs.springrain.site)
# 
# 
# ########################################################3
# #effects of summer rainfall
# # Specify model in JAGS language
# set.seed(100)
# sink("cjs-summer-rain-site.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
# int ~ dnorm(0,0.01)
# 
# site.beta[1] <- 0
# 
# for (u in 2:3){
#   site.beta[u] ~ dnorm(0,0.01)
# }
# 
# rain.beta ~ dnorm(0,0.01)
# 
# 
# rain.site.beta[1] <-0
# 
# for (u in 2:3){
#   rain.site.beta[u] ~ dnorm(0,0.01)
# }
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
#             logit(phi[i,t-1]) <- int + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1]
# 
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
#       } #t
#    } #i
# 
# #derived parameters
#   for (i in 1:3){ #site and rain.site.beta
#     for (j in 1:1000){ #rain
#       survival[j,i] <- exp(int + site.beta[i] + rain.beta*summer.rain.sim[j] + rain.site.beta[i]*summer.rain.sim[j])/
#                             (1 + exp(int + site.beta[i] + rain.beta*summer.rain.sim[j] + rain.site.beta[i]*summer.rain.sim[j]))
#       }                     #delta method to convert from logit back to probability Powell et al. 2007
#     }
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
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = summer.rainfall, summer.rain.sim = summer.rain.sim)
# 
# # Initial values
# inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1),
#                          site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)))} #
# 
# parameters <- c('int', 'site.beta', 'rain.beta', 'rain.site.beta', 'p','survival' )
# 
# # MCMC settings
# ni <- 1000
# nt <- 1
# nb <- 500
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# cjs.summer.rain.site <- jagsUI(jags.data, inits, parameters, "cjs-summer-rain-site.jags", n.chains = nc,
#                               n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
# 
# print(cjs.summer.rain.site)


#model looking at rainfall, site, and AGE, and the interaction between site and rain, no random effects
# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site-age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
int ~ dnorm(0,0.01)

age.beta[1] <- 0
site.beta[1] <- 0

for (u in 2:14) {
  age.beta[u] ~ dnorm(0, 0.01)
}

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1]
          
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      for (k in 1:11){
        survival[j,i,k] <- exp(int + age.beta[k] + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + age.beta[k] + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
    }
  }
  
   for (i in 1:3){ #site
    survival_diff_rain1[i] <- survival[1,1,1] - survival[1,i,1]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain500[i] <- survival[500,1,1] - survival[500,i,1]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain1000[i] <- survival[1000,1,1] - survival[1000,i,1]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim, ageclass=ageclass)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1), age.beta = c(NA, rnorm(13, 0,1)),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta', 'age.beta', 'rain.site.beta', 'p','survival_diff_rain1','survival_diff_rain500','survival_diff_rain1000')#survival, 'survival

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site.age <- jagsUI(jags.data, inits, parameters, "cjs-rain-site-age.jags", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site.age)
write.csv(cjs.rain.site.age$summary, './output/rain.site.age.csv')
# 
# #create a tibble of the posterior draws
# gather<- cjs.rain.site.age %>% gather_draws(survival[rain,site,age]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 49500 values of rain 1 

#unscale and uncenter rain.sim
# rain.sim1 <- (rain.sim * sd(data$annual)) + mean(data$annual)

#create vector containing simulated rainfall data but in the format to sync up with gather 
vector <- numeric(0)
rain.sim2 <- for (i in rain.sim1) {
                rep_i <- rep(i, times = 99000)
                vector <- c(vector,rep_i)
  
}

gather$rain1 <- vector

#plot for all ages and facet wrap
phi.plot<- (gather %>%  
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) + #color equals line, fill equals ribbon
  facet_wrap(vars(age))+
  stat_lineribbon(.width = 0.95, alpha = 1/4)) #.width is the CRI, alpha is opacity

phi.plot.1<- 
  subset(gather, age %in% '1') %>%  
    ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
    stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI 
    scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
    scale_color_viridis_d()+ #color of line but no opacification
    labs(x = "RAINFALL (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Phi ~ int + age(1-2) + site + rain + rain*site")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.position = c(0.5,0.3),
          legend.title = element_blank(),
          legend.text = element_text(size = 28),
          plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
          axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
          axis.text = element_text(face='bold',size = 28),
          # axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
#plot 
phi.plot.ran.age<- gather1 %>% 
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI 
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Phi ~ int + site + rain + rain*site + (1|ageclass)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

phi.plot.6<- subset(gather, age %in% '6') %>%  
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age))+
  stat_lineribbon(.width = 0.95)+
  scale_fill_viridis_d(alpha = .2) +
  scale_color_viridis_d()

phi.plot.10<- subset(gather, age %in% '10') %>%  
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age))+
  stat_lineribbon(.width = 0.95)+
  scale_fill_viridis_d(alpha = .2) +
  scale_color_viridis_d()

ggsave('./figures/phi_rain_site_age1.jpg', phi.plot.1, width = 15, height = 10)
ggsave('./figures/phi_rain_site_age6.jpg', phi.plot.6, width = 10, height = 8)

#other questions would be how does rainfall during birth year influence survival
#how does birthdate affect survival

#include random effects 

#model looking at rainfall, site, and age, random effect of capture year 
# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site-age-ran.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
int ~ dnorm(0,0.01)

age.beta[1] <- 0
site.beta[1] <- 0
eps.capyear[1] <- 0

for (u in 2:14) {
  age.beta[u] ~ dnorm(0, 0.01)
}

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

for (u in 2:14){
  eps.capyear[u] ~ dnorm(sigma, tau.capyear)
}

  sigma ~ dunif(0,10) #standard dev
  tau.capyear <- 1/(sigma*sigma) #precision = 1/variance

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1] + eps.capyear[capyear[t-1]]
          
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      for (k in 1:11){
        survival[j,i,k] <- exp(int + age.beta[k] + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + age.beta[k] + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim, ageclass=ageclass,
                  capyear = capyear)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1), age.beta = c(NA, rnorm(13, 0,1)),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)), eps.capyear = c(NA, rnorm(13,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta', 'age.beta', 'rain.site.beta', 'p','survival')

# MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site.age.ran <- jagsUI(jags.data, inits, parameters, "cjs-rain-site-age-ran.jags", n.chains = nc, 
                            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site.age.ran)


#model looking at rainfall, site, and interaction between rain and site, random effect of age 
# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site-ran-age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
int ~ dnorm(0,0.01)

site.beta[1] <- 0
eps.ageclass[1] <- 0

for (u in 2:14) {
  eps.ageclass[u] ~ dnorm(sigma, tau.ageclass)
}
  sigma ~ dunif(0,10) #standard dev
  tau.ageclass <- 1/(sigma*sigma) #precision = 1/variance

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int  + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1] + eps.ageclass[ageclass[i,t-1]]
          
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      
        survival[j,i] <- exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
  }
  
  for (i in 1:3){ #site
    survival_diff_rain1[i] <- survival[1,1] - survival[1,i]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain500[i] <- survival[500,1] - survival[500,i]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain1000[i] <- survival[1000,1] - survival[1000,i]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim, ageclass=ageclass,
                  capyear = capyear)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)), eps.ageclass = c(NA, rnorm(13,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta','rain.site.beta', 'p', 'survival_diff_rain1','survival_diff_rain500','survival_diff_rain1000','survival') #
                

# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site.ran.age <- jagsUI(jags.data, inits, parameters, "cjs-rain-site-ran-age.jags", n.chains = nc, 
                                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site.ran.age)
write.csv(cjs.rain.site.ran.age$summary, './output/rainsite.randomage.csv')

#create a tibble of the posterior draws
gather1<- cjs.rain.site.ran.age %>% gather_draws(survival[rain,site]) #this creates a dataframe in long format with indexing
gather1$site <- as.factor(gather1$site)

#find first row for 2nd rain value
first_idx <- which(gather1$rain == 2)[1] # 4500 values of rain 1 

#unscale and uncenter rain.sim
# rain.sim1 <- (rain.sim * sd(data$annual)) + mean(data$annual)

#create vector containing simulated rainfall data but in the format to sync up with gather 
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 9000)
  vector1 <- c(vector1,rep_i)
  
}

gather1$rain1 <- vector1


#plot 
phi.plot.ran.age<- gather1 %>% 
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI 
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Phi ~ int + site + rain + rain*site + (1|ageclass)")+
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.position = c(0.5,0.3),
          legend.title = element_blank(),
          legend.text = element_text(size = 28),
          plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
          axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
          axis.text = element_text(face='bold',size = 28),
          # axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/phi_rain_site_ran_age.jpg', phi.plot.ran.age, width = 15, height = 10)


#rainfall, site, interaction, random id and age 
# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site-ran-age-id.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
int ~ dnorm(0,0.01)

site.beta[1] <- 0
eps.ageclass[1] <- 0
eps.id[1] <- 0

for (u in 2:14) {
  eps.ageclass[u] ~ dnorm(sigma, tau.ageclass)
}
  sigma ~ dunif(0,10) #standard dev
  tau.ageclass <- 1/(sigma*sigma) #precision = 1/variance

for (u in 2:489) {
  eps.id[u] ~ dnorm(sigma, tau.id)
}
  tau.id <- 1/(sigma*sigma) #precision = 1/variance

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int  + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1] + 
                                      eps.ageclass[ageclass[i,t-1]] + eps.id[id[i]]
          
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      
        survival[j,i] <- exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
  }
  
  for (i in 1:3){ #site
    survival_diff_rain1[i] <- survival[1,1] - survival[1,i]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain500[i] <- survival[500,1] - survival[500,i]
  }
  
  for (i in 1:3){ #site
    survival_diff_rain1000[i] <- survival[1000,1] - survival[1000,i]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim, 
                  ageclass=ageclass, id = id)#capyear = capyear

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1), sigma = runif(1,0,1),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)), 
                         eps.ageclass = c(NA, rnorm(13,0,1)), eps.id = c(NA, rnorm(488,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta','rain.site.beta', 'p', 'survival_diff_rain1','survival_diff_rain500',
                'survival_diff_rain1000','survival') #


# MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site.ran.age.id <- jagsUI(jags.data, inits, parameters, "cjs-rain-site-ran-age-id.jags", n.chains = nc, 
                                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site.ran.age.id)
write.csv(cjs.rain.site.ran.age.id$summary, './output/rain.ran.age.id.csv')

#create a tibble of the posterior draws
gather2<- cjs.rain.site.ran.age.id %>% gather_draws(survival[rain,site]) #this creates a dataframe in long format with indexing
gather2$site <- as.factor(gather2$site)

#find first row for 2nd rain value
first_idx <- which(gather2$rain == 2)[1] # 9000 values of rain 1 

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual)) + mean(data$annual)

#create vector containing simulated rainfall data but in the format to sync up with gather 
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 9000)
  vector1 <- c(vector1,rep_i)
  
}

gather2$rain1 <- vector1



#plot 
phi.plot.ran.age<- gather2 %>% 
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI 
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Phi ~ int + site + rain + rain*site + (1|ageclass) + (1|id)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/phi_rain_site_ran_age_id.jpg', phi.plot.ran.age, width = 15, height = 10)

