#dummy weight survival model
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

data<- read.csv('cleaned/final.ch1.csv', header = T)
data$weight.sc <- scale(data$weight)
data$bcs.sc <- scale(data$bcsin)

#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'year', values_from = 'status', id_cols = 'animal_id' )
ch<-ch[,-1]
known.fate <- ch # 34 individuals
ch[ch == 2] <- 1 #convert any known fates to 1s 

ch<-as.matrix(ch)

#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight.sc', id_cols = 'animal_id' )
weight<- as.matrix(weight[,-1])# 

antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcs.sc', id_cols = 'animal_id' )
antlers<- as.matrix(antlers[,-1])

# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

# # create vector with last occasion for each individual, marked by 2, 15 for end of study
# # rework h to only include capture myopathy or harvest, do not censor natural mortality
# get.last<- function(x) min(which(x>1))
# h <- apply(known.fate,1,get.last)
# h <- replace(h, is.infinite(h), 15)
# h
# f-h #check for zero
# 

#Function for latent state   #is it still necessary to account for known fates of the 34 individuals in h?
known.state.cjs <- function(ch){
  state <- ch
  for ( i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[ i , ] == 1))
    n2 <- max(which(ch[ i , ] == 1))
    state[ i , n1:n2 ] <- 1
    state[ i , n1 ] <- NA
  }
  state[ state == 0] <- NA
  return(state)
}

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

#create birthsite vector
id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass matrix
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

#create birth year vector
birthyear <- as.numeric(as.factor(id.bs.by$birth_year))

#create capture year vector
capyear <- f

#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

# 
# 
# #####################################
# #model looking at phi as a function of weight 
# 
# set.seed(100)
# sink("cjs-weight.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta( 1 , 1 )
# 
# #priors
# int ~ dnorm( 0 , 0.0001 )
# 
# bs.beta[1] <- 0
# ageclass.beta[1] <- 0
# # eps.capyear[1] <- 0
# 
# for (u in 1:nind){
#   for (j in 1:occasions[u]){  #prior for missing weights
#   weight[u,NA_indices[u,j]] ~ dnorm( 0 , 0.0001 )
#      }
# }
# 
# weight.beta ~ dnorm( 0 , 0.0001 )
# 
# for (u in 2:15){                              #prior for ageclass
#     ageclass.beta[u] ~ dnorm( 0 , 0.0001 )
#   }
# 
# for (u in 2:3){                               #prior for birth site
#   bs.beta[u] ~ dnorm( 0 , 0.0001 )
# }
# 
# # for (u in 2:15){
# #   eps.capyear[u] ~ dnorm(0, tau)
# # }
# 
# tau <- 1/(sigma*sigma)
# sigma ~ dunif(0,100)
# 
# 
# # Likelihood
# for (i in 1:nind){
#    # Define latent state at first capture, we know for sure the animal is alive
#       z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):nocc){
#         # State process
#             z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
#             logit(phi[i,t-1]) <- int + weight.beta*weight[i,t-1] + ageclass.beta[ageclass[i,t-1]] + bs.beta[bs[i]] #+ eps.capyear[capyear[i]] 
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
#                                             
# 
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
#       } #t
#    } #i
# 
# #Derived parameters
# for (i in 1:100){ #weight.sim
#   for (j in 1:12) { #ageclass
#     for (k in 1:3){ # birthsites
#       survival[i,j,k] <- exp(int + weight.beta*weight.sim[i] + ageclass.beta[j] + bs.beta[k]) / 
#                                 (1+exp(int + weight.beta*weight.sim[i] + ageclass.beta[j] + bs.beta[k]))
#     }
# }
# }
#       
#   
# }
# ",fill = TRUE)
# sink()
# 
# 
# 
# # Bundle data
# jags.data <- list(ch = ch, f = f, nind = nrow(ch), nocc = ncol(ch), weight = weight, bs = bs, #capyear = capyear,
#                   occasions=occasions, NA_indices=NA_indices, ageclass = ageclass, weight.sim = weight.sim)#
# 
# # Initial values
# inits <- function(){list(weight = weight.init, weight.beta = rnorm(1,0,1), z=known.state.cjs(ch), #eps.capyear = c(NA, rnorm(14,0,1)),
#                          ageclass.beta =  c(NA, rnorm(14,0,1)), bs.beta = c(NA, rnorm(2,0,1)), int = rnorm(1,0,1) )} #, 
# 
# parameters <- c('int', 'bs.beta', 'weight.beta', 'ageclass.beta', 'survival')#, 'eps.capyear'
# 
# # MCMC settings
# ni <- 40000
# nt <- 10
# nb <- 30000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
#                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
# # traceplot(cjs.weight)
# print(cjs.weight)
# MCMCtrace(cjs.weight)
# 
# 
# write.csv(cjs.weight$summary, 'weight.csv', row.names = T)
# 
# #########################################################################################
# 
# #create a tibble of the posterior draws
# gather<- cjs.weight %>% gather_draws(survival[weight,age, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$weight == 2)[1] # 108000 values of rain 1
# 
# #unscale and uncenter rain.sim
# weight.sim1 <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
# 
# #create vector containing simulated rainfall data but in the format to sync up with gather
# vector <- numeric(0)
# weight.sim2 <- for (i in weight.sim1) {
#   rep_i <- rep(i, times = 108000)
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$weight1 <- vector
# 
# #plot for ageclass 7
# 
# phi<-
#   subset(gather, age %in% '7') %>%
#   ggplot(aes(x=weight1, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
#   labs(x = "WEIGHT (POUNDS)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Survival by weight")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.15,0.1),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# ggsave('./figures/weight.jpg', phi, width = 15, height = 10)
# 
# 


#####################################
#model looking at phi as a function of weight and the interaction of weight with age

set.seed(100)
sink("cjs-weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta( 1 , 1 )

#priors

int ~ dnorm( 0 , 0.1 )
bs.beta[1] <- 0
bs.weight.beta[1] <- 0

for (u in 2:3) { #ageclass and weight interaction
  bs.weight.beta[u] ~ dnorm(0, 0.01)
}

for (u in 1:nind){
  for (j in 1:occasions[u]){  #prior for missing weights
  weight[u,NA_indices[u,j]] ~ dnorm( 0 , 0.0001 )
     }
}

weight.beta ~ dnorm( 0 , 1 )

for (u in 2:3){                               #prior for birth site
  bs.beta[u] ~ dnorm( 0 , 0.0001 )
}


tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):nocc){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            logit(phi[i,t-1]) <- int + weight.beta*weight[i,t-1] + bs.beta[bs[i]] 
                                        + bs.weight.beta[bs[i]]*weight[i,t-1]
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
                                            

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i
# 
# #Derived parameters
# for (i in 1:100){ #weight.sim
#   for (j in 1:11) { #ageclass
#     for (k in 1:3){ # birthsites
#       survival[i,j,k] <- exp(int + weight.beta*weight.sim[i] + ageclass.beta[j] + bs.beta[k] + age.weight.beta[j]*weight.sim[i]) /
#                                 (1+exp(int + weight.beta*weight.sim[i] + ageclass.beta[j] + bs.beta[k] + age.weight.beta[j]*weight.sim[i]))
#     }
# }
# }
#       
  
}
",fill = TRUE)
sink()



# Bundle data
jags.data <- list(ch = ch, f = f, nind = nrow(ch), nocc = ncol(ch), weight = weight, bs = bs, #capyear = capyear,ageclass = ageclass
                  occasions=occasions, NA_indices=NA_indices)#, weight.sim = weight.sim

# Initial values
inits <- function(){list(weight = weight.init, weight.beta = rnorm(1,0,1), z=known.state.cjs(ch), bs.weight.beta = c(NA, rnorm(2,0,1)), #eps.capyear = c(NA, rnorm(14,0,1)),
                         bs.beta = c(NA, rnorm(2,0,1)), int = rnorm(1,0,1))} #, 

parameters <- c('int', 'bs.beta', 'weight.beta', 'bs.weight.beta')#, 'eps.capyear', 

# MCMC settings
ni <- 60000
nt <- 10
nb <- 30000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

print(cjs.weight)
MCMCtrace(cjs.weight)
# 
# 
# write.csv(cjs.weight$summary, 'weight.csv', row.names = T)

#########################################################################################
# 
# #create a tibble of the posterior draws
# gather<- cjs.weight %>% gather_draws(survival[weight,age, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$weight == 2)[1] # 108000 values of rain 1
# 
# #unscale and uncenter rain.sim
# weight.sim1 <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
# 
# #create vector containing simulated rainfall data but in the format to sync up with gather
# vector <- numeric(0)
# weight.sim2 <- for (i in weight.sim1) {
#   rep_i <- rep(i, times = 108000)
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$weight1 <- vector
# 
# #plot for ageclass 7
# 
# phi<-
#   subset(gather, age %in% '7') %>%
#   ggplot(aes(x=weight1, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
#   labs(x = "WEIGHT (POUNDS)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Survival by weight")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.15,0.1),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# ggsave('./figures/weight.jpg', phi, width = 15, height = 10)
# 

