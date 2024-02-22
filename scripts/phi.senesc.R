#phi senescence
#Survival model with effect of birth site and age on survival, excludes fawn to 1.5 year old

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
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

#create birthsite vector
id.bs.by <- unique(data[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass matrix
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

# #create animal id vector
# id <- as.numeric(factor(id.bs.by$animal_id))

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


#model eval the effect of site and age on surv and the interaction between the two
# Specify model in JAGS language
set.seed(100)
sink("cjs-senesc.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

for (u in 1:14){
  age.beta[u] ~ dnorm(0,0.01)
}

for (u in 1:14){
  eps.capyear[u] ~ dnorm(0, sigma)
}

for (u in 1:14){
  eps.birthyear[u] ~ dnorm(0, sigma)
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
            logit(phi[i,t-1]) <-  age.beta[ageclass[i,t-1]] + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
    for (j in 1:10){ #ageclass
      survival[j] <- exp(age.beta[j])/
                            (1 + exp(age.beta[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
  
      
#   for (i in 1:3){
#     for (j in 1:10){
#     survival_diff[i,j] <- survival[1,j] - survival[i,j]  
#     }
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, ageclass = ageclass, capyear = capyear,
                  birthyear = birthyear)

# Initial values
inits <- function(){list(z = z.init, eps.capyear = rnorm(14,0,1), eps.birthyear = rnorm(14,0,1),
                         age.beta = rnorm(14,0,1))} #

parameters <- c('age.beta','eps.capyear', 'eps.birthyear', 'p','survival' )

# MCMC settings
ni <- 1500
nt <- 1
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.senesc <- jagsUI(jags.data, inits, parameters, "cjs-senesc.jags", n.chains = nc,
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.senesc)
# 
# write.csv(cjs.rain.site.age$summary, './output/rain.site.age.csv')
# # 
#create a tibble of the posterior draws
gather<- cjs.senesc %>% gather_draws(survival[age]) #this creates a dataframe in long format with indexing
gather$age <- as.factor(gather$age)


plot2<- gather %>% 
  ggplot(aes(x=age, y=.value)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  labs(x = "AGE", y = "SURVIVAL (%)", title = "SURVIVAL BY AGE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/phi_senesc.jpg', plot2, width = 10, height = 7)


plot2<- gather %>% filter(age== c("7", "8", "9", "10")) %>% #filter(age== c("1", "2", "3", "4", "5")) %>% filter(age== c("5", "6", "7")) %>% #
  ggplot(aes(x=age, y=.value)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  labs(x = "AGE", y = "SURVIVAL (%)", title = "SURVIVAL BY AGE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/phi_senesc78910.jpg', plot2, width = 10, height = 7)


# #find first row for 2nd rain value
# first_idx <- which(gather$rain == 2)[1] # 49500 values of rain 1 
# 
# #unscale and uncenter rain.sim
# # rain.sim1 <- (rain.sim * sd(data$annual)) + mean(data$annual)
# 
# #create vector containing simulated rainfall data but in the format to sync up with gather 
# vector <- numeric(0)
# rain.sim2 <- for (i in rain.sim1) {
#   rep_i <- rep(i, times = 99000)
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$rain1 <- vector
# 
# #plot for all ages and facet wrap
# phi.plot<- (gather %>%  
#               ggplot(aes(x=rain1, y=.value, color = site, fill = site)) + #color equals line, fill equals ribbon
#               facet_wrap(vars(age))+
#               stat_lineribbon(.width = 0.95, alpha = 1/4)) #.width is the CRI, alpha is opacity
# 
# phi.plot.1<- 
#   subset(gather, age %in% '1') %>%  
#   ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI 
#   scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d()+ #color of line but no opacification
#   labs(x = "RAINFALL (in)", y = "ANNUAL SURVIVAL PROBABILITY", title = "Phi ~ int + age(1-2) + site + rain + rain*site")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.5,0.3),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)