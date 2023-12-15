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


#model eval the effect of site and age on surv
# Specify model in JAGS language
set.seed(100)
sink("cjs-age-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
int ~ dnorm(0,0.01)

site.beta[1] <- 0
age.beta[1]<- 0

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

for (u in 2:14){
  age.beta[u] ~ dnorm(0,0.01)
}

# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + site.beta[bs[i]] + age.beta[ageclass[i,t-1]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site 
    for (j in 1:10){ #ageclass
      survival[i,j] <- exp(int + site.beta[i] + age.beta[j])/
                            (1 + exp(int + site.beta[i] + age.beta[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
  }
      
  for (i in 1:3){
    for (j in 1:10){
    survival_diff[i,j] <- survival[1,j] - survival[i,j]  
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, ageclass = ageclass)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, site.beta = c(NA, rnorm(2,0,1)), 
                         age.beta = c(NA, rnorm(13,0,1)))} #

parameters <- c('int', 'site.beta', 'age.beta', 'p','survival_diff','survival' )

# MCMC settings
ni <- 1500
nt <- 1
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age.site <- jagsUI(jags.data, inits, parameters, "cjs-age-site.jags", n.chains = nc,
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.age.site)
# 
# write.csv(cjs.rain.site.age$summary, './output/rain.site.age.csv')
# # 
# # #create a tibble of the posterior draws
# # gather<- cjs.rain.site.age %>% gather_draws(survival[rain,site,age]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# # gather$age <- as.factor(gather$age)
# 
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