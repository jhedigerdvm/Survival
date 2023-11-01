#Survival model for age and birthsite, capyear as a random effect
#treatment and control survival are not different, but TGT survive less 
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

data<- read.csv('./cleaned/caphx.rainfall.long.csv', header = T)
data<- data[data$birth_year != '2021',] #remove 2021 cohort

#remove individuals with no capture ocassions
data1<-
  data %>% 
  group_by(animal_id) %>%
  summarise(status_sum = sum(status)) %>%
  filter(status_sum < 1) -> to_remove

data2 <- data[!data$animal_id %in% to_remove$animal_id,]

ch<- pivot_wider(data2, names_from = 'year', values_from = 'status', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)
data2$annual <-scale(data2$annual)
annual.rainfall<-pivot_wider(data2, names_from = 'year', values_from = 'annual', id_cols = 'animal_id' )
annual.rainfall<-annual.rainfall[,-1]
annual.rainfall<-as.matrix(annual.rainfall)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates 
ch[indices] <- 1


# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0))
f <- apply(ch, 1, get.first) 

#create birthsite vector
id.bs.by <- unique(data2[, c("animal_id", "bs",'birth_year')])
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = dmp, 2 = ey, 3 = wy

#create ageclass matrix
ageclass<- pivot_wider(data2, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
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
f-h #need to remove birth year 2021 



# Specify model in JAGS language
sink("cjs-age-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
# age.beta[1] <- 0
site.beta[1] <- 0
rain.beta[1] <- 0
int ~ dnorm(0, 0.001)

# for (u in 2:14){     #15 age classes but 14 survival periods                        
#    age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
# }

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

for (u in 2:489){
  rain.beta[u] ~ dnorm(0,0.01)
}


# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + site.beta[bs[i]] + rain.beta[rain[i]]
            
        # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i


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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall)#, 

# Initial values
inits <- function(){list(int = rnorm(1, 0, 1), z = z.init, rain.beta = c(NA, rnorm(488, 0,1)),
                                                              site.beta = c(NA, rnorm(2,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta', 'p')
# 'survival', 'site_diff' 'survival',, 'site_diff','eps.capyear' 'int','site.beta',

# MCMC settings
ni <- 1000
nt <- 10
nb <- 100
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age.site <- jagsUI(jags.data, inits, parameters, "cjs-age-site.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.age.site)
