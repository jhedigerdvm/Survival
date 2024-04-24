#dummy weight survival model
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

data<- read.csv('./cleaned/dummy.weight.csv', header = T)

#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'year', values_from = 'status', id_cols = 'id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 


#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'id' )
weight<- as.matrix(weight[,-1])


#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:7] <- 1
  z.init[i,f[i]] <- NA
}
# 
#function for weight matrix
weight.init <- weight
weight.init[is.na(weight.init)]<-mean(data$weight, na.rm = T)
weight.init[!is.na(weight)]<-NA




#####################################33
set.seed(100)
sink("cjs-weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)

#priors

for (j in 1:nocc){
 weight.beta[j] ~ dnorm(0,0.0001)
  for (u in 1:nind){  #prior for missing weights
     weight[u,j] ~ dnorm(0,0.0001)
     }
}


tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood
for (i in 1:nind){

  z[i,f[i]] <- 1      # Define latent state at first capture, we know for sure the animal is alive

      for (t in (f[i]+1):nocc){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- weight.beta[ t-1 ] * weight[ i, t-1 ]
                                                           

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

      
  
}
",fill = TRUE)
sink()



# Bundle data
jags.data <- list(ch = ch, f = f, nind = nrow(ch), nocc = ncol(ch), weight = weight)#capyear=capyear, birthyear = birthyear  bs = bs,weight.sim = weight.sim,ageclass = ageclass

# Initial values
inits <- function(){list(z = z.init, weight = weight.init, weight.beta = rnorm(7,0,1))} 

parameters <- c('weight.beta')#'int', 'age.beta', , 'p','survival'

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.weight)
MCMCtrace(cjs.weight)
