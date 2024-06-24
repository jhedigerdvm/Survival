library(jagsUI)
library(tidyverse)

ch <- as.matrix(read.csv('./cleaned/caphist.csv', header = T))
weight <- as.matrix(read.csv('./cleaned/weight.csv', header = T))
weight.init <- as.matrix(read.csv('./cleaned/weight.init.csv', header = T))
bs <- ch[,16] #create birthsite vector
ch <- ch[,-16] #remove birthsite column from cap history

# Create vector with the occasion each indiv is marked, 
get.first <- function(x) min(which(x!=0)) 
f <- apply(ch, 1, get.first) 


#Function for latent state  
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


#####################################
#model looking at phi as a function of weight and the interaction of weight with site

set.seed(100)
sink("cjs-weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta( 1 , 1 )

#priors

for (u in 1:3) { #birth site and weight interaction
  bs.weight.beta[u] ~ dnorm(0, 0.001)
}

for (u in 1:nind){
  for (j in 1:occasions[u]){  #prior for missing weights
  weight[u,NA_indices[u,j]] ~ dnorm( 0 , 0.0001 )
     }
}

weight.beta ~ dnorm( 0 , 0.001 )

for (u in 1:3){                               #prior for birth site
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
            logit(phi[i,t-1]) <- weight.beta*weight[i,t-1] + bs.beta[bs[i]]+ bs.weight.beta[bs[i]]*weight[i,t-1]
                                       # int +   #+ eps.capyear[capyear[i]] + ageclass.beta[ageclass[i,t-1]] 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
                                            

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i
      
  
}
",fill = TRUE)
sink()



# Bundle data
jags.data <- list(ch = ch, f = f, nind = nrow(ch), nocc = ncol(ch), weight = weight, bs = bs, 
                  occasions=occasions, NA_indices=NA_indices)#

# Initial values
inits <- function(){list(weight = weight.init, weight.beta = rnorm(1,0,1), z=known.state.cjs(ch), 
                         bs.weight.beta = rnorm(3,0,1), 
                         bs.beta = rnorm(3,0,1) )} 

parameters <- c('weight.beta','bs.beta', 'bs.weight.beta' )

# MCMC settings
ni <- 20000
nt <- 10
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
# traceplot(cjs.weight)
print(cjs.weight)
MCMCtrace(cjs.weight)


