#dummy weight survival model
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
dim(ch)

ch <- ch[c(1:3),c(1:15)]


ch1 <- matrix(c(1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,
                1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,
                1,1,1,1,1,0,1,1,0,0,0,0,0,0,0),
              nrow = 3, ncol = 15, byrow = T)

weight1 <- matrix(c(135,145,156,167,178,188,NA,NA,200,NA,NA,NA,NA,NA,NA,
                    124,134,145,156,167,178,189,199,199,187,NA,NA,NA,NA,NA,
                    122,134,134,145,166,NA,188,167,NA,NA,NA,NA,NA,NA,NA),
                  nrow = 3, ncol = 15, byrow = T)



# Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
#at birth, but we are starting at the second capture occasion
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch1, 1, get.first) 


#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'animal_id' )
weight<- as.matrix(weight[c(1:3),c(2:16)])


#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch1), ncol = ncol(ch1))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:ncol(ch1)] <- 1
  z.init[i,f[i]] <- NA
}

# 
#function for weight matrix
weight.init <- weight1
weight.init[is.na(weight.init)]<-mean(data$weight, na.rm = T)
weight.init[!is.na(weight1)]<-NA

for (i in 1:dim(weight.init)[1]){
  weight.init[,f[i]]<- NA
}



occasions <- rowSums(is.na(weight1)) # number of NA occasions for individual

# # weight <- as.data.frame(weight)
# weight <- as.data.frame(weight)
# indices <- as.data.frame(which(is.na(weight), arr.ind=T))
# indices <- indices %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup() #arrange orders the "rows" , group by may not be necessary, mutate creates a new column inserting the number for that individual 
NA_indices <- matrix(c(7,8,10,11,12,13,14,15,NA, NA,NA,NA,NA,NA,NA,
                       12,13,14,15,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       6,9,10,11,12,13,14,15,NA,NA,NA,NA,NA,NA,NA), 
                     nrow=3, ncol=15, byrow = T) #rows of individuals, columns of capture occasions 
# for(i in 1:nrow(indices)){
#   NA_indices[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
# }
# 
# weight <- as.data.frame(weight)
# indices <- as.data.frame(which(is.na(weight), arr.ind=T))
# indices <- indices %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
# NA_indices <- matrix(NA, nrow=5, ncol=15)
# for(i in 1:40){
#   NA_indices[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
# }


# 
# NA_indices <- as.data.frame(which(is.na(weight), arr.ind = TRUE))
# na_matrix <- matrix(NA, nrow = 493, ncol = 15, byrow = TRUE)
# 
# NA_indices1<- pivot_wider(NA_indices, names_from = 'col', values_from = 'col', id_cols = 'row' )
# NA_indices1 <- NA_indices1[order(NA_indices1$row),]
# 
# NA_indices <- NA_indices[,-1]
# NA_indices <- as.matrix(NA_indices)

# 
# known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest
# 
# #create capture history with just 1s and 0s, remove 'known fates'
# indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates
# ch[indices] <- 1
# 
# 
# # 
# # create vector with last occasion for each individual, marked by 2, 15 for end of study
# # rework h to only include capture myopathy or harvest, do not censor natural mortality
# get.last<- function(x) min(which(x>1))
# h <- apply(known.fate,1,get.last)
# h <- replace(h, is.infinite(h), 15)
# h
# f-h #check for zero

#####################################33
set.seed(100)
sink("cjs-weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)

#priors

for (u in 1:nind){                 #for u individuals 1 to 489
  for (j in 1:occasions[u]){  #for number of NA occasions for individual u
  weight[u,NA_indices[u,j]] ~ dnorm(0,0.0001)   #in weight.init, row u, column with NA, this code is just trying to find the column with NAs
     }
}


weight.beta ~ dnorm(0,0.0001)


tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):nocc){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            logit(phi[i,t-1]) <- weight.beta*weight[i,t-1] 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
                                            

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

      
  
}
",fill = TRUE)
sink()



# Bundle data
jags.data <- list(ch = ch1, f = f,  weight = weight1, nind = nrow(ch1), nocc = ncol(ch1), 
                  occasions=occasions, NA_indices=NA_indices)#h = h,capyear=capyear, birthyear = birthyear  bs = bs,weight.sim = weight.sim,ageclass = ageclass

# Initial values
inits <- function(){list(weight = weight.init, weight.beta = rnorm(1,0,10), z=z.init)} #, z.init

parameters <- c('weight.beta')#'int', 'age.beta', , 'p','survival'

# MCMC settings
ni <- 3000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

print(cjs.weight)
