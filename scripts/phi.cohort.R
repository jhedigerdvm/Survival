#phi cohorts
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
sink("cjs-cohort.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

int ~ dnorm(0,0.001)
age.beta[1] <- 0
birthyear.beta[1] <- 0
# eps.capyear[1] <- 0

for (u in 2:14){
  age.beta[u] ~ dnorm(0,0.01)
}
# 
# for (u in 2:14){
#   eps.capyear[u] ~ dnorm(0, sigma)
# }

for (u in 2:14){
  birthyear.beta[u] ~ dnorm(0, 0.001)
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
            logit(phi[i,t-1]) <-  int + age.beta[ageclass[i,t-1]]  + birthyear.beta[birthyear[i]] 
                                          #+ eps.capyear[capyear[i]]

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
    for (j in 1:10){ #ageclass
      for (i in 1:14){ #birthyear
      survival[j, i] <- exp(int + age.beta[j] + birthyear.beta[i])/
                            (1 + exp(int + age.beta[j]+ birthyear.beta[i]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
    }
      
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
inits <- function(){list(int = rnorm(1,0,1), z = z.init, eps.capyear = c(NA, rnorm(13,0,1)), 
                          birthyear.beta = c(NA, dnorm(13,0,1)), age.beta = c(NA, rnorm(13,0,1)))} #

parameters <- c('int','age.beta','eps.capyear', 'birthyear.beta', 'p','survival' )

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.cohort<- jagsUI(jags.data, inits, parameters, "cjs-cohort.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.cohort)
MCMCtrace(cjs.cohort)
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
