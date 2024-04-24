#Survival as a function of weight and age

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

#create ageclass matrix
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

#create ageclass with ages 1-5 as 1, 6-9 as 2, 10-11 as 3
ageclass2<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass2 <- ageclass2[,-1]

ageclass2[] <- ifelse(ageclass2[] >= 1 & ageclass2[] <= 5, 1, 
                      ifelse(ageclass2[] >= 6 & ageclass2[] <= 9, 2,
                             3))
ageclass2 <- as.matrix(ageclass2)

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

#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'animal_id' )
weight<- as.matrix(weight[,-1])

antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcsin', id_cols = 'animal_id' )
antlers<- as.matrix(antlers[,-1])

#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

#survival as a function of weight, age and the interaction
set.seed(100)
sink("cjs-weight.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors

for (u in 1:14){
  age.beta[u] ~ dnorm(0,0.01)
}

weight.beta ~ dnorm(0,0.001)

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
            logit(phi[i,t-1]) <- age.beta[ageclass[i,t-1]] + weight.beta*weight[i,t-1] 
                                            

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i
# 
# #derived parameters
#   for (j in 1:100){ #weight sim 
#     for (i in 1:10){ #ageclass
#       survival[j,i] <- exp(weight.beta*weight.sim[j] + age.beta[i])/
#                             (1 + exp(weight.beta*weight.sim[j] + age.beta[i]))
#       }                     #delta method to convert from logit back to probability Powell et al. 2007
#   }
      
  
}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}

#function for weight matrix
weight.init <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))

for (i in 1:nrow(weight.init)) {
  for (j in 1:ncol(weight.init)) {
    if (is.na(weight[i, j])) {
      weight.init[i, j] <- 1
    }
  }
}

# Bundle data
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), weight = weight, 
                  ageclass = ageclass)#capyear=capyear, birthyear = birthyear  bs = bs,weight.sim = weight.sim,

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, weight=weight.init, 
                         age.beta = rnorm(14,0,1),
                         weight.beta = rnorm(1,0,1))} 
#eps.capyear = c(NA, rnorm(13,0,1)), eps.birthyear = c(NA, rnorm(13,0,1age.site.beta = c(NA, rnorm(13,0,1)),site.beta = c(NA, rnorm(2,0,1)), 

parameters <- c('int', 'age.beta', 'weight.beta', 'p','survival' )

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.weight <- jagsUI(jags.data, inits, parameters, "cjs-weight.jags", n.chains = nc,
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.weight)
# write.csv(cjs.age.site$summary, "./output/age.site2.csv")
# MCMCtrace(cjs.age.site)
# 
# write.csv(cjs.rain.site.age$summary, './output/rain.site.age.csv')
# # 
#create a tibble of the posterior draws
# gather<- cjs.age.site %>% gather_draws(survival[site,age]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)
# 
# 
# plot2<- gather %>% 
#   ggplot(aes(x=age, y=.value, color = site, fill = site)) +
#   stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
#                      position = position_dodge(width = 0.5)) +
#   scale_fill_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
#   scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
#   
#   # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
#   #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
#   # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
#   #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
#   labs(x = "AGE", y = "SURVIVAL (%)", title = "SURVIVAL BY AGE AND SITE")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.5,0.3),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 24),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# plot2
# ggsave('./figures/phi.age.site.jpg', plot2, width = 10, height = 6)

