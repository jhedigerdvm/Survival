#survival model looking at senescence 
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

# 
# annual.rainfall<-pivot_wider(data, names_from = 'year', values_from = 'annual.sc', id_cols = 'animal_id' )
# annual.rainfall<-annual.rainfall[,-1]
# annual.rainfall<-as.matrix(annual.rainfall)
# # write.csv(annual.rainfall, 'rain.jan.dec.csv', row.names = F)
# 
# nvalues <- 1000
# rain.sim <- seq(from = min(annual.rainfall, na.rm = T), to = max(annual.rainfall, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# rain.sim


# Specify model in JAGS language
sink("cjs-senescence.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
age.beta[1] <- 0
eps.capyear[1] <- 0
# eps.id[1] <- 0

int ~ dnorm(0, 0.001)

for (u in 2:14){     #15 age classes but 14 survival periods                        
   age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
}

for (u in 2:15){      #capyear
  eps.capyear[u] ~ dnorm(0, tau.capyear)
}

  sigma ~ dunif(0,10) #standard dev
  tau.capyear <- 1/(sigma*sigma) #precision = 1/variance
# 
# for (u in 2:489){ #animal_id
#   eps.id[u] ~ dnorm(0, 0.01)
# }
#   tau.id <-1/(sigma*sigma)

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + eps.capyear[capyear[t-1]] #+ eps.id[id[i]]
        # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i


#derived parameter
    for (j in 1:11){ #age class
      survival[j] <- exp(int+ age.beta[j] )/ (1 + exp(int+ age.beta[j]))
    }                     #delta method to convert from logit back to probability Powell et al. 2007
  
    for (j in 1:11){
      surv_diff[j] <- survival[j] - survival[6]
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  ageclass = ageclass, capyear = capyear, id = id)#, , 

# Initial values
inits <- function(){list(int = rnorm(1, 0, 1), z = z.init, age.beta = c(NA, rnorm(13,0,1)),
                          sigma = runif(1,0,10), eps.capyear = c(NA, rnorm(14,0,1))
                        )} #,,  eps.id = c(NA, rnorm(488, 0, 1))

parameters <- c('int', 'age.beta', 'p', 'survival', 'surv_diff')
# 'survival', 'site_diff' 'survival',, 'site_diff','eps.capyear' 'int','site.beta',

# MCMC settings
ni <- 4000
nt <- 10
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.senescence <- jagsUI(jags.data, inits, parameters, "cjs-senescence.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
MCMCtrace(cjs.senescence)
print(cjs.senescence)

# 
# summary<- cjs.age.site$summary
# write.csv(summary,'surv_output_site_age_capyear.csv')

# 
# 
# #create a tibble of the posterior draws
# posterior<- tidy_draws(cjs.age.site)
# posterior<- posterior[,-c(1:22)]
# posterior<- posterior[,c(1:33)]
# 
# #create dataframe with posteriors of just survival age1 across the three sites
# #pivot longer puts them in a tibble format
# posterior_long <- posterior %>% pivot_longer(everything())
# #make a new column 'site', rep 1:3 assigns 1 to site 1, 2 to site 2 for the number of rows divided by 3
# posterior_long$ageclass <- rep(c('1', '1', '1','2','2','2','3','3','3','4','4','4','5','5','5','6','6','6','7','7','7',
#                                  '8','8','8','9','9','9','10','10','10','11','11','11'), nrow(posterior_long)/33)
# 
# posterior_long$site<- rep(c('1','2','3'), nrow(posterior_long)/3)
# 
# ##GGPLOT
# plot_base_phi <- 
#   ggplot(data = posterior_long, aes(x=ageclass, y=value, group = site))+
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.9,0.93),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 24),
#         plot.title = element_text(face = 'bold', size = 40, hjust = 0.5 ),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 24),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# 
# phi.plot<- plot_base_phi +
#   stat_pointinterval(aes(color = site), alpha = 1, .width = c(0.5, 0.95), 
#                      position = position_dodge(width = 0.5)) +
#   scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
#                      values=c("orange", "black", "darkorchid"))+
#   scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'), 
#                    labels = c('1.5-2.5', '2.5-3.5', '3.5-4.5' ,'4.5-5.5' ,'5.5-6.5' ,'6.5-7.5','7.5-8.5',
#                               '8.5-9.5','9.5-10.5','10.5-11.5', '11.5-12.5'))+
#   labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", 
#        title = "ANNUAL SURVIVAL BY AGE CLASS")
# 
# ggsave('phi_age_site.png', phi.plot, bg='transparent', width = 15, height = 10)
# ggsave('phi_age_site.jpg', phi.plot, width = 15, height = 10)
