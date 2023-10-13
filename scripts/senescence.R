#survival model looking at senescence 

library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

# 
# 
# #load data 
# data<-read.csv('./cleaned/capture_cleaned_nofawns_wide.csv', header = T)
# 
# #rename column names to year
# data<-data %>% rename('2008' = X2008, '2009' = X2009, '2010' = X2010, '2011' = X2011, '2012' = X2012, 
#                       '2013' = X2013, '2014' = X2014, '2015' = X2015, '2016' = X2016, '2017' = X2017, 
#                       '2018' = X2018, '2019' = X2019, '2020' = X2020, '2021' = X2021, '2022' = X2022)
# 
# CH<-data
# CH<- CH[,-c(1:3)] #removing animal id, birth site, and birth year
# CH<-as.matrix(CH) # this will become y (i.e. our observations)
# 
# 
# known.fate<-data  #known deaths marked with 2
# known.fate <- known.fate[, -c(1:3)]
# known.fate<-as.matrix(known.fate) 
# 
# #create capture history with just 1s and 0s
# indices <- which(CH == 2, arr.ind = TRUE)
# CH[indices] <- 1 #replaces all 2s with a 1
# 
# #create a vector with 1 for treatment, 2 for control, 3 for tgt
# bs<- as.numeric(factor(data$bs))
# unique(bs) 
# 
# id <- as.numeric(factor(data$animal_id))
# 
# #make age matrix, age class 1 is a 1.5 year old, age class 2 is a 2.5 year old 
# ageclass<-data[,-c(1,3)]
# for (i in 1:dim(ageclass)[1]){
#   ageclass[i,2] <- 2007 - data$birth_year[i]  +1
#   ageclass[i,3] <- 2007 - data$birth_year[i]  +2
#   ageclass[i,4] <- 2007 - data$birth_year[i]  +3
#   ageclass[i,5] <- 2007 - data$birth_year[i]  +4
#   ageclass[i,6] <- 2007 - data$birth_year[i] +5
#   ageclass[i,7] <- 2007 - data$birth_year[i] +6
#   ageclass[i,8] <- 2007 - data$birth_year[i] +7
#   ageclass[i,9] <- 2007 - data$birth_year[i] +8
#   ageclass[i,10] <- 2007 - data$birth_year[i] +9
#   ageclass[i,11] <- 2007 - data$birth_year[i] +10
#   ageclass[i,12] <- 2007 - data$birth_year[i] +11
#   ageclass[i,13] <- 2007 - data$birth_year[i] +12
#   ageclass[i,14] <- 2007 - data$birth_year[i] +13
#   ageclass[i,15] <- 2007 - data$birth_year[i] +14
#   ageclass[i,16] <- 2007 - data$birth_year[i] +15
#   # ageclass[i,17] <- 2007 - data$birth_year[i] +16
#   
#   }
# 
# ageclass<- ageclass[,-c(1)]
# 
# #Make any zeros or negatives NA
# ageclass[ageclass<=0] <- NA
# ageclass<-as.matrix(ageclass) 
# 
# # Create vector with the occasion each indiv is marked, this gets weird because we know each individual was caught
# #at birth, but we are starting at the second capture occasion
# get.first <- function(x) min(which(x!=0))
# f <- apply(ageclass, 1, get.first) 
# f #could use this to evaluate cohort effects by year
# 
# #create vector with last occasion for each individual, marked by 2, 15 for end of study 
# #rework h to only include capture myopathy or harvest, do not censor natural mortality 
# get.last<- function(x) min(which(x>1))
# h <- apply(known.fate,1,get.last)
# h <- replace(h, is.infinite(h), 15)
# h
# f-h #make sure there are no first and last capture occasions that overlap
# 
# # Assuming you have two matrices named matrix1 and matrix2
# # # Subtract matrix2 from matrix1
# result_matrix <- f - h
# 
# # Create a logical matrix indicating where result_matrix is zero
# zero_indices <- result_matrix == 0
# 
# # Use the logical matrix to subset matrix1 (or matrix2, they should be equivalent)
# f <- f[!zero_indices]
# h<-h[!zero_indices]
# 
# data1<-cbind(zero_indices, data)
# data1<- data1 %>% filter(data1$zero_indices=='FALSE')
# data<-data1
# CH <- data[,-c(1:4)]
# CH<-as.matrix(CH) # this will become y (i.e. our observations)
# # 
# # known.fate<-data  #known deaths marked with 2
# # known.fate <- known.fate[, -c(1:3)]
# # known.fate<-as.matrix(known.fate) 
# 
# #create capture history with just 1s and 0s
# indices <- which(CH == 2, arr.ind = TRUE)
# CH[indices] <- 1 #replaces all 2s with a 1
# 
# #create a vector with 1 for treatment, 2 for control, 3 for tgt
# bs<- as.numeric(factor(data$bs))
# unique(bs) 
# 
# id <- as.numeric(factor(data$animal_id))
# 
# #make age matrix, age class 1 is a 1.5 year old, age class 2 is a 2.5 year old 
# ageclass<-data[,-c(1,2,4)]
# for (i in 1:dim(ageclass)[1]){
#   ageclass[i,2] <- 2007 - data$birth_year[i]  +1
#   ageclass[i,3] <- 2007 - data$birth_year[i]  +2
#   ageclass[i,4] <- 2007 - data$birth_year[i]  +3
#   ageclass[i,5] <- 2007 - data$birth_year[i]  +4
#   ageclass[i,6] <- 2007 - data$birth_year[i] +5
#   ageclass[i,7] <- 2007 - data$birth_year[i] +6
#   ageclass[i,8] <- 2007 - data$birth_year[i] +7
#   ageclass[i,9] <- 2007 - data$birth_year[i] +8
#   ageclass[i,10] <- 2007 - data$birth_year[i] +9
#   ageclass[i,11] <- 2007 - data$birth_year[i] +10
#   ageclass[i,12] <- 2007 - data$birth_year[i] +11
#   ageclass[i,13] <- 2007 - data$birth_year[i] +12
#   ageclass[i,14] <- 2007 - data$birth_year[i] +13
#   ageclass[i,15] <- 2007 - data$birth_year[i] +14
#   ageclass[i,16] <- 2007 - data$birth_year[i] +15
#   # ageclass[i,17] <- 2007 - data$birth_year[i] +16
#   
# }
# 
# ageclass<- ageclass[,-c(1)]
# 
# #Make any zeros or negatives NA
# ageclass[ageclass<=0] <- NA
# ageclass<-as.matrix(ageclass) 
# 
# f-h #make sure there are no first and last capture occasions that overlap
# 
# #create capture year vector
# #make age matrix with 1 as fawn, 2 as immature, 3 as mature
# capyear<-f
# capyear


# Specify model in JAGS language
sink("cjs-senescence.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
age.beta[1] <- 0
eps.bs[1] <- 0
eps.capyear[1] <- 0
eps.id[1] <- 0

int ~ dnorm(0, 0.001)

for (u in 2:14){     #15 age classes but 14 survival periods                        
   age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
}

for (u in 2:15){      #capyear
  eps.capyear[u] ~ dnorm(0, tau.capyear)
}

  sigma ~ dunif(0,10) #standard dev
  tau.capyear <- 1/(sigma*sigma) #precision = 1/variance

for (u in 2:3){ #birth site
  eps.bs[u] ~ dnorm(0,tau.bs)
}
  
  tau.bs <- 1/(sigma*sigma)

for (u in 2:490){ #animal_id
  eps.id[u] ~ dnorm(0, tau.id)
}
  tau.id <-1/(sigma*sigma)

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + eps.capyear[capyear[t-1]] + eps.bs[bs[i]] + eps.id[id[i]]
        # Observation process
            CH[i,t] ~ dbern(mu2[i,t])
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
z.init <- matrix(NA, nrow = nrow(CH), ncol = ncol(CH))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}


# Bundle data
jags.data <- list(h = h, CH = CH, f = f, nind = nrow(CH),  ageclass = ageclass, capyear = capyear, bs = bs, id = id)#, , 

# Initial values
inits <- function(){list(int = rnorm(1, 0, 1), z = z.init, age.beta = c(NA, rnorm(13,0,1)),
                          sigma = runif(1,0,10), eps.capyear = c(NA, rnorm(14,0,1)),
                         eps.bs = c(NA, rnorm(2,0,1)), eps.id = c(NA, rnorm(489, 0, 1)))} #,, 

parameters <- c('int', 'age.beta', 'p', 'survival', 'surv_diff')
# 'survival', 'site_diff' 'survival',, 'site_diff','eps.capyear' 'int','site.beta',

# MCMC settings
ni <- 15000
nt <- 1
nb <- 10000
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
