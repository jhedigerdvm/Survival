#code to look at senescence, NO birthsite effects , will not converge derived parameters 

library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#load data 
data<-read.csv('./cleaned/caphx2022.csv', header = T)

#rename column names to year
data<-data %>% rename('2007' = X2007, '2008' = X2008, '2009' = X2009, '2010' = X2010, '2011' = X2011, '2012' = X2012, 
                      '2013' = X2013, '2014' = X2014, '2015' = X2015, '2016' = X2016, '2017' = X2017, 
                      '2018' = X2018, '2019' = X2019, '2020' = X2020, '2021' = X2021, '2022' = X2022)

CH<-data
known.fate<-data  #known deaths marked with 2
known.fate <- known.fate[, -c(1,18,19)]
known.fate<-as.matrix(known.fate) # 

#dataframe with only capture histories
CH <- CH[, -c(1,18,19)]
CH<-as.matrix(CH) # this will become y (i.e. our observations)

#create capture history with just 1s and 0s
indices <- which(CH == 2, arr.ind = TRUE)
CH[indices] <- 1

#create a vector with 1 for treatment, 2 for control, 3 for tgt
# bs<- as.numeric(factor(data$bs))
# unique(bs) 
# bs  #dmp is 1, e yana 2, w yana 3

#make age matrix, age class 1 is a 1.5 year old, age class 2 is a 2.5 year old 
ageclass<-data
for (i in 1:dim(ageclass)[1]){
  ageclass[i,2] <- 2007 - data$birth_year[i]  +1
  ageclass[i,3] <- 2007 - data$birth_year[i]  +2
  ageclass[i,4] <- 2007 - data$birth_year[i]  +3
  ageclass[i,5] <- 2007 - data$birth_year[i]  +4
  ageclass[i,6] <- 2007 - data$birth_year[i] +5
  ageclass[i,7] <- 2007 - data$birth_year[i] +6
  ageclass[i,8] <- 2007 - data$birth_year[i] +7
  ageclass[i,9] <- 2007 - data$birth_year[i] +8
  ageclass[i,10] <- 2007 - data$birth_year[i] +9
  ageclass[i,11] <- 2007 - data$birth_year[i] +10
  ageclass[i,12] <- 2007 - data$birth_year[i] +11
  ageclass[i,13] <- 2007 - data$birth_year[i] +12
  ageclass[i,14] <- 2007 - data$birth_year[i] +13
  ageclass[i,15] <- 2007 - data$birth_year[i] +14
  ageclass[i,16] <- 2007 - data$birth_year[i] +15
  ageclass[i,17] <- 2007 - data$birth_year[i] +16}


ageclass<- ageclass[,-c(1,19, 18)]


#Make any zeros or negatives NA
ageclass[ageclass<=0] <- NA
ageclass<-as.matrix(ageclass) 


# Create vector with the occasion each indiv is marked
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first) 
f #could use this to evaluate cohort effects by year

#create vector with last occasion for each individual, marked by 2, 15 for end of study 
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 15)
h

#create capture year vector
#make age matrix with 1 as fawn, 2 as immature, 3 as mature
capyear<-f
capyear


# Specify model in JAGS language, model runs and converges without derived parameters, but does not with the derived param
sink("cjs-age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
age.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate

for (u in 2:15){     #16 age classes but 15 survival periods                        
   age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
}

int~dnorm(0,0.001)



for(capyear in 1:16){
  eps.capyear[capyear] ~ dnorm(0, tau.capyear)
  }
  tau.capyear <- 1/(sigma.capyear * sigma.capyear)
  sigma.capyear ~ dunif(0,5) #random effect SD youre saying that each year-specific error term is coming
                                # from the same distributionm i.e. same mean (0) and standard
                                # deviation (sigma.period). The standard deviation is what connects all the
                                # year-specific random effects to the same distribution.

   
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
   z[i,f[i]] <- 1

      for (t in (f[i]+1):n.occasions){ #could put 15 for number of occasions, f[i]+1 because we condition on the first capture
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]] + eps.capyear[capyear[t-1]]

        # Observation process
            CH[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameter
  for (j in 1:15){ #age class
    survival[j] <- exp(int + age.beta[j] )/ (1 + exp(int + age.beta[j]))
  }                     #delta method to convert from logit back to probability Powell et al. 2007

# # # # for (i in 2:3){
# # #   site_diff[i] <- site.beta[1] - site.beta[i]
# # # }


  

}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(CH), ncol = ncol(CH))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:dim(z.init)[2]] <- 1
  z.init[i,f[i]] <- NA
  # z.init[i,h[i]:dim(z.init)[2]]<- NA #update z.init to account for known death or end of study 
  # z.init[i,h[i]]<- 1
}


# Bundle data
jags.data <- list(CH = CH, f = f, nind = nrow(CH), ageclass = ageclass, n.occasions = ncol(CH), capyear=capyear)# ,h = h,bs = bs , 


# Initial values
inits <- function(){list(int= rnorm(1,0,1), z = z.init, age.beta = c(NA, rnorm(14,0,1)), eps.capyear = rnorm(16,0,1), sigma.capyear = runif(1,0,5))}
#,site.beta = c(NA, rnorm(2,0,1)),int   

parameters <- c('age.beta', 'p', 'survival' )
# 'survival_age1_site1','survival_age2_site1', 'survival_age3_site1', 'survival_age4_site1','site_diff'
# 'survival_age1_site2','survival_age2_site2', 'survival_age3_site2','survival_age4_site2',
# 'survival_age1_site3','survival_age2_site3', 'survival_age3_site3', 'survival_age4_site3') 'survival',, 'site_diff','eps.capyear', 'survival'

# MCMC settings
ni <- 3000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age <- jagsUI(jags.data, inits, parameters, "cjs-age.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
MCMCtrace(cjs.age)
print(cjs.age)
summary<- cjs.age$summary
write.csv(summary,'surv_age_output1.csv')

#look into cohort effects, we are able to follow survival of a deer over its lifetime
#whats the effect of body mass and antlers on survival , include all animals including DMP
#fawn is 4 months to 1.5 yo 

#create a tibble of the posterior draws
posterior<- tidy_draws(cjs.age)
posterior<- posterior[,c(20:34)]

#create dataframe with posteriors of just survival age1 across the three sites
#pivot longer puts them in a tibble format
posterior_long <- posterior %>% pivot_longer(everything())
#make a new column 'site', rep 1:3 assigns 1 to site 1, 2 to site 2 for the number of rows divided by 3
posterior_long$ageclass <- rep(c('1', '2', '3','4','5','6','7','8','9','10','11','12','13','14','15'), nrow(posterior_long)/15)

##GGPLOT
plot_base_phi <- 
  ggplot(data = posterior_long, aes(x=ageclass, y=value))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.9,0.93),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        plot.title = element_text(face = 'bold', size = 40, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

phi.plot<- plot_base_phi +
  stat_pointinterval(alpha = 1, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  # scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
  #                    values=c("orange", "black", "darkorchid"))+
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                   labels = c('0.5-1.5', '1.5-2.5', '2.5-3.5', '3.5-4.5' ,'4.5-5.5' ,'5.5-6.5' ,'6.5-7.5','7.5-8.5',
                              '8.5-9.5','9.5-10.5','10.5-11.5', '11.5-12.5'))+
  labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", 
       title = "ANNUAL SURVIVAL BY AGE CLASS")

ggsave('phi_age.png', phi.plot, bg='transparent', width = 15, height = 10)
ggsave('phi_age.jpg', phi.plot, width = 15, height = 10)

