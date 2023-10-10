library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#load data 
data<-read.csv('./cleaned/data_cleaned.csv', header = T)

#rename column names to year
data<-data %>% rename('2007' = X2007, '2008' = X2008, '2009' = X2009, '2010' = X2010, '2011' = X2011, '2012' = X2012, 
                        '2013' = X2013, '2014' = X2014, '2015' = X2015, '2016' = X2016, '2017' = X2017, 
                        '2018' = X2018, '2019' = X2019, '2020' = X2020, '2021' = X2021)
# 
# data<-data[,-c(3)]
# data$sum<- rowSums(data[,c(3:16)])
# data2 <- subset(data, sum >1) 
# data2<- data2[,-c(17)]
# rownames(data2)<-NULL

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$bs))
unique(bs) 
bs  #dmp is 1, e yana 2, w yana 3

#dataframe with only capture histories
CH <- data[, c(4:5)]
CH<-as.matrix(CH) # this will become y (i.e. our observations)

#make age matrix
ageclass<-data
for (i in 1:nrow(ageclass)){
  ageclass[i,3] <- 2007 - ageclass$birth_year[i]+1
  
}
for (i in 1:dim(ageclass)[1]){
  ageclass[i,3] <- 2007 - ageclass$birth_year[i]  +2
  ageclass[i,4] <- 2007 - ageclass$birth_year[i]  +3
  ageclass[i,5] <- 2007 - ageclass$birth_year[i]  +4
  ageclass[i,6] <- 2007 - ageclass$birth_year[i]  +5
  ageclass[i,7] <- 2007 - ageclass$birth_year[i] +6
  ageclass[i,8] <- 2007 - ageclass$birth_year[i] +7
  ageclass[i,9] <- 2007 - ageclass$birth_year[i] +8
  ageclass[i,10] <- 2007 - ageclass$birth_year[i] +9
  ageclass[i,11] <- 2007 - ageclass$birth_year[i] +10
  ageclass[i,12] <- 2007 - ageclass$birth_year[i] +11
  ageclass[i,13] <- 2007 - ageclass$birth_year[i] +12
  ageclass[i,14] <- 2007 - ageclass$birth_year[i] +13
  ageclass[i,15] <- 2007 - ageclass$birth_year[i] +14
  ageclass[i,16] <- 2007 - ageclass$birth_year[i] +15
  ageclass[i,17] <- 2007 - ageclass$birth_year[i] +16
  
}

ageclass<- ageclass[,-c(1:2)]

#Make any zeros or negatives NA
ageclass[ageclass<=1] <- NA 
ageclass<- as.matrix(ageclass)
unique(ageclass)

# Create vector that states the column at which the first capture occasion occurs
get.first <- function(x) min(which(x!=0))
f <- apply(ageclass, 1, get.first) #change to second capture occasion CH
f #could use this to evaluate cohort effects by year

#found first capture occasion, now we want the second one
# for (i in 1:length(f)){
#   CH[i,f[i]]<- 0 #f[i] indicates the column of the first capture
# }

#Remove any rows with all zeroes
# CH<-CH[as.logical(rowSums(CH[,c(3:17)] != 0)), ]

# update f with new CH 
# get.first <- function(x) min(which(x!=0))
# f <- apply(CH, 1, get.first) #change to second capture occasion CH
# f #could use this to evaluate cohort effects by year

#create capture year vector
capyear<-f
capyear
# 
# ageclass<-ageclass[,-c(1)]
# CH <- CH[,-c(1)]

# Specify model in JAGS language
sink("cjs-age-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
age.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate 
site.beta[1] <- 0     #need to set baseline of zero for >1 categorical covariate

for (u in 2:13){   #dropped to 14 age groups from 15                          
   age.beta[u] ~ dnorm(0,0.01)              # Priors for age-specific survival
}

for (u in 2:3){                             
   site.beta[u] ~ dnorm(0, 0.01)               #priors for birthsite
}
int~dnorm(0,0.001)


# 
# for(capyear in 1:14){                               ##ERROR WITH DIMENSION MISMATCH
#   eps.capyear[capyear] ~ dnorm(0, tau.capyear)
#   }
#   tau.capyear <- 1/(sigma.capyear * sigma.capyear)
#   sigma.capyear ~ dunif(0,10) #random effect SD youre saying that each year-specific error term is coming
#                                 # from the same distributionm i.e. same mean (0) and standard
#                                 # deviation (sigma.period). The standard deviation is what connects all the
#                                 # year-specific random effects to the same distribution.

   
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
   z[i,f[i]] <- 1
   
      for (t in (f[i]+1):n.occasions){ #when adding in known death, for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t])
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
            logit(phi[i,t-1]) <- int + age.beta[ageclass[i,t-1]-1] + site.beta[bs[i]] #+ eps.capyear[capyear[t-1]]

        # Observation process
            CH[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameter
for (i in 1:3){ #bs
  for (j in 1:13){ #age class
    survival[i,j] <- exp(int + age.beta[j] + site.beta[i])/ (1 + exp(int + age.beta[j] + site.beta[i]))
  }
}
for (i in 1:2){
  site_diff[i] <- site.beta[1] - site.beta[i]
}


  

}
",fill = TRUE)
sink()


#Function for latent state
z.init <- matrix(NA, nrow = nrow(CH), ncol = ncol(CH))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:dim(z.init)[2]] <- 1
  z.init[i,f[i]] <- NA
}


# Bundle data
jags.data <- list(CH = CH, f = f, nind = nrow(CH), n.occasions = ncol(CH),  ageclass = ageclass, 
                 bs = bs)#capyear=capyear,eps.capyear = rnorm(14,0,1),  sigma.capyear = runif(1,0,10),


# Initial values
inits <- function(){list(int = rnorm(1,0,1),z = z.init, age.beta = c(NA, rnorm(12,0,1)),  site.beta = c(NA, rnorm(2,0,1)))}

parameters <- c(
  'int', 'age.beta', 'site.beta', 'p','survival', 'site_diff' )#, 'eps.capyear'
    # 'survival_age1_site1','survival_age2_site1', 'survival_age3_site1', 'survival_age4_site1','site_diff'
    # 'survival_age1_site2','survival_age2_site2', 'survival_age3_site2','survival_age4_site2',
    # 'survival_age1_site3','survival_age2_site3', 'survival_age3_site3', 'survival_age4_site3') 'survival',, 'site_diff'

# MCMC settings
ni <- 3000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age.site <- jagsUI(jags.data, inits, parameters, "cjs-age-site.jags", n.chains = nc, 
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
MCMCtrace(cjs.age.site)
print(cjs.age.site)
summary<- cjs.age.site$summary
write.csv(summary,'surv_output1.csv')

#create a tibble of the posterior draws
posterior<- tidy_draws(cjs.age.site)
posterior<- posterior[,-c(1:23)]
posterior<- posterior[,c(1:33)]

#create dataframe with posteriors of just survival age1 across the three sites
#pivot longer puts them in a tibble format
posterior_long <- posterior %>% pivot_longer(everything())
#make a new column 'site', rep 1:3 assigns 1 to site 1, 2 to site 2 for the number of rows divided by 3
posterior_long$bs <- rep(c('dmp','ey','wy'), nrow(posterior_long)/3)
# survival_res <- survival_res %>% rename_at('name', ~'ageclass')
posterior_long$ageclass <- rep(c('1','1','1','2','2','2','3','3','3','4','4','4',
                                 '5','5','5','6','6','6','7','7','7','8','8','8',
                                 '9','9','9','10','10','10','11','11','11'), nrow(posterior_long)/33)

##GGPLOT
plot_base_phi <- 
  ggplot(data = posterior_long, aes(x=ageclass, y=value, group = bs))+
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
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

phi.plot<- plot_base_phi +
  stat_pointinterval(aes(color = bs), alpha = 1, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                     values=c("royalblue2", "green", "darkorchid"))+
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'), 
                   labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                              '9.5','10.5','11.5'))+
  labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", 
       title = "ANNUAL SURVIVAL BY AGE CLASS")

ggsave('phi.png', phi.plot, bg='transparent', width = 15, height = 10)

