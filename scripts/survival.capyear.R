#code to look annual survival across years, NO birthsite effects 

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
sink("cjs-capyear.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors

for(t in 1:n.occasions){ # 16 columns in CH, should I be modeling survival for each capyear?
  alpha[t] ~ dunif(0, 5)
}
  
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
   z[i,f[i]] <- 1

      for (t in (f[i]+1):n.occasions){ #
      # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- alpha[t-1]

        # Observation process
            CH[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameter
  for (j in 1:16){ #cap year
    survival[j] <- exp(alpha[j] )/ (1 + exp(alpha[j]))
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
 
}


# Bundle data
jags.data <- list(CH = CH, f = f, nind = nrow(CH),  n.occasions = ncol(CH))# ,h = h,bs = bs , ageclass = ageclass,capyear=capyear


# Initial values
inits <- function(){list(z = z.init, alpha = runif(16,0,5))}
#,site.beta = c(NA, rnorm(2,0,1)),int   int= rnorm(1,0,1),  age.beta = c(NA, rnorm(14,0,1)), , sigma.capyear = runif(1,0,5)

parameters <- c('alpha', 'p', 'survival' )
# 'survival_age1_site1','survival_age2_site1', 'survival_age3_site1', 'survival_age4_site1','site_diff'
# 'survival_age1_site2','survival_age2_site2', 'survival_age3_site2','survival_age4_site2',
# 'survival_age1_site3','survival_age2_site3', 'survival_age3_site3', 'survival_age4_site3') 'survival',, 'site_diff','eps.capyear', 'survival'

# MCMC settings
ni <- 3000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.capyear <- jagsUI(jags.data, inits, parameters, "cjs-capyear.jags", n.chains = nc, 
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
MCMCtrace(cjs.capyear)
print(cjs.capyear)
# summary<- cjs.age$summary
# write.csv(summary,'surv_age_output1.csv')


#create a tibble of the posterior draws
posterior<- tidy_draws(cjs.capyear)
posterior<- posterior[,c(21:36)]

#create dataframe with posteriors of just survival age1 across the three sites
#pivot longer puts them in a tibble format
posterior_long <- posterior %>% pivot_longer(everything())
#make a new column 'site', rep 1:3 assigns 1 to site 1, 2 to site 2 for the number of rows divided by 3
posterior_long$capyear <- rep(c('2007', '2008', '2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020','2021', '2022'), nrow(posterior_long)/16)

##GGPLOT
plot_base_phi <- 
  ggplot(data = posterior_long, aes(x=capyear, y=value))+
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
  stat_pointinterval(alpha = 1, .width = c(0.5, 0.95)) +
  scale_x_discrete(limits=c('2007', '2008', '2009','2010','2011','2012','2013',
                            '2014','2015','2016','2017','2018','2019','2020','2021'))+
                   
  labs(x = "YEAR", y = "ANNUAL SURVIVAL PROBABILITY", 
       title = "ANNUAL SURVIVAL BY YEAR")

ggsave('phi_capyear.png', phi.plot, bg='transparent', width = 15, height = 10)
ggsave('phi_capyear.jpg', phi.plot, width = 15, height = 10)
