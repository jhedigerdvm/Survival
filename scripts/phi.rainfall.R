#Survival model with interaction between rainfall and birth site 
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(here)

data<- read.csv('./cleaned/caphx.rainfall.long.csv', header = T)


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


annual.rainfall<-pivot_wider(data, names_from = 'year', values_from = 'annual.sc', id_cols = 'animal_id' )
annual.rainfall<-annual.rainfall[,-1]
annual.rainfall<-as.matrix(annual.rainfall)


nvalues <- 1000
rain.sim <- seq(from = min(annual.rainfall, na.rm = T), to = max(annual.rainfall, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
rain.sim


# Specify model in JAGS language
set.seed(100)
sink("cjs-rain-site.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)
   

#priors
int ~ dnorm(0,0.01)

site.beta[1] <- 0

for (u in 2:3){
  site.beta[u] ~ dnorm(0,0.01)
}

rain.beta ~ dnorm(0,0.01)


rain.site.beta[1] <-0

for (u in 2:3){
  rain.site.beta[u] ~ dnorm(0,0.01)
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1
      
      for (t in (f[i]+1):h[i]){ 
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <- int + site.beta[bs[i]] + rain.beta*rain[i,t-1]  + rain.site.beta[bs[i]]*rain[i,t-1]
          
          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]
      } #t
   } #i

#derived parameters
  for (i in 1:3){ #site and rain.site.beta
    for (j in 1:1000){ #rain
      survival[j,i] <- exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j])/
                            (1 + exp(int + site.beta[i] + rain.beta*rain.sim[j] + rain.site.beta[i]*rain.sim[j]))
      }                     #delta method to convert from logit back to probability Powell et al. 2007
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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch),  bs = bs, rain = annual.rainfall, rain.sim = rain.sim)

# Initial values
inits <- function(){list(int = rnorm(1,0,1), z = z.init, rain.beta = rnorm(1, 0, 1),
                         site.beta = c(NA, rnorm(2,0,1)), rain.site.beta = c(NA, rnorm(2,0,1)))} #

parameters <- c('int', 'site.beta', 'rain.beta', 'rain.site.beta', 'p','survival' )

# MCMC settings
ni <- 1000
nt <- 1
nb <- 100
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.rain.site <- jagsUI(jags.data, inits, parameters, "cjs-rain-site.jags", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)

print(cjs.rain.site)
summary<- cjs.rain.site$summary
write.csv(summary, './output/rain.bs.csv', row.names = F)



#
#create a tibble of the posterior draws
posterior<- tidy_draws(cjs.rain.site)
posterior<- posterior[,-c(1:12)]
posterior <- posterior[,-3001]

#create dataframe with posteriors of just survival age1 across the three sites
#pivot longer puts them in a tibble format
posterior_long <- posterior %>% pivot_longer(everything())
#make a new column 'site', rep 1:3 assigns 1 to site 1, 2 to site 2 for the number of rows divided by 3
posterior_long$rain <- rep(rain.sim, nrow(posterior_long)/1000)
# 
# #need to unscale and uncenter rainfall data
# mean(data$annual, na.rm = T) #23.62
# min(data$annual, na.rm = T) # 12
# max(data$annual, na.rm = T) # 41.8
# sd(data$annual) #7

# To undo center and scaling:
posterior_long$rain1<- (posterior_long$rain * 7) + 23.62

# #make a new column for site
# posterior_long$site <- '1'

# Set total rows
num_rows <- 8100000

# Create empty values vector
site <- vector(length = num_rows)

# Index for filling
index <- 1

# Repeat till last row
while(index <= num_rows){

  # Insert 1
  site[index:(index+999)] <- 1
  index <- index + 1000

  # Insert 2
  site[index:(index+999)] <- 2
  index <- index + 1000

  # Insert 3
  site[index:(index+999)] <- 3
  index <- index + 1000

}
posterior_long<- cbind(posterior_long,site)

posterior_long$site<- as.factor(posterior_long$site)

##GGPLOT
plot_base_phi <-
  ggplot(data = posterior_long, aes(x=rain1, y=value, group = site))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        plot.title = element_text(face = 'bold', size = 36, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

phi.plot<- plot_base_phi +
  geom_smooth(method = lm, aes(color = site)) +
  scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                     values=c("orange", "black", "darkorchid"))+
  # scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'),
  #                  labels = c('1.5-2.5', '2.5-3.5', '3.5-4.5' ,'4.5-5.5' ,'5.5-6.5' ,'6.5-7.5','7.5-8.5',
  #                             '8.5-9.5','9.5-10.5','10.5-11.5', '11.5-12.5'))+
  labs(x = "ANNUAL RAINFALL (IN)", y = "ANNUAL SURVIVAL PROBABILITY",
       title = "RAINFALL AND ANNUAL SURVIVAL BY SITE")

# ggsave('phi_age_site.png', phi.plot, bg='transparent', width = 15, height = 10)
ggsave('./figures/phi_rain_site2.jpg', phi.plot, width = 12, height = 10)

# compute wAIC for model with covariate
library(R2jags)
samples.m1 <- jags.samples(cjs.rain.site$model, 
                           c("WAIC","deviance"), 
                           type = "mean", 
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)

samples.m1$p_waic <- samples.m1$WAIC
samples.m1$waic <- samples.m1$deviance + samples.m1$p_waic
tmp <- sapply(samples.m1, sum)
waic.m1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

# The difference in wAIC tells us that the covariate has some effect
# wAIC of m1 the model with covariate << wAIC of m0 intercept only
data.frame(m1 = waic.m1)
