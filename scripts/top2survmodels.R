# ---- Load Packages ----
#final run of survival analyses
library(jagsUI)
library(ggplot2)
library(ggdist)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

# ---- Load data ----
# data <- read.csv('./cleaned/ch.pmdi.csv', header = T)
data <- read.csv('./cleaned/ch.pmdi.dens.csv', header = T)

# data <- data %>% filter(!birth_year < 2011 & !year < 2011) #filter for years where we have all three sites
data <- data %>%  mutate(bs = recode(bs, ey = "control")) #rename ey to control to serve as a reference class 
data <- data %>%  mutate(bs = recode(bs, dmp = "wy")) #rename dmp to WY to merge wy and dmp into one bs 

#how many capture histories do we have
count <- sum(data$status == '1')
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
bs <- as.numeric(factor(id.bs.by$bs)) # 1 = control, 2  = dmp + wy

sum(bs==1) #control
sum(bs==2) #tgt + dmp
unique(bs)

#create ageclass matrix treating age as categorical
data$age.sc <- scale(data$ageclass)
ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
ageclass<- ageclass[,-1]
ageclass<-as.matrix(ageclass)

age.sc <- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc<- as.matrix(age.sc[,-1])

#create birth year vector
birthyear <- as.numeric(as.factor(id.bs.by$birth_year))

#create capture year vector
capyear <- f

# 
# create vector with last occasion for each individual, marked by 2, 15 for end of study
# rework h to only include capture myopathy or harvest, do not censor natural mortality
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 15) #change to equal number of columns/years
h
f-h #check for zero

#age simulation for continuous model
nvalues <- 15
age.sim <- seq(from = min(age.sc, na.rm = T), to = max(age.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 


# 
#---- Model3: phi ~ int+ age + age2 + site  ----


# Specify model in JAGS language
set.seed(100)
sink("model3.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)

  beta3[1] <- 0 #ey
  beta3[2] ~ dnorm(0, 0.01)  #wy



  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]   #age
                                  + (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
                                  + beta3[ bs[ i ]] #birthsite
                                 
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = c(NA, rnorm(1,0,1)),#birth site
  
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'beta3','eps1')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model3<- jagsUI(jags.data, inits, parameters, "model3.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model3)

#---- Model4: int + age + age2 ----

# Specify model in JAGS language
set.seed(100)
sink("model4.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)

  eps1[1] <- 0 #capture year RE
  for (u in 2:14){  #prior for year effect
   eps1[u] ~ dnorm(0,tau.year)
  }

  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year  ~ dunif(0,100)

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
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]   #age
                                  + (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
                                  + eps1[year[ i ]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
          for (j in 1:15) { #simulated age
            phi.age[j] <- exp( int + beta1*age.sim[j]
                                      + beta2*age.sim[j]*age.sim[j]
                                      ) /
                               (1 + exp( int+ beta1*age.sim[j]
                                              + beta2*age.sim[j]*age.sim[j]
                                             ) )

          }


            for ( j in 2:15 ) {
              phi.age.diff[j] <- phi.age[j] - phi.age[j-1]

            }
            
             for ( j in 3:15 ) {
              phi.age.2diff[j] <- phi.age[j] - phi.age[j-2]

            }
            for ( j in 4:15 ) {
              phi.age.3diff[j] <- phi.age[j] - phi.age[j-3]

            }
            
            for ( j in 5:15 ) {
              phi.age.4diff[j] <- phi.age[j] - phi.age[j-4]

            }
            
            for ( j in 6:15 ) {
              phi.age.5diff[j] <- phi.age[j] - phi.age[j-5]

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'phi.age', 'phi.age.diff', 'phi.age.2diff',
                'phi.age.3diff', 'phi.age.4diff', 'phi.age.5diff')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model4<- jagsUI(jags.data, inits, parameters, "model4.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model4)


# #---- Model 4 Plots ----

#create a tibble of the posterior draws
gather <- model4 %>%
  spread_draws(phi.age[age])

age_lookup <- tibble(
  age = 1:length(age.sim), #change to account for age.sim
  ageclass = (age.sim * sd(ageclass, na.rm = TRUE)) +
    mean(ageclass, na.rm = TRUE)
)

gather <- gather %>%
  left_join(age_lookup, by = "age")


phi.plot<- gather %>%
  ggplot(aes(x=ageclass, y=phi.age)) +
  stat_lineribbon(.width = 0.95)+
  guides(fill = "none")+ #remove legend from ribbon
  scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
  labs(x = "Age", y = "Annual Survival Probability", title = "")+
 # scale_x_continuous(breaks = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
        axis.text = element_text(face='bold',size = 16),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/phi.model4age.jpg', phi.plot, width = 8, height = 8)


# #senescence rate plot
# gather <- model7 %>%
#   spread_draws(phi.age.diff[site,age])
# 
# age_lookup <- tibble(
#   age = 1:15, #change to account for age.sim
#   ageclass = (age.sim * sd(ageclass, na.rm = TRUE)) + 
#     mean(ageclass, na.rm = TRUE)
# )
# 
# gather <- gather %>%
#   left_join(age_lookup, by = "age")
# 
# gather$site <- as.factor(gather$site)
# gather$site <- factor(gather$site,
#                       levels = c(1, 2),
#                       labels = c("East Yana", "West Yana"))
# gather$age <- as.factor(gather$age)
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=phi.age.diff, color = site)) +
#   stat_pointinterval(position = position_dodge(.5))+
#   guides(fill = "none")+ #remove legend from ribbon
#   scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
#   labs(x = "Age", y = "Change in Survival", title = "")+
#   #scale_x_continuous(breaks = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5)) +  
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.8),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 16),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.age.diff.jpg', phi.plot, width = 8, height = 8)
