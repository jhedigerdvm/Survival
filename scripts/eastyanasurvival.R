#East Yana fawn survival 

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
data <- read.csv('./cleaned/fawncaphx.csv', header = T)

data <- data %>% filter(bs %in% "ey")

#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'cap_year', values_from = 'status_cam', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #83 individuals with known fates 
ch[indices] <- 1


# Create vector with the occasion each indiv is marked
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

#create ageclass matrix treating age as categorical
# ageclass<- pivot_wider(data, names_from = 'year', values_from = 'ageclass', id_cols = 'animal_id' )
# ageclass<- ageclass[,-1]
# ageclass<-as.matrix(ageclass)

#scale and center age
age.mean <- mean(data$age, na.rm = TRUE)
age.sd <- sd(data$age, na.rm = TRUE)
data$age.sc <- (data$age - age.mean) / age.sd

age.sc <- pivot_wider(data, names_from = 'cap_year', values_from = 'age.sc', id_cols = 'animal_id' )
age.sc<- as.matrix(age.sc[,-1])

age.pred <- seq(0.5, 15.5, by = 1)

age.pred.sc <- (age.pred - age.mean) / age.sd


#create birth year vector
birthyear <- as.numeric(as.factor(data$birth_year))

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

#identify 0 and +1 indiv
problem <- which((h - f) <= 1)

length(problem)

data.frame(
  id = data$animal_id[problem],
  f = f[problem],
  h = h[problem]
)

known.fate[problem, ]

#identify the individuals to keep and adjust lists and vectors
keep <- which((h - f) > 0)

known.fate <- known.fate[keep, ]
ch <- ch[keep, ]
f <- f[keep]
h <- h[keep]
birthyear <- birthyear[keep]
capyear <- capyear[keep]
age.sc <- age.sc[keep, ]

#verify the fix
which((h - f) <= 0)

#age simulation for continuous model
age.sim <- age.pred.sc
nvalues <- length(age.sim)


#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}

# 
#---- Model1: phi ~ int + age + age2   ----


# Specify model in JAGS language
set.seed(100)
sink("model1.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001) #age
  beta2 ~ dnorm(0,0.001) #age ^2
  
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


}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc,
                   year = capyear, age.sim = age.sim)


# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1)

  
  
  # 
  # eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2', 'phi.age')

# MCMC settings
ni <- 2000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
model1<- jagsUI(jags.data, inits, parameters, "model1.jags", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(model1)



# #---- Model 1 Plots ----

#create a tibble of the posterior draws
gather <- model1 %>%
  spread_draws(phi.age[age])
# 
# age_lookup <- tibble(
#   age = 1:length(age.sim), #change to account for age.sim
#   ageclass = (age.sim * sd(ageclass, na.rm = TRUE)) +
#     mean(ageclass, na.rm = TRUE)
# )
# 
# gather <- gather %>%
#   left_join(age_lookup, by = "age")


phi.plot<- gather %>%
  ggplot(aes(x=age, y=phi.age)) +
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

#---- Model2: piece-wise linear regression analysis    ----
#create knot/inflexion point
# knot <- 7:13
#store models for different knot points
models <- list()

dic.df <- data.frame(
  knot = 7:13,
  DIC = NA
)

#create age data
  age.raw <- pivot_wider(
    data,
    names_from = cap_year,
    values_from = age,
    id_cols = animal_id
  )
  
  age.raw <- as.matrix(age.raw[,-1])

# keep same individuals as ch
  age.raw <- age.raw[keep, ]
  
  # age.old <- pmax(age.raw - knot, 0)
  



#scale and center age
  age.mean <- mean(data$age, na.rm = TRUE)
  age.sd <- sd(data$age, na.rm = TRUE)
  data$age.sc <- (data$age - age.mean) / age.sd
  
  age.sc <- pivot_wider(data, names_from = 'cap_year', values_from = 'age.sc', id_cols = 'animal_id' )
  age.sc<- as.matrix(age.sc[,-1])

  age.pred <- seq(0.5, 15.5, by = 1)
  
  age.pred.sc <- (age.pred - age.mean) / age.sd
  
  age.old.pred <- pmax(age.pred - knot, 0)
  age.old.sim <- age.old.pred
  
  
#age simulation for continuous model
  age.sim <- age.pred.sc
  nvalues <- length(age.sim)
  
# Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
    beta2=rnorm(1,0,1)
     )
  }

#Parameters to monitor
  parameters <- c('int', 'beta1', 'beta2', 'phi.age')
  
# MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3

# Specify model in JAGS language
set.seed(100)
sink("model2.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001) #age
  beta2 ~ dnorm(0,0.001) #age.old
  
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
                                  + beta2*age.old[ i , t-1 ]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
    for (j in 1:15) { #simulated age
            phi.age[j] <- exp( int + beta1*age.sim[j]
                                      + beta2*age.old.sim[j]

                                      ) /
                               (1 + exp( int + beta1*age.sim[j]
                                                  + beta2*age.old.sim[j]


                                             ) )

          }
   
}
",fill = TRUE)
sink()

#start loop for knot points
for(k in 7:12){

  cat("\nFitting knot =", k, "\n")

  age.old <- pmax(age.raw - k, 0)
  age.old.sim <- pmax(age.pred - k, 0)
  
  #bundle data
    jags.data <- list(
      h = h,
      ch = ch,
      f = f,
      nind = nrow(ch),
      ageclass = age.sc,
      age.old = age.old,
      age.sim = age.sim,
      age.old.sim = age.old.sim
    )
    
  # Call JAGS from R (BRT 3 min)
    mod <- jagsUI(
      jags.data,
      inits,
      parameters,
      "model2.jags",
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
  
  #paste models into list
    models[[paste0("knot_",k)]] <- mod

    dic.df$DIC[dic.df$knot == k] <- mod$DIC

}

#review results from competing models to identify senescence
dic.df$deltaDIC <- dic.df$DIC - min(dic.df$DIC, na.rm = T)

dic.df



# #---- Model 2 Plots ----

#create a tibble of the posterior draws
gather <- models$knot_9 %>%
  spread_draws(phi.age[age])
# 
# age_lookup <- tibble(
#   age = 1:length(age.sim), #change to account for age.sim
#   ageclass = (age.sim * sd(ageclass, na.rm = TRUE)) +
#     mean(ageclass, na.rm = TRUE)
# )
# 
# gather <- gather %>%
#   left_join(age_lookup, by = "age")


phi.plot<- gather %>%
  ggplot(aes(x=age, y=phi.age)) +
  stat_lineribbon(.width = 0.95)+
  guides(fill = "none")+ #remove legend from ribbon
  scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
  labs(x = "Age", y = "Annual Survival Probability", title = "")+
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13)) +
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
