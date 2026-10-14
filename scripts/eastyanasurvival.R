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

data <- data %>% filter(bs %in% "ey" & !birth_year == "2022")

#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'cap_year', values_from = 'status_cam', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch)

known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #32 individuals with known fates 
ch[indices] <- 1


# Create vector with the occasion each indiv is marked
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 

#scale and center age
  age.mean <- mean(data$age, na.rm = TRUE)
  age.sd <- sd(data$age, na.rm = TRUE)
  data$age.sc <- (data$age - age.mean) / age.sd
  
  age.sc <- pivot_wider(
    data, 
    names_from = 'cap_year', 
    values_from = 'age.sc', 
    id_cols = 'animal_id' 
    )
  age.sc<- as.matrix(age.sc[,-1])
  
  age.raw <- pivot_wider(
    data,
    names_from = cap_year,
    values_from = age,
    id_cols = animal_id
  )
  age.raw <- as.matrix(age.raw[,-1])
  
  age.pred <- seq(0.5, 15.5, by = 1)
  age.pred.sc <- (age.pred - age.mean) / age.sd


#create age.old list
  age.old.list <- lapply(
    7:13,
    function(k) pmax(age.raw - k, 0)
  )

  names(age.old.list) <- paste0("knot_", 7:13)

#age simulation
  age.sim <- age.pred.sc
  nvalues <- length(age.sim)

#create birth year vector
birthyear <- as.numeric(as.factor(data$birth_year))

#create capture year vector
capyear <- f


#create pmdi data
  pmdi.vars <- c(     #create list of pmdi variables
    "pmdi_annual",
    "pmdi_spring",
    "pmdi_summer",
    "pmdi_fall",
    "pmdi_winter"
  )
  
  pmdi.stats <- data.frame(   #create dataframe to store means and SD of pmdi variables
    variable = pmdi.vars,
    mean = NA,
    sd = NA
  )
  
  for(i in seq_along(pmdi.vars)){ #for loop to scale and center pmdi values and create new SC columns
    
    v <- pmdi.vars[i]
    
    mu <- mean(data[[v]], na.rm = TRUE)
    sig <- sd(data[[v]], na.rm = TRUE)
    
    pmdi.stats$mean[i] <- mu
    pmdi.stats$sd[i] <- sig
    
    data[[paste0(v, "_sc")]] <-
      (data[[v]] - mu) / sig #scaled and centered calculation
  }

#now create list for pmdi matrices that are scaled and centered
  pmdi.vars <- c(
    "pmdi_annual_sc",
    "pmdi_spring_sc",
    "pmdi_summer_sc",
    "pmdi_fall_sc",
    "pmdi_winter_sc"
  )
  
  pmdi.list <- list()
  
  for(v in pmdi.vars){
    
    tmp <- pivot_wider(
      data,
      names_from = cap_year,
      values_from = all_of(v),
      id_cols = animal_id
    )
    
    tmp <- as.matrix(tmp[,-1])
    
    pmdi.list[[v]] <- tmp
    
  }


# create vector with last occasion for each individual, marked by 2, 15 for end of study
# rework h to only include capture myopathy or harvest, do not censor natural mortality
get.last<- function(x) min(which(x>1))
h <- apply(known.fate,1,get.last)
h <- replace(h, is.infinite(h), 15) #change to equal number of columns/years
h

#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}

# # 
# #---- Model1: phi ~ int + age + age2   ----
# 
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("model1.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#   int ~ dnorm(0, 0.001)
# 
#   beta1 ~ dnorm(0, 0.001) #age
#   beta2 ~ dnorm(0,0.001) #age ^2
#   
#   tau <- 1/(sigma*sigma)
#   sigma ~ dunif(0,100)
# 
# 
# # Likelihood
# for (i in 1:nind){
#    # Define latent state at first capture, we know for sure the animal is alive
#       z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):h[i]){
#         # State process
#             z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
#             logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]   #age
#                                   + (beta2*ageclass[ i , t-1 ]*ageclass[ i , t-1 ] )
# 
# 
#           # Observation process
#             ch[i,t] ~ dbern(mu2[i,t])
#             mu2[i,t] <- p * z[i,t]
# 
# 
# 
#       } #t
#    } #i
# 
#    #derived parameters
#    
#        # for (j in 1:15) { #simulated age
#        #      phi.age[j] <- exp( int + beta1*age.sim[j]
#        #                                + beta2*age.sim[j]*age.sim[j]
#        # 
#        #                                ) /
#        #                         (1 + exp( int+ beta1*age.sim[j]
#        #                                        + beta2*age.sim[j]*age.sim[j]
#        # 
#        # 
#        #                                       ) )
#        # 
#        #    }
# 
# 
# }
# ",fill = TRUE)
# sink()
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc,
#                    year = capyear, age.sim = age.sim)
# 
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1),
#   z = z.init,
#   beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
#   beta2=rnorm(1,0,1)
# 
#   
#   
#   # 
#   # eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
# )
# }
# 
# 
# parameters <- c('int', 'beta1', 'beta2', 'phi.age')
# 
# # MCMC settings
# ni <- 2000
# nt <- 10
# nb <- 1000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# model1<- jagsUI(jags.data, inits, parameters, "model1.jags", n.chains = nc,
#                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(model1)
# 
# 
# 
# # #---- Model 1 Plots ----
# 
# #create a tibble of the posterior draws
# gather <- model1 %>%
#   spread_draws(phi.age[age])
# # 
# # age_lookup <- tibble(
# #   age = 1:length(age.sim), #change to account for age.sim
# #   ageclass = (age.sim * sd(ageclass, na.rm = TRUE)) +
# #     mean(ageclass, na.rm = TRUE)
# # )
# # 
# # gather <- gather %>%
# #   left_join(age_lookup, by = "age")
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=phi.age)) +
#   stat_lineribbon(.width = 0.95)+
#   guides(fill = "none")+ #remove legend from ribbon
#   scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
#   labs(x = "Age", y = "Annual Survival Probability", title = "")+
#   # scale_x_continuous(breaks = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5)) +
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

#---- Model2: piece-wise linear regression analysis    ----

#store models for different knot points
models <- list()

dic.df <- data.frame(
  knot = 7:13,
  DIC = NA
)
# 
# #create beta spline with single knot point
#   age.old.pred <- pmax(age.pred - knot, 0)
#   age.old.sim <- age.old.pred
  
#age simulation for continuous model
  age.sim <- age.pred.sc
  nvalues <- length(age.sim)
  
# Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), 
    beta2=rnorm(1,0,1)
     )
  }

#Parameters to monitor
  parameters <- c('int', 
                  'beta1', 
                  'beta2', 
                  'phi.age')
  
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
            z[i,t] ~ dbern(mu1[i,t]) 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]  
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
for(k in 7:13){

  cat("\nFitting knot =", k, "\n")

  age.old <- age.old.list[[paste0("knot_", k)]]
  age.old.sim <- pmax(age.pred - k, 0)
  
  # age.old <- pmax(age.raw - k, 0)
  # age.old.sim <- pmax(age.pred - k, 0)
  
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


#---- Model3: piece-wise linear regression analysis, PMDI    ----

#keep simulated age.old to 9 
age.old <- age.old.list$knot_9
age.old.sim <- pmax(age.pred - 9, 0) #knot point of 9 from previous model selection

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  beta3 = rnorm(1,0,1)
)
}

#Parameters to monitor
parameters <- c('int', 'beta1', 'beta2', 'beta3', 'phi.age')

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
  beta3 ~ dnorm(0,0.001)
  
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
                                  + beta3*drought[i,t-1]


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


#create dataframe to hold pmdi data
  results.df <- data.frame(
    pmdi = pmdi.vars,
    DIC = NA
  )
  
  models <- list()
  
#run for loop with all pmdi variables
for(i in seq_along(pmdi.vars)){
  
  v <- pmdi.vars[i]
  
  cat("\nFitting", v, "\n")
  
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    ageclass = age.sc,
    age.old = age.old,   # age.old for knot = 9
    drought = pmdi.list[[v]],
    age.sim = age.sim,
    age.old.sim = age.old.sim
  )
  
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
  
  models[[v]] <- mod
  
  results.df$DIC[i] <- mod$DIC
}

#review results from competing models to identify pmdi
  results.df$deltaDIC <- results.df$DIC - min(results.df$DIC, na.rm = T)
  
  results.df
  
  
  
  
  # #---- Model 3 Plots ----
  print(models$pmdi_summer_sc)
  #create a tibble of the posterior draws
  gather <- models$pmdi_summer_sc %>%
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
    geom_vline(xintercept = 9, linetype = "dashed", color = "red", linewidth = 1)+
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
  
  

#---- Model3: simulation    ----

  set.seed(123)
  
  # TRUE PARAMETERS

  int.true   <- 2
  beta1.true <- 0.2
  beta2.true <- -0.8
  
  p.true <- 0.85
  
  knot <- 9
  
  nind  <- 500
  nyear <- 15
  
  
  # AGE MATRICES

  age.raw <- matrix(
    rep(seq(0.5, 14.5, by = 1),
        each = nind),
    nrow = nind,
    byrow = FALSE
  )
  
  age.old <- pmax(age.raw - knot, 0)
  
  
  # TRUE SURVIVAL

  phi.true <- plogis(
    int.true +
      beta1.true * age.raw +
      beta2.true * age.old
  )
  
  
  # SIMULATE LATENT STATES

  z.true <- matrix(NA, nind, nyear)
  
  z.true[,1] <- 1
  
  for(i in 1:nind){
    
    for(t in 2:nyear){
      
      z.true[i,t] <- rbinom(
        1,
        1,
        phi.true[i,t-1] * z.true[i,t-1]
      )
      
    }
    
  }
  z.init <- z.true
  z.init[,1] <- NA
  

  # SIMULATE ENCOUNTER HISTORIES

  ch <- matrix(0, nind, nyear)
  
  ch[,1] <- 1
  
  for(i in 1:nind){
    
    for(t in 2:nyear){
      
      ch[i,t] <- rbinom(
        1,
        1,
        p.true * z.true[i,t]
      )
      
    }
    
  }
  
  
  # FIRST CAPTURE OCCASION

  f <- rep(1, nind)
  
  # run through last occasion
  
  h <- rep(nyear, nind)
  
  
  # JAGS DATA

  jags.data <- list(
    ch = ch,
    f = f,
    h = h,
    nind = nind,
    ageclass = age.raw,
    age.old = age.old
  )
  
  
  # INITS

  inits <- function(){
    
    list(
      int   = rnorm(1,0,1),
      beta1 = rnorm(1,0,1),
      beta2 = rnorm(1,0,1),
      z = z.init
    )
    
  }
  
  
  # PARAMETERS

  parameters <- c(
    "p",
    "int",
    "beta1",
    "beta2"
  )
  
  
  # WRITE MODEL

  writeLines(
    "
model {

  p ~ dbeta(1,1)

  int   ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)

  for(i in 1:nind){

    z[i,f[i]] <- 1

    for(t in (f[i]+1):h[i]){

      z[i,t] ~ dbern(mu1[i,t])

      mu1[i,t] <- phi[i,t-1] * z[i,t-1]

      logit(phi[i,t-1]) <-
           int
         + beta1*ageclass[i,t-1]
         + beta2*age.old[i,t-1]

      ch[i,t] ~ dbern(mu2[i,t])

      mu2[i,t] <- p*z[i,t]

    }

  }

}
",
"sim_model.jags"
  )
  
  
  # FIT MODEL

  mod <- jags(
    data = jags.data,
    inits = inits,
    parameters.to.save = parameters,
    model.file = "sim_model.jags",
    n.chains = 3,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 10,
    parallel = TRUE
  )
  
  print(mod)
  
  
  
  # ---- simulation 2 ----
  data <- read.csv('./cleaned/fawncaphx.csv', header = T)
  
  data <- data %>% filter(bs %in% "ey" & !birth_year == "2022")
  
  #take long form and convert into wide for CH matrix
  ch<- pivot_wider(data, names_from = 'cap_year', values_from = 'status_cam', id_cols = 'animal_id' )
  ch<-ch[,-1]
  ch<-as.matrix(ch)
  
  known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest
  
  #create capture history with just 1s and 0s, remove 'known fates'
  indices <- which(ch == 2, arr.ind = TRUE) #32 individuals with known fates 
  ch[indices] <- 1
  
  
#input true data
  set.seed(123)
  
  nind  <- 120
  nyear <- 16
  p.true <-0.85 #capture probability

  int.true   <- 2
  beta1.true <- 0.2
  beta2.true <- -0.8
  knot <- 9

#build age matrix
  age.raw <- matrix(
    rep(seq(0.5, 16.5, by = 1),
        each = nind),
    nrow = nind,
    byrow = FALSE
  )
  
  age.old <- pmax(age.raw - knot, 0)
  
  
# TRUE SURVIVAL
  phi.true <- plogis(
    int.true +
      beta1.true * age.raw +
      beta2.true * age.old
  )
  
  
#calculate z.true, latent state, we know whether they are alive or dead
  z.true <- matrix(NA, nind, nyear)
  z.true[,1] <- 1 #individual alive at t
  for(i in 1:nind){
    for(t in 2:nyear){
      z.true[i,t] <- rbinom( 1, 1, phi.true[i,t-1] * z.true[i,t-1] #if z.true is 0 at any point, it stays 0
      )
    }  }
  
  # z.init <- z.true
  # z.init[,1] <- NA
  
#create simulated capture histories   
  ch.sim <- matrix(0, nind, nyear)
  ch.sim[,1] <- 1 #everyone is caught the first year
  
  for(i in 1:nind){
    for(t in 2:nyear){
      ch.sim[i,t] <- rbinom(1,1, p.true * z.true[i,t] #once z.true is zero, the animal is DEAD and cannot be recaught regardles of p
      )
      }  }
  
# Create vector with the occasion each indiv is marked
  get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
  f <- apply(ch.sim, 1, get.first) 
  
# create vector with last occasion for each individual, marked by 2, 15 for end of study
  get.last<- function(x) min(which(x>1))
  h <- apply(ch.sim,1,get.last)
  h <- replace(h, is.infinite(h), nyear) #change to equal number of columns/years
  h
  
#Function for latent state
  z.init <- matrix(NA, nrow = nrow(ch.sim), ncol = ncol(ch.sim))
  
  for(i in 1:dim(z.init)[1]){
    z.init[i, f[i]:h[i]] <- 1
    z.init[i,f[i]] <- NA
  }
  
# Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), 
    beta2=rnorm(1,0,1)
  )
  }
  
#Parameters to monitor
  parameters <- c('int', 
                  'beta1', 
                  'beta2', 
                  'p')
  
# MCMC settings
  ni <- 2000
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
            z[i,t] ~ dbern(mu1[i,t]) 
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  
            logit(phi[i,t-1]) <- int + beta1*ageclass[ i , t-1 ]  
                                  + beta2*age.old[ i , t-1 ]
                                  


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   
   
}
",fill = TRUE)
  sink()
    
#bundle data
    jags.data <- list(
      h = h,
      ch = ch.sim,
      f = f,
      nind = nrow(ch.sim),
      ageclass = age.raw,
      age.old = age.old
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
    
print(mod)
    
  