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
# data <- read.csv('./cleaned/fawncaphx.csv', header = T)
# 
# data <- data %>% filter(bs %in% "ey" & !birth_year == "2022")
# write.csv(data, './cleaned/eydatawfawns.csv')

data<- read.csv('./cleaned/eydatawfawns.csv', header = T)
#remove individual 90-2015-0-560 because first and last capture occasions are the same
data <- data %>%
  filter(animal_id != "90-2015-0-560")

#remove individuals with known fates to not confound senescence data
    # data_clean <- data %>%
    #   group_by(animal_id) %>%
    #   filter(
    #     !any(status %in% c(
    #       "Harvested",
    #       "Nat Mortality",
    #       "Dead-Capture"
    #     ))
    #   ) %>%
    #   ungroup()
    # 
    # removed_ids <- data %>%
    #   group_by(animal_id) %>%
    #   filter(any(status %in% c(
    #     "Harvested",
    #     "Nat Mortality",
    #     "Dead-Capture"
    #   ))) %>%
    #   distinct(animal_id)
    # 
    # removed_ids #32 individuals removed

    # data <- data_clean
#take long form and convert into wide for CH matrix
ch<- pivot_wider(data, names_from = 'cap_year', values_from = 'status_cam', id_cols = 'animal_id' )
ch<-ch[,-1]
ch<-as.matrix(ch) #
# 
known.fate <- ch #known fate matrix with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #32 individuals with known fates
# ch[indices] <- 1


# Create vector with the occasion each indiv is marked
get.first <- function(x) min(which(x!=0)) #x! identifies when x is not equal to zero
f <- apply(ch, 1, get.first) 


# create vector with last occasion for each individual, marked by 2, 16 for end of study
get.last<- function(x) min(which(x>1))
h <- apply(ch,1,get.last) #changed from known fate
h <- replace(h, is.infinite(h), ncol(ch)) #change to equal number of columns/years
h
diff<- (f-h)  #remove any whose first capture was also their last, row 150 in ch
sum(diff == 0)

#Function for latent state
z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))

for(i in 1:dim(z.init)[1]){
  z.init[i, f[i]:h[i]] <- 1
  z.init[i,f[i]] <- NA
}
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }




#age data
    age.raw <- pivot_wider(
      data,
      names_from = cap_year,
      values_from = age,
      id_cols = animal_id
    )
    age.raw <- as.matrix(age.raw[,-1])
    
    age.sim <- seq(0.5, 15.5, by = 1) # updated from 15.5
    
    #age simulation
    nvalues <- length(age.sim)
    
    #create age.old list
    age.old.list <- lapply(
      6:12, #changed from 7:13, we are interested in identifying the onset of senescence
      function(k) pmax(age.raw - k, 0)
    )
    
    names(age.old.list) <- paste0("knot_", 6:12) #changed from 7:13
    age.old.list$knot_6
    
    #create ageclasses
    data$ageclass <- data$age + 0.5 # 1 are fawns, 2 are yearlings, 3 and 2.5 yos

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


#Add density data
data$deerkm2.sc <- scale(data$adult_density)
density.sc <- pivot_wider(
  data, 
  names_from = 'cap_year', 
  values_from = 'deerkm2.sc', 
  id_cols = 'animal_id' 
)
density.sc<- as.matrix(density.sc[,-1])



#---- Model_piecewise: piece-wise linear regression analysis    ----


# Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1),
    beta2= rnorm(1,0,1)
     )
  }

#Parameters to monitor
  parameters <- c('int',
                  'beta1',
                  'beta2')

# MCMC settings
  ni <- 8000
  nt <- 10
  nb <- 5000
  nc <- 3

# Specify model in JAGS language
set.seed(100)
sink("model.jags")
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

      for (t in (f[i]+1):16){
        # State process
            z[i,t] ~ dbern(mu1[i,t])
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]
                                  + beta2*age.old[ i , t-1 ]



          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

    

}
",fill = TRUE)
sink()



#bundle data
jags.data <- list(
  #h = h,
  ch = ch,
  f = f,
  nind = nrow(ch),
  age = age.raw,
  age.old = age.old.list$knot_6
)

# Call JAGS from R
mod6 <- jagsUI(
  jags.data,
  inits,
  parameters,
  "model.jags",
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE
)

print(mod6)
MCMCtrace(mod6)
#
# 
# 
# #bundle data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old.list$knot_7
#     )
# 
# # Call JAGS from R
#     mod7 <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
# 
#     print(mod7)
#     MCMCtrace(mod7)
# 
# 
# #bundle data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old.list$knot_8
#     )
# 
# # Call JAGS from R
#     mod8 <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
# 
#     print(mod8)
#     MCMCtrace(mod8)
# 
# 
# #bundle data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old.list$knot_9
#     )
# 
# # Call JAGS from R
#     mod9 <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
# 
#     print(mod9)
#     MCMCtrace(mod9)

# 
# #Knot 10
# #bundle data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old.list$knot_10
#     )
# 
# # Call JAGS from R
#     mod10 <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
# 
#     print(mod10)
#     MCMCtrace(mod10)
#
# #knot 11
# #bundle data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old.list$knot_11
#     )
# 
# # Call JAGS from R
#     mod11 <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
# 
#     print(mod11)
#     MCMCtrace(mod11)

  #knot 12
  #bundle data
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    age = age.raw,
    age.old = age.old.list$knot_12
  )
  
  # Call JAGS from R
  mod12 <- jagsUI(
    jags.data,
    inits,
    parameters,
    "model.jags",
    n.chains = nc,
    n.thin = nt,
    n.iter = ni,
    n.burnin = nb,
    parallel = TRUE
  )
  
  print(mod12)
  MCMCtrace(mod12) 
    
#store models for different knot points
mods <- list(
  
  mod6 = mod6,
  mod7 = mod7,
  mod8 = mod8,
  mod9 = mod9,
  mod10 = mod10,
  mod11 = mod11,
  mod12 = mod12
)

dic.df <- data.frame(
  knot = 6:12,
  pD = sapply(mods, function(x) x$pD),
  Deviance = sapply(mods, function(x) x$mean$deviance),
  DIC = sapply(mods, function(x) x$DIC)
)

dic.df$deltaDIC <- dic.df$DIC - min(dic.df$DIC)

dic.df <- dic.df[order(dic.df$DIC), ]

dic.df
write.csv(dic.df, './output/knotpointDIC.csv')

# # #---- Piece-wise regression for top age model ----

#keep simulated age.old to 6 
age.old <- age.old.list$knot_6
age.old.sim <- pmax(age.sim - 6, 0) #knot point of 9 from previous model selection


# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1),
  beta2= rnorm(1,0,1)
)
}

#Parameters to monitor
parameters <- c('int',
                'beta1',
                'beta2',
                'phi.age')

# MCMC settings
ni <- 8000
nt <- 10
nb <- 5000
nc <- 3

# Specify model in JAGS language
set.seed(100)
sink("model.jags")
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]
                                  + beta2*age.old[ i , t-1 ]



          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
      for (j in 1:14) { #simulated age
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


#bundle data
jags.data <- list(
  h = h,
  ch = ch,
  f = f,
  nind = nrow(ch),
  age = age.raw,
  age.old = age.old.list$knot_6,
  age.sim = age.sim,
  age.old.sim = age.old.sim
)

# Call JAGS from R 
mod6 <- jagsUI(
  jags.data,
  inits,
  parameters,
  "model.jags",
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE
)

print(mod6)
MCMCtrace(mod6)


# # #---- Piece-wise regression  age model with random effect of year ----
# 
# #keep simulated age.old to 6 
# age.old <- age.old.list$knot_6
# age.old.sim <- pmax(age.sim - 6, 0) #knot point of 9 from previous model selection
# 
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1),
#   z = z.init,
#   beta1 = rnorm(1,0,1),
#   beta2= rnorm(1,0,1)#,
#   #eps.year = c(NA, rnorm(14,0,1))
# )
# }
# 
# #Parameters to monitor
# parameters <- c('int',
#                 'beta1',
#                 'beta2',
#                 'phi.age',
#                 'eps.year')
# 
# # MCMC settings
# ni <- 8000
# nt <- 10
# nb <- 5000
# nc <- 3
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("model.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#   int ~ dnorm(0, 0.001)
#   eps.year[1] <- 0
#   
#   beta1 ~ dnorm(0, 0.001) #age
#   beta2 ~ dnorm(0,0.001) #age.old
# 
#   tau <- 1/(sigma*sigma)
#   sigma ~ dunif(0,100)
#   
#   for (u in 2:15){
#     eps.year[u] ~ dnorm(0,eps.tau)
#   }
# 
#   eps.tau <- 1/(eps.sigma*eps.sigma)
#   eps.sigma ~ dunif(0,100)
#   
# 
# 
# # Likelihood
# for (i in 1:nind){
#    # Define latent state at first capture, we know for sure the animal is alive
#       z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):h[i]){
#         # State process
#             z[i,t] ~ dbern(mu1[i,t])
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]
#             logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]
#                                   + beta2*age.old[ i , t-1 ]
#                                   + eps.year[capyear[i]]
# 
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
#       for (j in 1:14) { #simulated age
#             phi.age[j] <- exp( int + beta1*age.sim[j]
#                                       + beta2*age.old.sim[j]
#                                       ) /
#                                (1 + exp( int + beta1*age.sim[j]
#                                                   + beta2*age.old.sim[j]
#                                              ) )
# 
#     }
#     
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #bundle data
# jags.data <- list(
#   h = h,
#   ch = ch,
#   f = f,
#   nind = nrow(ch),
#   age = age.raw,
#   age.old = age.old.list$knot_6,
#   age.sim = age.sim,
#   age.old.sim = age.old.sim,
#   capyear = capyear
# )
# 
# # Call JAGS from R 
# mod6capyear <- jagsUI(
#   jags.data,
#   inits,
#   parameters,
#   "model.jags",
#   n.chains = nc,
#   n.thin = nt,
#   n.iter = ni,
#   n.burnin = nb,
#   parallel = TRUE
# )
# 
# print(mod6capyear)
# MCMCtrace(mod6)
# 
# 
#
#create a tibble of the posterior draws
gather <- mod6 %>%
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
  geom_vline(xintercept = 6, linetype = "dashed", color = "red", linewidth = 1)+
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


#---- Model_pmdi: piece-wise linear regression analysis with a single knot point, PMDI    ----

#keep simulated age.old to 6 
age.old <- age.old.list$knot_6
# age.old.sim <- pmax(age.sim - 6, 0) #knot point of 9 from previous model selection

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), #age
  beta2=rnorm(1,0,1),   #age.old
  beta3 = rnorm(1,0,1)#, #pmdi
  #eps.year = c(NA, rnorm(15, 0, 1))
)
}

#Parameters to monitor
parameters <- c('int', 'beta1', 'beta2', 'beta3')

# MCMC settings
ni <- 8000
nt <- 10
nb <- 5000
nc <- 3

# Specify model in JAGS language
set.seed(100)
sink("model_pmdi.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)
  eps.year[1] <- 0
  
  beta1 ~ dnorm(0, 0.01) #age
  beta2 ~ dnorm(0,0.01) #age.old
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                  + beta2*age.old[ i , t-1 ]
                                  + beta3*drought[i,t-1]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
}
",fill = TRUE)
sink()
# 
# #bundle annual pmdi data
#   jags.data <- list(
#     h = h,
#     ch = ch,
#     f = f,
#     nind = nrow(ch),
#     age = age.raw,
#     age.old = age.old,   # age.old for knot = 6
#     drought = pmdi.list$pmdi_annual_sc#,
#     # age.sim = age.sim,
#     # age.old.sim = age.old.sim
#     # capyear = capyear
#   )
#   
#   mod_annual <- jagsUI(
#     jags.data,
#     inits,
#     parameters,
#     "model_pmdi.jags",
#     n.chains = nc,
#     n.thin = nt,
#     n.iter = ni,
#     n.burnin = nb,
#     parallel = TRUE
#   )
#   
#   print(mod_annual)
# 
# 
# #bundle spring pmdi data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old,   # age.old for knot = 9
#       drought = pmdi.list$pmdi_spring_sc
#     )
#     
#     mod_spring <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model_pmdi.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
#     print(mod_spring)
# 
# #bundle summer data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old,   
#       drought = pmdi.list$pmdi_summer_sc
#     )
#     
#     mod_summer <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model_pmdi.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = ni,
#       n.burnin = nb,
#       parallel = TRUE
#     )
#     print(mod_summer)
# 
# #bundle fall data
#     jags.data <- list(
#       h = h,
#       ch = ch,
#       f = f,
#       nind = nrow(ch),
#       age = age.raw,
#       age.old = age.old,   
#       drought = pmdi.list$pmdi_fall_sc
#     )
#     
#     mod_fall <- jagsUI(
#       jags.data,
#       inits,
#       parameters,
#       "model_pmdi.jags",
#       n.chains = nc,
#       n.thin = nt,
#       n.iter = 8000,
#       n.burnin = 5000,
#       parallel = TRUE
#     )
#     print(mod_fall)
#     MCMCtrace(mod_fall)

#bundle winter data
    jags.data <- list(
      h = h,
      ch = ch,
      f = f,
      nind = nrow(ch),
      age = age.raw,
      age.old = age.old,   # age.old for knot = 9
      drought = pmdi.list$pmdi_winter_sc
    )
    
    mod_winter <- jagsUI(
      jags.data,
      inits,
      parameters,
      "model_pmdi.jags",
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
    print(mod_winter)
    MCMCtrace(mod_winter)

#compile DIC results from pmdi models
results.dic.pmdi <- data.frame(
  pmdi = c(
    "annual",
    "spring",
    "summer",
    "fall",
    "winter"
  ),
  pD = c(
    mod_annual$pD,
    mod_spring$pD,
    mod_summer$pD,
    mod_fall$pD,
    mod_winter$pD
  ),
  Deviance = c(
    mod_annual$mean$deviance,
    mod_spring$mean$deviance,
    mod_summer$mean$deviance,
    mod_fall$mean$deviance,
    mod_winter$mean$deviance
  ),
  DIC = c(
    mod_annual$DIC,
    mod_spring$DIC,
    mod_summer$DIC,
    mod_fall$DIC,
    mod_winter$DIC
  )
)

results.dic.pmdi$deltaDIC <- results.dic.pmdi$DIC - min(results.dic.pmdi$DIC)

results.dic.pmdi <- results.dic.pmdi[order(results.dic.pmdi$DIC), ]

results.dic.pmdi

write.csv(results.dic.pmdi, './output/pmdi_dic.csv') #annual and spring are competing models


# #---- Model PMDI Plots ----
# 
# #create a tibble of the posterior draws
# gather <- mod_summer %>%
#   spread_draws(phi.age[age])
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=phi.age)) +
#   stat_lineribbon(.width = 0.95)+
#   geom_vline(xintercept = 9, linetype = "dashed", color = "red", linewidth = 1)+
#   guides(fill = "none")+ #remove legend from ribbon
#   scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
#   labs(x = "Age", y = "Annual Survival Probability", title = "")+
#   scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13)) +
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
# 
# 
# 
# 

#---- Full model: age + drought + density     ----

# Specify model in JAGS language
set.seed(100)
sink("model_full.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.01) #age
  beta2 ~ dnorm(0,0.001) #age.old
  beta3 ~ dnorm(0,0.001)
  beta4 ~ dnorm(0, 0.001)
  
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                 + beta2*age.old[ i , t-1 ]
                                  + beta3*drought[i,t-1]
                                   + beta4*density[i,t-1]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
sink()

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), 
  beta2=rnorm(1,0,1),
  beta3 = rnorm(1,0,1),
  beta4=rnorm(1,0,1)
)
}

#Parameters to monitor
parameters <- c('int', 'beta1', 'beta2', 'beta3', 'beta4')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3


#full model
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    age = age.raw,
    age.old = age.old.list$knot_6,
    drought = pmdi.list$pmdi_summer_sc,
    density = density.sc
  )
  
  model_full <- jagsUI( #full model
    jags.data,
    inits,
    parameters,
    "model_full.jags",
    n.chains = nc,
    n.thin = nt,
    n.iter = ni,
    n.burnin = nb,
    parallel = TRUE
  )
  
print(model_full)
MCMCtrace(model_full)
  
#---- Model: age + drought     ----
  
# Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.01) #age
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                 + beta2*age.old[ i , t-1 ]
                                  + beta3*drought[i,t-1]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), 
    beta2=rnorm(1,0,1),
    beta3 = rnorm(1,0,1)
  )
  }
  
  #Parameters to monitor
  parameters <- c('int', 'beta1', 'beta2', 'beta3')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
#age + drought  
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    age = age.raw,
    age.old = age.old,   # age.old for knot = 9
    drought = pmdi.list$pmdi_summer_sc#,
    # age.sim = age.sim,
    # age.old.sim = age.old.sim,
    # density = density.sc
  )
  
  model_agedrought <- jagsUI( #full model
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
  
print(model_agedrought)
  
#---- model: age + density     ----
  
  # Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.01) #age
  beta2 ~ dnorm(0,0.001) #age.old
  beta4 ~ dnorm(0, 0.001)
  
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                 + beta2*age.old[ i , t-1 ]
                                   + beta4*density[i,t-1]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), 
    beta2=rnorm(1,0,1),
    beta4=rnorm(1,0,1)
  )
  }
  
  #Parameters to monitor
  parameters <- c('int', 'beta1', 'beta2', 'beta4')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
#age + density  
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    age = age.raw,
    age.old = age.old,   # age.old for knot = 9
    # drought = pmdi.list$pmdi_spring_sc,
    density = density.sc
  )
  
  model_agedensity <- jagsUI( #full model
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

  print(model_agedensity)
  
  #---- model: age     ----
  
  # Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.01) #age
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                 + beta2*age.old[ i , t-1 ]
                                  


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta1 = rnorm(1,0,1), 
    beta2=rnorm(1,0,1)
    
  )
  }
  
  #Parameters to monitor
  parameters <- c('int', 'beta1', 'beta2')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
  
#age   
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    age = age.raw,
    age.old = age.old#,   # age.old for knot = 9
    
  )
  
  model_age <- jagsUI( #full model
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

print(model_age)  
  
  
#----  model: drought      ----
  
  # Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  
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
            logit(phi[i,t-1]) <- int 
                                 
                                  + beta3*drought[i,t-1]
                                   


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    beta3 = rnorm(1,0,1)
  )
  }
  
  #Parameters to monitor
  parameters <- c('int', 'beta3')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
#drought   
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    
    drought = pmdi.list$pmdi_summer_sc
  )
  
  model_drought <- jagsUI( #full model
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
  
  print(model_drought)
  
  #----  model: density     ----
  
  # Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  
  beta4 ~ dnorm(0, 0.001)
  
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
            logit(phi[i,t-1]) <- int 
                                   + beta4*density[i,t-1]


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init,
    
    beta4=rnorm(1,0,1)
  )
  }
  
  #Parameters to monitor
  parameters <- c('int', 'beta4')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
#density   
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch),
    
    density = density.sc
  )
  
  model_density <- jagsUI( #full model
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
  
  print(model_density)
  
  #----  model: null     ----
  
  # Specify model in JAGS language
  set.seed(100)
  sink("model2.jags")
  cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  
  
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
            logit(phi[i,t-1]) <- int 


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
   
}
",fill = TRUE)
  sink()
  
  # Initial values
  inits <- function(){list(
    int = rnorm(1,0,1),
    z = z.init
    
  )
  }
  
  #Parameters to monitor
  parameters <- c('int')
  
  # MCMC settings
  ni <- 5000
  nt <- 10
  nb <- 1000
  nc <- 3
  
#null   
  jags.data <- list(
    h = h,
    ch = ch,
    f = f,
    nind = nrow(ch)
  )
  
  model_null <- jagsUI( #full model
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
  
  print(model_null)

# #---- Model Selection Tables ----
  
  models.df <- data.frame(
    model = c(
      "full",
      "age + drought",
      "age + density",
      "age",
      "drought",
      "density",
      "null"
    ),
    pD = c(
      model_full$pD,
      model_agedrought$pD,
      model_agedensity$pD,
      model_age$pD,
      model_drought$pD,
      model_density$pD,
      model_null$pD
    ),
    Deviance = c(
      model_full$mean$deviance,
      model_agedrought$mean$deviance,
      model_agedensity$mean$deviance,
      model_age$mean$deviance,
      model_drought$mean$deviance,
      model_density$mean$deviance,
      model_null$mean$deviance
    ),
    DIC = c(
      model_full$DIC,
      model_agedrought$DIC,
      model_agedensity$DIC,
      model_age$DIC,
      model_drought$DIC,
      model_density$DIC,
      model_null$DIC
    )
  )
  
  models.df$deltaDIC <-
    models.df$DIC -
    min(models.df$DIC)
  
  models.df <-
    models.df[order(models.df$deltaDIC), ]
  
  models.df
  
  write.csv(models.df, './output/modelselection.csv')

# ---- simulated data ----

#input true data
set.seed(123)

nind  <- 250
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
ni <- 4000
nt <- 10
nb <- 2000
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

# ---- Top Model: age + age.old + density---- 


#keep simulated age.old to 6 
age.old <- age.old.list$knot_6
age.old.sim <- pmax(age.sim - 6, 0) #knot point of 9 from previous model selection


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
  beta2 ~ dnorm(0, 0.001) #age.old
  beta3 ~ dnorm(0, 0.001)
  
  
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
            logit(phi[i,t-1]) <- int + beta1*age[ i , t-1 ]   #age
                                 + beta2*age.old[ i , t-1 ] #age old
                                 + beta3*density[i , t-1 ]  #density
                                  


          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters
   
    for (j in 1:14) { #simulated age
            phi.age[j] <- exp( int + beta1*age.sim[j]
                                      + beta2*age.old.sim[j]
                                      ) /
                               (1 + exp( int + beta1*age.sim[j]
                                                  + beta2*age.old.sim[j]
                                             ) )

    }
    
    phi.slope.6to7 <- (phi.age[7] - phi.age[6]) 
    phi.slope.7to8 <- (phi.age[8] - phi.age[7]) 
    phi.slope.8to9 <- (phi.age[9] - phi.age[8])
    phi.slope.9to14 <- (phi.age[14] - phi.age[9])/5

    # phi.diff <- phi.age[9] - phi.age[1]
    # phi.mean <- (phi.age[1] +phi.age[2] +phi.age[3] +phi.age[4] +
    #             phi.age[5] +phi.age[6] +phi.age[7]) / 7
    # 
    # phi.9to10 <- phi.age[10] - phi.age[9]
    # phi.10to13 <- (phi.age[13] - phi.age[10] )/ 3
   
}
",fill = TRUE)
sink()

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = rnorm(1,0,1), 
  beta2=rnorm(1,0,1)
  
)
}

#Parameters to monitor
parameters <- c('int', 'beta1', 'beta2', 'beta3', 'p', 'phi.slope.6to7', 
                'phi.slope.7to8', 'phi.slope.8to9', 'phi.slope.9to14', 'phi.age')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3


#age   
jags.data <- list(
  h = h,
  ch = ch,
  f = f,
  nind = nrow(ch),
  age = age.raw,
  age.old = age.old,   # age.old for knot = 6
  age.sim = age.sim,
  age.old.sim = age.old.sim,
  density = density.sc
  
)

model_top <- jagsUI( #full model
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

print(model_top)  
write.csv(model_top$summary, './output/topmodelagedens.csv')
# # #---- Top model plot ----

#create a tibble of the posterior draws
gather <- model_top %>%
  spread_draws(phi.age[age])

gather$ageclass <- as.factor(gather$age - 0.5)




phi.plot<- gather %>%
  ggplot(aes(x= ageclass, y=phi.age)) +
  stat_lineribbon(.width = 0.95)+
  geom_vline(xintercept = 9, linetype = "dashed", color = "red", linewidth = 1)+
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1)+
  
  guides(fill = "none")+ #remove legend from ribbon
  scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
  labs(x = "Age", y = "Annual Survival Probability", title = "")+
  scale_x_discrete(breaks = c(.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5)) +
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
        plot.background = element_rect(fill='transparent', color=NA),
        axis.text.x = element_text(
          angle = -45,
          hjust = 0)) #transparent plot bg)
phi.plot
ggsave('./figures/topmodelplot8to9.png', phi.plot)


# ---- summary stats ----
library(dplyr)
data <- read.csv('./cleaned/eydatawfawns.csv', header = T)

#ensure no duplicates
data %>%
  count(animal_id, cap_year) %>%
  filter(n > 1)

#how many capture occasions including cameras
sum(data$status_cam != 0) #how many rows do not contain 0
sum(data$status_nocam !=0) #how many were not supplemented with cameras

#how many unique individuals
data %>% summarise(n = n_distinct(animal_id))

#count individuals in each age class
agecount <- data %>% filter(status_cam %in% c("1",'2')) %>%  
  distinct(animal_id, age) %>%   
  count(age) 
data %>% filter(status_cam %in% c("1",'2')) %>%  
  distinct(animal_id, age) %>%   
  filter(age == 12.5) %>% 
  pull(animal_id)
