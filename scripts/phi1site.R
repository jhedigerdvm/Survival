#survival one sites
# ---- Load Packages ----
#final run of survival analyses
library(jagsUI)
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
# bs.2 <- bs
# bs.2[bs.2 %in% c(3)] <-2 #binary variable of dmp versus pasture
# sum(bs.2==1) #182 individuals in DMP
# sum(bs.2==2) #307 individuals in pasture

sum(bs==1) #control
sum(bs==2) #tgt + dmp
unique(bs)
# 
# #create ageclass matrix with continuous age cov
# data$age.sc <- scale(data$ageclass)
# age.sc<- pivot_wider(data, names_from = 'year', values_from = 'age.sc', id_cols = 'animal_id' )
# age.sc <- as.matrix(age.sc[,-1])


#create ageclass matrix treating age as categorical
data$age.sc <- scale(data$ageclass)
hist(data$age.sc)
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

#add weight and antler vectors
weight<- pivot_wider(data, names_from = 'year', values_from = 'weight', id_cols = 'animal_id' )
weight<- as.matrix(weight[,-1])
weight <- scale(weight) # scale and center


#function for weight matrix
weight.init <- weight
weight.init[is.na(weight.init)]<-0 #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
weight.init[!is.na(weight)]<-NA
for (i in 1:dim(weight.init)){ #cant have mean weight for years before animal was first captured
  weight.init[i,1:f[i]]<-NA
}

#need to find where NA values are in weight matrix, will use this information to build priors
weight <- as.data.frame(weight)
indices <- as.data.frame(which(is.na(weight), arr.ind=T))
indices <- indices %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
NA_indices_weight <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
for(i in 1:nrow(indices)){
  NA_indices_weight[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
}
weight<-as.matrix(weight)

#how many occasions does each individual have of an NA weight
occasions_weight <- rowSums(is.na(weight))


#add simulated weight values
nvalues <- 100
weight.sim <- seq(from = min(weight, na.rm = T), to = max(weight, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 


##now do the same thing we did with weight with antlers
antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcsin', id_cols = 'animal_id' )
antlers<- as.matrix(antlers[,-1])

antlers <- scale(antlers) # scale and center


#function for weight matrix
antlers.init <- antlers
antlers.init[is.na(antlers.init)]<-0 #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
antlers.init[!is.na(antlers)]<-NA
for (i in 1:dim(antlers.init)){ #cant have mean weight for years before animal was first captured
  antlers.init[i,1:f[i]]<-NA
}

#need to find where NA values are in weight matrix, will use this information to build priors
antlers <- as.data.frame(antlers)
indices_antlers <- as.data.frame(which(is.na(antlers), arr.ind=T))
indices_antlers <- indices_antlers %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
NA_indices_antlers <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
for(i in 1:nrow(indices_antlers)){
  NA_indices_antlers[indices_antlers[[i,1]],indices_antlers[[i,3]]] <- indices_antlers[[i,2]]
}
antlers<-as.matrix(antlers)

#how many occasions does each individual have of an NA weight
occasions_antlers <- rowSums(is.na(antlers))
#


#add simulated antler values
nvalues <- 100
antler.sim <- seq(from = min(antlers, na.rm = T), to = max(antlers, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data


#create spring pmdi data
pmdi.spring.sc<-pivot_wider(data, names_from = 'year', values_from = 'pmdi_spring.sc', id_cols = 'animal_id' )
pmdi.spring.sc<-as.matrix(pmdi.spring.sc[,-1])

##create simulated birthyear rainfall vector
nvalues <- 100
pmdi.spring.sc.sim <- seq(from = min(pmdi.spring.sc, na.rm = T), to = max(pmdi.spring.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

unique(bs) # two sites

#add density data
data$densitykm2.sc <- scale(data$densitykm2)
density <-pivot_wider(data, names_from = 'year', values_from = 'densitykm2.sc', id_cols = 'animal_id' )
density <- density[,-1]
density <- as.matrix(density)

##create simulated densities
nvalues <- 100
density.sim <- seq(from = min(density, na.rm = T), to = max(density, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

#age simulation for continuous model
nvalues <- 15
age.sim <- seq(from = min(age.sc, na.rm = T), to = max(age.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
# 

# # ---- Model1: phi ~ site + age + weight + spring pmdi + density ----
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("phi.age.jags")
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
#   beta1[1] <- 0 #age
#   for ( u in 2:15) {
#     beta1[u] ~ dnorm(0, 0.01)  #age
#   }
# 
#   beta2[1] <- 0 # birth site, control
#   beta2[2] ~ dnorm(0,0.001) #birthsite, west yana
#   
#   beta3 ~ dnorm(0,0.001)  #capture year spring pmdi
# 
#   beta4 ~ dnorm(0, 0.001) #density
#   
#   beta5 ~ dnorm(0, 0.001) #weight
#   
# 
#   eps1[1] <- 0 #capture year RE
#    for (u in 2:14){  #prior for year effect
#     eps1[u] ~ dnorm(0,tau)
#   }
# 
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
#             logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
#                                       + beta2[bs[i]]            #birth site, 2 sites
#                                       + beta3*pmdi[i, t-1]   #capture year pmdi spring
#                                       + beta4*density[i,t-1] #population density
#                                       + beta5*morpho[i,t-1]   #morphometric measurement
#                                       + eps1[year[i]]           #capture year random effect
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
#         # # for (i in 1:2){ #birthsite
#         #   for (j in 1:12) { #age
#         #     phi.age[ j] <- exp( int + beta1[j] ) / #+ beta3[i]
#         #                              (1 + exp( int + beta1[j]  ) )#+ beta3[i]
#         # 
#         #     }
# 
#         # for (i in 1:100){ #density sim
#         #   for (j in 1:2) { #site
#         #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
#         #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
#         #     
#         #     }}
#         # for (k in 1:14) { #cap year effect eps1
#         # 
#         #   phi.year[k] <- exp( int + eps1[k]   )/
#         #                     (1 + exp( int + eps1[k] ))
#         # 
#         #     }
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #Function for latent state
# z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }
# 
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
#                   bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
#                   NA_indices = NA_indices_weight, occasions = occasions_weight,
#                   morpho = weight, year = capyear, density = density, density.sim = density.sim)
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1), 
#   z = z.init,
#   beta1 = c(NA, rnorm(14,0,1)),     #age beta
#   beta2 = c(NA, rnorm(1, 0, 1)),    # site beta
#   beta3 = rnorm(1,0,1),             #spring pmdi
#   beta4 =  rnorm(1,0,1),            # density
#   beta5 =  rnorm(1,0,1),            #morpho
#   eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
# )
# }
# 
# 
# parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'beta5', 'eps1')
# 
# # MCMC settings
# ni <- 5000
# nt <- 10
# nb <- 1000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
#                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.age)
# MCMCtrace(phi.age)
# 
# write.csv(phi.age$summary, './output/model1.csv')
# 
# #create a tibble of the posterior draws
# gather<- phi.age %>% gather_draws(phi.year[year]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# gather$year<- as.factor(gather$year)
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=year, y=.value)) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("EAST YANA", "WEST YANA") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("EAST YANA", "WEST YANA"))+ #color of line but no opacification
#   labs(x = "CAPTURE YEAR", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   scale_x_discrete(labels = c(
#     "1" = "2008",
#     "2" = "2009",
#     "3" = "2010",
#     "4" = "2011",
#     "5" = "2012",
#     "6" = "2013",
#     "7" = "2014",
#     "8" = "2015",
#     "9" = "2016",
#     "10" = "2017",
#     "11" = "2018",
#     "12" = "2019",
#     "13" = "2020",
#     "14" = "2021"
#   ))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.year.jpg', phi.plot, width = 10, height = 5)
# 
# 
# 
# #create a tibble of the posterior draws
# gather<- phi.age %>% gather_draws(phi.dens[density, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# # gather$year<- as.factor(gather$year)
# # 
# # #find first row for 2nd rain value
# # first_idx <- which(gather$density == 2)[1] # 4500 values of antler 1
# # # 
# # # #unscale and uncenter weight
# # # morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
# # 
# # #create vector containing simulated morpho data but in the format to sync up with gather
# # vector <- numeric(0)
# # morpho.sim.usc1 <- for (i in density.sim) {
# #   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
# #   vector <- c(vector,rep_i)
# # 
# # }
# 
# # gather$bodymass <- vector
# 
# # gather$bodymass <- gather$bodymass/2.2
# 
# #plot for average age individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=density, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   # stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "DENSITY", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.densityxsite.noagenopmdi.jpg', phi.plot, width = 10, height = 10)


# #Prepare to plot phi.rain
# gather<- phi.weight %>% gather_draws(phi.drought[pmdi, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$pmdi == 2)[1] # 4500 values of antler 1
# 
# #unscale and uncenter weight
# pmdi.sim.usc <- (pmdi.spring.sc.sim * sd(data$pmdi_spring, na.rm = T)) + mean(data$pmdi_spring, na.rm = T)
# 
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# rain.sim.usc1 <- for (i in pmdi.sim.usc) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$pmdispring <- vector
# 
# 
# #plot for 1.5 YO individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=pmdispring, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "SPRING DROUGHT INDEX", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         #axis.text.x = element_text(margin = margin(t = 5)),
#         panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.PMDI.2sites.jpg', phi.plot, width = 10, height = 10)
# 
# 
# 
# #Prepare to plot phi.age
# gather<- phi.weight %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
# 
# gather$age <- as.factor(gather$age)
# 
# gather <- gather %>% filter(!age %in% c(13,14,15))
# 
# #plot for average age individual
# 
# phi.plot<- gather %>% 
#   ggplot(aes(x=age, y=.value )) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.1,0.1),          # x, y inside the plot area
#         legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         #axis.text.x = element_text(margin = margin(t = 5)),
#         panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.AGE.1site.jpg', phi.plot, width = 10, height = 10)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ---- Model2: phi ~ antlers + site + age + spring pmdi ----
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("phi.weight.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#   int ~ dnorm(0, 0.001)
#   beta1[1] <- 0 #age
#   beta2[1] <- 0  #site
#   eps1[1] <- 0 #capture year RE
#   
#   beta2[2] ~ dnorm(0, 0.001)  #site #only two sites for now
#   beta3 ~ dnorm(0,0.001)  #capture year spring pmdi
#   beta4 ~ dlnorm(0, 0.01)    # morpho beta
# 
#   
#   for ( u in 2:15) {
#     beta1[u] ~ dnorm(0, 0.01)  #age
#   }
#   
# 
# 
#   
#   for (u in 1:nind){      #prior for missing morphometrics
#     for (j in 1:occasions[u]){
#     morpho[u,NA_indices[u,j]] ~ dnorm( 0, 0.01)
#        }
#   }
# 
#   for (u in 2:14){  #prior for year effect
#     eps1[u] ~ dnorm(0,tau)
#   }
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
#             logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
#                                       + beta2[bs[i]]            #birth site
#                                       + beta3*pmdi[i, t-1]   #capture year pmdi spring
#                                       + beta4* morpho[i, t-1]   #morphology
#                                       + eps1[year[i]]           #capture year random effect
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
#    # #derived parameters
#       for (i in 1:100 ) { #weight simulation, beta3
#       for (j in 1:2){ #site, beta2
# 
# 
#       phi.weight[i, j] <- exp( int+ beta4*morpho.sim[i]  + beta2[j]  )/
#                             (1 + exp( int + beta4*morpho.sim[i]  + beta2[j]))
# 
#     } # for j
#     } # for l
# 
#       for (i in 1:100 ) { #rain simulation, beta3
#       for (j in 1:2){ #site, beta2
# 
#       phi.drought[i, j] <- exp( int + beta3*pmdi.sim[i]  + beta2[j]  )/
#                             (1 + exp( int + beta3*pmdi.sim[i]  + beta2[j]))
# 
#     } # for j
#     } # for l
# 
#     for (i in 1:12 ) { #age beta1
#     for (j in 1:2){ #site, beta2
# 
#       phi.age[i, j] <- exp( int+ beta1[i]  + beta2[j]  )/
#                             (1 + exp( int + beta1[i]  + beta2[j]))
# 
#     } # for j
#     } # for l
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #Function for latent state
# z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }
# 
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
#                   bs = bs, morpho.sim = antler.sim, pmdi.sim = pmdi.spring.sc.sim,
#                   NA_indices = NA_indices_antlers, occasions = occasions_antlers,
#                   morpho = antlers, year = capyear)
# 
# # Initial values
# inits <- function(){list(
#   int = rnorm(1,0,1), 
#   z = z.init,
#   morpho = antlers.init, #initial values for NA morphos
#   beta1 = c(NA, rnorm(14,0,1)), #age beta
#   beta2 = c(NA, rnorm(1, 0, 1)),#,#site beta
#   beta3 = rnorm(1, 0, 1), # pmdi beta
#   beta4 = rlnorm(1,0,1),#,#morpho 
#   eps1 = c(NA, rnorm(13, 0, 1)) #capture year random effect
# )
# }
# 
# 
# parameters <- c('int', 'beta1','beta2', 'beta3', 'beta4', 'eps1', 'phi.weight', 'phi.drought', 'phi.age')
# 
# # MCMC settings
# ni <- 5000
# nt <- 10
# nb <- 1000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.weight <- jagsUI(jags.data, inits, parameters, "phi.weight.jags", n.chains = nc,
#                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.weight)
# MCMCtrace(phi.weight)
# write.csv(phi.weight$summary, './output/phi.antlers.2sites.csv')
# 
# #create a tibble of the posterior draws
# gather<- phi.weight %>% gather_draws(phi.weight[antlers, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$antlers == 2)[1] # 4500 values of antler 1
# 
# #unscale and uncenter weight
# morpho.sim.usc <- (antler.sim * sd(data$bcsin, na.rm = T)) + mean(data$bcsin, na.rm = T)
# 
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# morpho.sim.usc1 <- for (i in morpho.sim.usc) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$bcsin <- vector
# 
# # gather$bodymass <- gather$bodymass/2.2
# 
# #plot for average age individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=bcsin, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "ANTLER SCORE (IN)", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.bcs.in.2sites.jpg', phi.plot, width = 10, height = 10)
# 
# 
# #Prepare to plot phi.rain
# gather<- phi.weight %>% gather_draws(phi.drought[pmdi, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$pmdi == 2)[1] # 4500 values of antler 1
# 
# #unscale and uncenter weight
# pmdi.sim.usc <- (pmdi.spring.sc.sim * sd(data$pmdi_spring, na.rm = T)) + mean(data$pmdi_spring, na.rm = T)
# 
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# rain.sim.usc1 <- for (i in pmdi.sim.usc) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$pmdispring <- vector
# 
# 
# #plot for 1.5 YO individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=pmdispring, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "SPRING DROUGHT INDEX", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         #axis.text.x = element_text(margin = margin(t = 5)),
#         panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.PMDI.2sites.jpg', phi.plot, width = 10, height = 10)

# 
# 
# #Prepare to plot phi.age
# gather<- phi.weight %>% gather_draws(phi.age[age, site]) #this creates a dataframe in long format with indexing
# gather$site <- as.factor(gather$site)
# 
# gather$age <- as.factor(gather$age)
# 
# 
# #plot for average age individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=.value, color = site, fill = site)) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.1,0.1),          # x, y inside the plot area
#         legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         #axis.text.x = element_text(margin = margin(t = 5)),
#         panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.AGE.2sites.jpg', phi.plot, width = 10, height = 10)


# ---- Model5: phi ~  age (continuous) + capyear RE TOP MODEL ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)

  # beta1[1] <- 0 #age
  # for ( u in 2:15) {
  #   beta1[u] ~ dnorm(0, 0.01)  #age
  # }


  # eps1[1] <- 0 #capture year RE
  #  for (u in 2:14){  #prior for year effect
  #   eps1[u] ~ dnorm(0,tau.year)
  # }


  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)
  
  # tau.year <- 1/(sigma.year*sigma.year)
  # sigma.year ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <-  int + beta1*ageclass[i,t-1]  #age categorical
                                      #+ eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

          for (j in 1:15) { #age
            phi.age[j] <- exp( int + beta1*age.sim[j] ) / #+ beta3[i]
                                     (1 + exp( int + beta1*age.sim[j]  ) )#+ beta3[i]

            }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        # 
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        # 
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight, age.sim = age.sim,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  morpho = weight.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1',  'eps1','phi.age', 'phi.year')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
#
# write.csv(phi.age$summary, './output/model5.csv')



# #Posterior draws of survival by year
# gather<- phi.age %>% gather_draws(phi.year[year]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# gather$year<- as.factor(gather$year)
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=year, y=.value)) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("EAST YANA", "WEST YANA") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("EAST YANA", "WEST YANA"))+ #color of line but no opacification
#   labs(x = "CAPTURE YEAR", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   scale_x_discrete(labels = c(
#     "1" = "2008",
#     "2" = "2009",
#     "3" = "2010",
#     "4" = "2011",
#     "5" = "2012",
#     "6" = "2013",
#     "7" = "2014",
#     "8" = "2015",
#     "9" = "2016",
#     "10" = "2017",
#     "11" = "2018",
#     "12" = "2019",
#     "13" = "2020",
#     "14" = "2021"
#   ))+
#   scale_y_continuous(expand = expansion(mult = c(.1, 0.2)),
#                      breaks = seq(0.9, 1.4, by = 0.02))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.year.jpg', phi.plot, width = 10, height = 5)



#create a tibble of the posterior draws by age class
gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
gather$age <- as.factor(gather$age)
# gather$year<- as.factor(gather$year)
#
# #find first row for 2nd rain value
# first_idx <- which(gather$density == 2)[1] # 4500 values of antler 1
# #
# # #unscale and uncenter weight
# # morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
#
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# morpho.sim.usc1 <- for (i in density.sim) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#
# }

# gather$bodymass <- vector

# gather$bodymass <- gather$bodymass/2.2

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=age, y=.value)) +
  stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
  # scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
  scale_x_discrete(labels = c(
    "1" = "1.5",
    "2" = "2.5",
    "3" = "3.5",
    "4" = "4.5",
    "5" = "5.5",
    "6" = "6.5",
    "7" = "7.5",
    "8" = "8.5",
    "9" = "9.5",
    "10" = "10.5",
    "11" = "11.5",
    "12" = "12.5"

  ))+
  labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 22),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/PHI.AGE.JPG', phi.plot, width = 10, height = 10)

#

# ---- Model5: phi ~  age (categorical) + capyear RE TOP MODEL ----

# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  int ~ dnorm(0, 0.001)

  beta1[1] <- 0 #age
  for ( u in 2:15) {
    beta1[u] ~ dnorm(0, 0.001)  #age
  }


  eps1[1] <- 0 #capture year RE
   for (u in 2:14){  #prior for year effect
    eps1[u] ~ dnorm(0,tau.year)
  }


  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)
  
  tau.year <- 1/(sigma.year*sigma.year)
  sigma.year ~ dunif(0,100)


# Likelihood
for (i in 1:nind){
   # Define latent state at first capture, we know for sure the animal is alive
      z[i,f[i]] <- 1

      for (t in (f[i]+1):h[i]){
        # State process
            z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
            logit(phi[i,t-1]) <-  int + beta1[ageclass[i,t-1]]  #age categorical
                                      + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

          for (j in 1:15) { #age
            phi.age[j] <- exp( int + beta1[j] ) / #+ beta3[i]
                                     (1 + exp( int + beta1[j]  ) )#+ beta3[i]

          }
          
        
          for (j in 2:15) {
            phi.1diff[j] <- phi.age[j] - phi.age[j-1]
          }
          
          for (j in 3:15) {
            phi.2diff[j] <- phi.age[j] - phi.age[j-2]
          }
          
          phi.average <- (phi.age[1] + phi.age[2] + phi.age[3] + phi.age[4] + phi.age[5] + phi.age[6])/6

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        # 
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        # 
        #     }

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
jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = ageclass, pmdi = pmdi.spring.sc,
                  bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim,
                  NA_indices = NA_indices_weight, occasions = occasions_weight, age.sim = age.sim,
                  morpho = weight, year = capyear, density = density, density.sim = density.sim)

# Initial values
inits <- function(){list(
  int = rnorm(1,0,1),
  z = z.init,
  beta1 = c(NA, rnorm(14,0,1)),     #age beta
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1',  'phi.age', 'phi.average', 'phi.1diff', 'phi.2diff')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
#
# write.csv(phi.age$summary, './output/model5.csv')



# #Posterior draws of survival by year
# gather<- phi.age %>% gather_draws(phi.year[year]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# gather$year<- as.factor(gather$year)
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=year, y=.value)) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("EAST YANA", "WEST YANA") ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo', labels = c("EAST YANA", "WEST YANA"))+ #color of line but no opacification
#   labs(x = "CAPTURE YEAR", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   scale_x_discrete(labels = c(
#     "1" = "2008",
#     "2" = "2009",
#     "3" = "2010",
#     "4" = "2011",
#     "5" = "2012",
#     "6" = "2013",
#     "7" = "2014",
#     "8" = "2015",
#     "9" = "2016",
#     "10" = "2017",
#     "11" = "2018",
#     "12" = "2019",
#     "13" = "2020",
#     "14" = "2021"
#   ))+
#   scale_y_continuous(expand = expansion(mult = c(.1, 0.2)),
#                      breaks = seq(0.9, 1.4, by = 0.02))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.year.jpg', phi.plot, width = 10, height = 5)



#create a tibble of the posterior draws by age class
gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
gather$age <- as.factor(gather$age)
# gather$year<- as.factor(gather$year)
#
# #find first row for 2nd rain value
# first_idx <- which(gather$density == 2)[1] # 4500 values of antler 1
# #
# # #unscale and uncenter weight
# # morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
#
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# morpho.sim.usc1 <- for (i in density.sim) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#
# }

# gather$bodymass <- vector

# gather$bodymass <- gather$bodymass/2.2

#plot for average age individual

phi.plot<- gather %>%
  ggplot(aes(x=age, y=.value)) +
  stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
  # scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
  scale_x_discrete(labels = c(
    "1" = "1.5",
    "2" = "2.5",
    "3" = "3.5",
    "4" = "4.5",
    "5" = "5.5",
    "6" = "6.5",
    "7" = "7.5",
    "8" = "8.5",
    "9" = "9.5",
    "10" = "10.5",
    "11" = "11.5",
    "12" = "12.5"
    
  ))+
  labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 22),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
phi.plot
ggsave('./figures/PHI.AGE.JPG', phi.plot, width = 10, height = 10)

#


#---- Model: age quadratic ----


# Specify model in JAGS language
set.seed(100)
sink("phi.age.jags")
cat("
model {

#prior for recapture prob
p ~ dbeta(1, 1)


#priors
  # int ~ dnorm(0, 0.001)

  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0,0.001)
  # beta1[1] <- 0 #age
  # for ( u in 2:15) {
  #   beta1[u] ~ dnorm(0, 0.01)  #age
  # }


  #eps1[1] <- 0 #capture year RE
   # for (u in 1:14){  #prior for year effect
   #  eps1[u] ~ dnorm(0,tau.year)
   # }

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
            logit(phi[i,t-1]) <-  beta1*ageclass[i,t-1]   #age categorical
                                  + (beta2*ageclass[i,t-1]*ageclass[i,t-1])
                                     # + eps1[year[i]]           #capture year random effect

          # Observation process
            ch[i,t] ~ dbern(mu2[i,t])
            mu2[i,t] <- p * z[i,t]



      } #t
   } #i

   #derived parameters

          for (j in 1:100) { #age
            phi.age[j] <- exp(  beta1*age.sim[j] + beta2*age.sim[j]*age.sim[j] ) / #+ beta3[i]
                                     (1 + exp(  beta1*age.sim[j] + beta2*age.sim[j]*age.sim[j] ) )#+ beta3[i]

            }

        # for (i in 1:100){ #density sim
        #   for (j in 1:2) { #site
        #     phi.dens[i, j] <- exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) / #
        #                              (1 + exp( int + beta3[j] + beta4*density.sim[i] + beta5[j]*density.sim[i] ) ) #
        #
        #     }}
        # for (k in 1:14) { #cap year effect eps1
        #
        #   phi.year[k] <- exp( int + eps1[k]   )/
        #                     (1 + exp( int + eps1[k] ))
        #
        #     }

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
  morpho = weight.init,
  beta1 = rnorm(1,0,1), # c(NA, rnorm(14,0,1)),     #age beta
  beta2=rnorm(1,0,1),
  eps1 = c(NA, rnorm(13, 0, 1))     #capture year random effect
)
}


parameters <- c('int', 'beta1', 'beta2',  'eps1','phi.age', 'age.sim')

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(phi.age)
MCMCtrace(phi.age)
# # 
# # write.csv(phi.age$summary, './output/model5.csv')
# # 
# 
# 
# #create a tibble of the posterior draws
# gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# # gather$year<- as.factor(gather$year)
# 
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=.value)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   #   scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   #   scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   #     labs(x = "CAPTURE YEAR", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   # scale_x_discrete(labels = c(
#   #   "1" = "2008",
#   #   "2" = "2009",
#   #   "3" = "2010",
#   #   "4" = "2011",
#   #   "5" = "2012",
#   #   "6" = "2013",
#   #   "7" = "2014",
#   #   "8" = "2015",
#   #   "9" = "2016",
#   #   "10" = "2017",
#   #   "11" = "2018",
#   #   "12" = "2019",
#   #   "13" = "2020",
#   #   "14" = "2021"
#   # ))+
#   # scale_y_continuous(expand = expansion(mult = c(.1, 0.2)),
#   #                    breaks = seq(0.9, 1.4, by = 0.02))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# # ggsave('./figures/phi.year.jpg', phi.plot, width = 10, height = 5)
# # 
# # 
# # 
# # #create a tibble of the posterior draws
# # gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
# # gather$age <- as.factor(gather$age)
# # # gather$year<- as.factor(gather$year)
# # #
# # # #find first row for 2nd rain value
# # # first_idx <- which(gather$density == 2)[1] # 4500 values of antler 1
# # # #
# # # # #unscale and uncenter weight
# # # # morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
# # #
# # # #create vector containing simulated morpho data but in the format to sync up with gather
# # # vector <- numeric(0)
# # # morpho.sim.usc1 <- for (i in density.sim) {
# # #   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
# # #   vector <- c(vector,rep_i)
# # #
# # # }
# # 
# # # gather$bodymass <- vector
# # 
# # # gather$bodymass <- gather$bodymass/2.2
# # 
# # #plot for average age individual
# # 
# # phi.plot<- gather %>%
# #   ggplot(aes(x=age, y=.value)) +
# #   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
# #   # scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
# #   # scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
# #   scale_x_discrete(labels = c(
# #     "1" = "1.5",
# #     "2" = "2.5",
# #     "3" = "3.5",
# #     "4" = "4.5",
# #     "5" = "5.5",
# #     "6" = "6.5",
# #     "7" = "7.5",
# #     "8" = "8.5",
# #     "9" = "9.5",
# #     "10" = "10.5",
# #     "11" = "11.5",
# #     "12" = "12.5"
# #     
# #   ))+
# #   labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
# #   theme_bw() +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         panel.border = element_blank(),
# #         axis.line = element_line(),
# #         legend.position = "inside",
# #         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
# #         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
# #         legend.text = element_text(size = 28),
# #         legend.title = element_blank(),
# #         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
# #         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
# #         axis.text = element_text(face='bold',size = 22),
# #         # axis.text.x = element_text(angle = 45, hjust = 1),
# #         panel.background = element_rect(fill='transparent'), #transparent panel bg
# #         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# # phi.plot
# # ggsave('./figures/PHI.AGE.JPG', phi.plot, width = 10, height = 10)
# # 
# # # 
# 
# #---- Model: age gompertz ----
# 
# 
# # Specify model in JAGS language
# set.seed(100)
# sink("phi.age.jags")
# cat("
# model {
# 
# #prior for recapture prob
# p ~ dbeta(1, 1)
# 
# 
# #priors
#  
#   lambda ~ dunif(0,5)
#   #gamma ~ dnorm(0,0.001)
#   beta0 ~ dnorm(0, 0.001)
#   beta1[1] <- 0
#   beta1[2] ~ dnorm(0,0.001)
# 
# 
# # Likelihood
# for (i in 1:nind){
#       #gamma for site 
#           gamma[i] <- beta0 + beta1[bs[i]]
#             
#       # Define latent state at first capture, we know for sure the animal is alive
#           z[i,f[i]] <- 1
# 
#       for (t in (f[i]+1):h[i]){
#         
#         
#       
#         #Gompertz Hazard
#             y[i,t-1] <- lambda * exp(-gamma[i]*ageclass[i, t-1] )
#             
#         # State process
#             z[i,t] ~ dbern(mu1[i,t]) #toss of a coin whether individual is alive or not detected
#             mu1[i,t] <- phi[i,t-1] * z[i,t-1]  #t-1 because we are looking ahead to see if they survived from 1 to 2 based upon them being alive at 2
#             phi[i,t-1] <-  exp(-y[i, t-1])
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
# 
#    #derived parameters
#     
#     for (i in 1:2){
#         gamma.site[i] <- beta0 + beta1[bs[i]]
# 
#       for (j in 1:15){
#        phi.age[i,j] <- exp(-(lambda * exp(-gamma.site[i] * age.sim[j])))
#     }
#   
#   }
#    
#    
#         # for (i in 1:2){ #bs
#         #   for (j in 1:15) { #age
#         #     phi.age[i,j] <- exp(-(lambda * exp(-gamma[bs[i]]*age.sim[j])))
#         # 
#         #   }}
#           
#           #
#           # for (j in 2:15) {
#           #   phi_diff[j] <- phi.age[j] - phi.age[j-1]
#           # 
#           # }
#           # 
#           # 
#           # 
#           # age_decline <- log(lambda) / gamma
# 
#       
# 
# }
# ",fill = TRUE)
# sink()
# 
# 
# #Function for latent state
# z.init <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
# 
# for(i in 1:dim(z.init)[1]){
#   z.init[i, f[i]:h[i]] <- 1
#   z.init[i,f[i]] <- NA
# }
# 
# 
# # Bundle data
# jags.data <- list(h = h, ch = ch, f = f, nind = nrow(ch), ageclass = age.sc, pmdi = pmdi.spring.sc,
#                   bs = bs, morpho.sim = weight.sim, pmdi.sim = pmdi.spring.sc.sim, age.sim = age.sim,
#                   NA_indices = NA_indices_weight, occasions = occasions_weight,
#                   morpho = weight, year = capyear, density = density, density.sim = density.sim)
# 
# # Initial values
# inits <- function(){list(
#   beta0 = rnorm(1,0,1),
#   z = z.init,
#   beta1 = c(NA, rnorm(1,0,1)),
#   lambda = runif(1,0,5)
#   #gamma = rnorm(1,0,1)#,
#   
# )
# }
# 
# 
# parameters <- c('beta0', 'beta1', 'phi.age', 'age_decline', 'phi_diff', 'lambda'
#                 )
# 
# # MCMC settings
# ni <- 5000
# nt <- 10
# nb <- 1000
# nc <- 3
# 
# # Call JAGS from R (BRT 3 min)
# phi.age<- jagsUI(jags.data, inits, parameters, "phi.age.jags", n.chains = nc,
#                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# 
# print(phi.age)
# MCMCtrace(phi.age)
# # 
# write.csv(phi.age$summary, './output/phi.gompertz.csv')
# # 
# 
# 
# #create a tibble of the posterior draws
# gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
# # gather$site <- as.factor(gather$site)
# # gather$year<- as.factor(gather$year)
# 
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$age == 2)[1] # 4500 values of antler 1
# 
# #create vector containing simulated morpho data but in the format to sync up with gather
# vector <- numeric(0)
# age.sim1 <- for (i in age.sim) {
#   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
#   vector <- c(vector,rep_i)
#   
# }
# 
# gather$ageclass <- vector
# 
# #unscale and center
# gather$ageclass <-  (gather$ageclass * sd(ageclass, na.rm = T)) + mean(ageclass, na.rm = T)
# senescence <- 2.07 * sd(ageclass, na.rm = T) + mean(ageclass, na.rm = T)
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=.value)) +
#   stat_lineribbon(.width = 0.95)+
#   scale_fill_viridis_d(option = 'turbo', alpha = .2 ) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(option = 'turbo')+ #color of line but no opacification
#   labs(x = "Age", y = "Annual Survival Probability", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "none",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 18, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/phi.year.jpg', phi.plot, width = 10, height = 5)
# 
# 
# 
# #create a tibble of the posterior draws
# gather<- phi.age %>% gather_draws(phi.age[age]) #this creates a dataframe in long format with indexing
# gather$age <- as.factor(gather$age)
# # gather$year<- as.factor(gather$year)
# #
# # #find first row for 2nd rain value
# # first_idx <- which(gather$density == 2)[1] # 4500 values of antler 1
# # #
# # # #unscale and uncenter weight
# # # morpho.sim.usc <- (weight.sim * sd(data$weight, na.rm = T)) + mean(data$weight, na.rm = T)
# #
# # #create vector containing simulated morpho data but in the format to sync up with gather
# # vector <- numeric(0)
# # morpho.sim.usc1 <- for (i in density.sim) {
# #   rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
# #   vector <- c(vector,rep_i)
# #
# # }
# 
# # gather$bodymass <- vector
# 
# # gather$bodymass <- gather$bodymass/2.2
# 
# #plot for average age individual
# 
# phi.plot<- gather %>%
#   ggplot(aes(x=age, y=.value)) +
#   stat_pointinterval( position = position_dodge(width=0.5))+ #statline ribbon takes posterior estimates and calculates CRI
#   # scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("CONTROL", "TREATMENT") ) + #this allowed me to opacify the ribbon but not the line
#   # scale_color_viridis_d(option = 'turbo', labels = c("CONTROL", "TREATMENT"))+ #color of line but no opacification
#   scale_x_discrete(labels = c(
#     "1" = "1.5",
#     "2" = "2.5",
#     "3" = "3.5",
#     "4" = "4.5",
#     "5" = "5.5",
#     "6" = "6.5",
#     "7" = "7.5",
#     "8" = "8.5",
#     "9" = "9.5",
#     "10" = "10.5",
#     "11" = "11.5",
#     "12" = "12.5"
#     
#   ))+
#   labs(x = "AGE CLASS", y = "ANNUAL SURVIVAL PROBABILITY", title = "")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
#         legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         legend.title = element_blank(),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 22),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# phi.plot
# ggsave('./figures/PHI.AGE.JPG', phi.plot, width = 10, height = 10)
# 
# # 
