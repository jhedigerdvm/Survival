#Survival model for age and birthsite, capyear as a random effect
#treatment and control survival are not different, but TGT survive less 
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

data<- read.csv('./cleaned/caphx.rainfall.long.csv', header = T)


ch <- data$status #create vector with capture histories 

known.fate <- data$status #known fate vector with 2 indentifying deaths associated with capture or harvest

#create capture history with just 1s and 0s, remove 'known fates'
indices <- which(ch == 2, arr.ind = TRUE) #34 individuals with known fates 
ch[indices] <- 1

#create birthsite vector
bs <- as.numeric(factor(data$bs)) # 1 = dmp, 2 = ey, 3 = wy


