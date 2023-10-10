library(dplyr)
library(tidyr)
library(magrittr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ciTools)
library(ggpubr)
library(rstatix)
library(broom)
library(survival)
install.packages('survminer')
library(survminer)
library(lubridate)

data<-read.csv('.//Raw/cap_log.csv')
agr<- aggregate(data$cap_year, list(data$animal_id), FUN = max)
data$cap_year<-as.numeric(data$cap_year)

#add 1 year to capyear

for (i in 1:nrow(agr)) {
  x <- substr(agr$x[i], 1, 1)
  if (x == "2") {agr$event[i] <- "1"} 
}

agr$event <- as.numeric(agr$event)
agr$death_year <- rowSums(agr[ , c(2,3)], na.rm=TRUE)
agr$death_year<-as.character(agr$death_year)

agr<-rename(agr, animal_id = Group.1)
agr$animal_id<-as.factor(agr$animal_id)
mort_data <-agr[,c(1,4)]


#create birth_sites with for loop
for (i in 1:nrow(mort_data)) {
  x <- substr(mort_data$animal_id[i], 1, 1)
  if (x == "2") {mort_data$wyey[i] <- "west_yana"} else 
  { mort_data$wyey[i] <- "east_yana"}
}
##make a new column with discernment of wy and pasture born
for (i in 1:nrow(mort_data)) {
  x <- substr(mort_data$animal_id[i], 10, 10)
  if (x == "1") {mort_data$dmp_pasture[i] <- "dmp"} else 
  { mort_data$dmp_pasture[i] <- "pasture"}
}

##1 vector with 3 class variable 
mort_data$birth_site <- paste(mort_data$wyey, mort_data$dmp_pasture, sep = "_")
mort_data$birth_site<-as.factor(mort_data$birth_site)

#create birth yeaer
for (i in 1:nrow(mort_data)) {
  x <- substr(mort_data$animal_id[i], 5, 8)
  if (x == "2007") {mort_data$birth_year[i] <- "2007"} else 
  { mort_data$wyey[i] <- "east_yana"}
}

write.csv(mort_data, file = "mort_data.csv")
mort_data<-read.csv('.//mort_data.csv')

mort_data1<-mort_data[,c(2,3,6,7)]
write.csv(mort_data1, file = "mort_data1.csv")
