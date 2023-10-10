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

mort_data<- read.csv('.//cleaned/mort_data1.csv')
mort_data$animal_id<-as.factor(mort_data$animal_id)
mort_data$birth_site<-as.factor(mort_data$birth_site)
mort_data<- mort_data %>% select(-X)
write.csv(mort_data, file = 'mort_data.csv')
mort_data$death_year<- as.numeric(mort_data$death_year)
mort_data$birth_year<- as.numeric(mort_data$birth_year)

mort_data$age_death<- (mort_data$death_year - mort_data$birth_year)

data('GBSG2', package = "TH.data")
num_cens <- table(GBSG2$cens)
num_cens

sobj<-Surv(GBSG2$time, GBSG2$cens)
sobj[1:10]
summary(sobj)
str(sobj)

cap<-read.csv('.//Raw/cap_hx.csv')

mort_data<-rename(mort_data, last_cap_year = death_year)
mort_data<-rename(mort_data, time = age_death)

##censor individuals with last known capture >2020
mort_data$cens<- ifelse(mort_data$last_cap_year >= '2021', 0, 1)
mort_data$time<- as.integer(mort_data$time)
mort_data$cens<- as.integer(mort_data$cens)

km<- survfit(Surv(time,cens)~ birth_site, data = mort_data)
summary(km)
ggsurvplot(km, linetype = 0, conf.int = TRUE)
?ggsurvplot
ggsurvplot(km, risk.table = TRUE)

write.csv(mort_data, file = "mort_data_clean.csv")

num_cens <- table(mort_data$cens)
num_cens
