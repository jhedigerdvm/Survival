#total rainfall now starts in november and ends in october
library(tidyverse)
library(here)
#rainfall data
data<- read.csv('./raw/Rainfall.csv', header = T) 
data$surv.year<- data$surv.year - 1

#get into long form with column site
east.rain<- data[,-4]
west.rain<- data[,-3]
west.rain$site<- "wy" #add new column indicating site
west.rain <- rename(west.rain, rain=WEST.YANA) #renaming column to rain
east.rain$site<- 'ey'
east.rain <- rename(east.rain, rain=EAST.YANA)
data1<- rbind(west.rain,east.rain)
data1$surv.year <- as.factor(data1$surv.year)
data1$site <-as.factor(data1$site)
cor(west.rain$rain, east.rain$rain) #92% correlation between east and west yana

totalrain<- data1 %>% group_by(surv.year, site ) %>% summarise(total.rain=sum(rain)) #sum is total annual rainfall from nov to oct
totalrain<- totalrain[-c(43:44),]#remove 2022

#filter for only march april and may rainfall
march.apr.may <- data1 %>%
  filter(Month %in% c('Mar', 'Apr', 'May'))
#group by year to sum march april and may rainfall
total.march.apr.may<- march.apr.may %>% group_by(surv.year,site) %>% summarise(sum.march.apr.may=sum(rain))
  
#filter for only june july aug rainfall
jun.jul.aug <- data1 %>%
  filter(Month %in% c('Jun', 'Jul', 'Aug'))

#group by year to sum march april and may rainfall
total.jun.jul.aug<- jun.jul.aug %>% group_by(surv.year,site) %>% summarise(sum.jun.jul.aug=sum(rain))

#filter for only march april and may rainfall
sept.oct <- data1 %>%
  filter(Month %in% c('Sept', 'Oct'))
#group by year to sum march april and may rainfall
total.sept.oct<- sept.oct %>% group_by(surv.year,site) %>% summarise(sum.sept.oct=sum(rain))
#filter for only march april and may rainfall
oct <- data1 %>%
  filter(Month %in% c('Oct'))
#group by year to sum march april and may rainfall
total.oct<- oct %>% group_by(surv.year,site)
total.oct <- rename(total.oct, sum.oct = rain)

rainfall<- total.jun.jul.aug
rainfall$sum.march.apr.may <- total.march.apr.may$sum.march.apr.may
rainfall$sum.sept.oct <- total.sept.oct$sum.sept.oct
rainfall$oct <- total.oct$sum.oct
rainfall$annual<- totalrain$total.rain
rainfall <- rainfall[,-4]

write.csv(rainfall, './cleaned/rainfall_clean1.csv', row.names = F)
