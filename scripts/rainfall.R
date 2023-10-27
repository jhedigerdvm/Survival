library(tidyverse)
library(here)
#rainfall data
data<- read.csv('./raw/Rainfall.csv', header = T) #does not contain 2022 data

#remove JAnuary 2022 and remaining NA rows
data<- data[-c(241:252),]
#get into long form with column site
east.rain<- data[,-4]
west.rain<- data[,-3]
west.rain$site<- "wy" #add new column indicating site
west.rain <- rename(west.rain, rain=WEST.YANA) #renaming column to rain
east.rain$site<- 'ey'
east.rain <- rename(east.rain, rain=EAST.YANA)
data1<- rbind(west.rain,east.rain)
cor(west.rain$rain, east.rain$rain) #92% correlation between east and west yana

totalrain<- data1 %>% group_by(Year, site ) %>% summarise(total.rain=sum(rain),
                                                          monthly.mean=mean(rain)) #sum is total annual rainfall, mean is monthly avg
#filter for only march april and may rainfall
march.apr.may <- data1 %>%
  filter(Month %in% c('Mar', 'Apr', 'May'))
#group by year to sum march april and may rainfall
total.march.apr.may<- march.apr.may %>% group_by(Year,site) %>% summarise(sum.march.apr.may=sum(rain),
                                                                              mean.march.apr.may=mean(rain))
  
#filter for only june july aug rainfall
jun.jul.aug <- data1 %>%
  filter(Month %in% c('Jun', 'Jul', 'Aug'))

#group by year to sum march april and may rainfall
total.jun.jul.aug<- jun.jul.aug %>% group_by(Year,site) %>% summarise(sum.jun.jul.aug=sum(rain),
                                                                          mean.jun.jul.aug=mean(rain))

#filter for only march april and may rainfall
sept.oct <- data1 %>%
  filter(Month %in% c('Sept', 'Oct'))
#group by year to sum march april and may rainfall
total.sept.oct<- sept.oct %>% group_by(Year,site) %>% summarise(sum.sept.oct=sum(rain),
                                                                          mean.sept.oct=mean(rain))
#filter for only march april and may rainfall
oct <- data1 %>%
  filter(Month %in% c('Oct'))
#group by year to sum march april and may rainfall
total.oct<- oct %>% group_by(Year,site)
total.oct <- rename(total.oct, sum.oct = rain)

rainfall<- total.jun.jul.aug
rainfall$sum.march.apr.may <- total.march.apr.may$sum.march.apr.may
rainfall$sum.sept.oct <- total.sept.oct$sum.sept.oct
rainfall$oct <- total.oct$sum.oct
rainfall$annual<- totalrain$total.rain
rainfall <- rainfall[,-4]

write.csv(rainfall, './cleaned/rainfall_clean.csv', row.names = F)
