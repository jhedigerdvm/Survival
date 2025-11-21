library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

#clean temperature data
maxtemp <- read_csv("raw/temp_maxtemp.csv") #monthly max temperature

maxtemp <- maxtemp %>%  rename(maxtemp = "Texas Maximum Temperature", #rename columns
                               date = "#  Dimmit County")

maxtemp <- maxtemp[-c(1:2),] #remove row 1 and 2 filler columns

#current date format is yyyymm, break into two columns
maxtemp$year  <- substr(maxtemp$date, 1, 4)
maxtemp$month <- substr(maxtemp$date, 5, 6)

maxtemp_filter <- maxtemp %>% filter(year <= 2022 & year >= 2007)

averagetemp <- read_csv("raw/temp_monthlyaverage.csv") #monthly average temperature

averagetemp <- averagetemp %>%  rename(averagetemp = "Texas Average Temperature", #rename columns
                               date = "#  Dimmit County")

averagetemp <- averagetemp[-c(1:2),] #remove row 1 and 2 filler columns

#current date format is yyyymm, break into two columns
averagetemp$year  <- substr(averagetemp$date, 1, 4)
averagetemp$month <- substr(averagetemp$date, 5, 6)

averagetemp_filter <- averagetemp %>% filter(year <= 2022 & year >= 2007)

#add average temp column to max temp dataframe
maxtemp_filter$avgtemp <- averagetemp_filter$averagetemp

temperature <- maxtemp_filter[,c(3:5,2)] #reorder columns 

temperature$maxtemp <- as.numeric(temperature$maxtemp)
temperature$avgtemp <- as.numeric(temperature$avgtemp)

#add celsius columns
temperature$maxtempC <- (temperature$maxtemp - 32) * 5/9
temperature$avgtempC <- (temperature$avgtemp - 32) * 5/9

write.csv(temperature, './cleaned/temperaturedataNOAA.csv', row.names = F)
