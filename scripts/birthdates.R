#code to clean DMP birthdates
library(tidyverse)
library(here)

#clean dmp birthdates
data<- read.csv('./raw/dmp_birthdates.csv', header = T)
data$date<- as.Date(data$date, "%m/%d/%Y") #convert dates from 10/1/2007 to standard format


data <- data[!is.na(data$date),] #remove NA from date

data$monthday <- format(data$date, '%m-%d') #create column with just month and date, exclude year
barplot(table(data$monthday)) #plot birthdates to see distribution

#create column with 1 for fawns born before august, and 2 for fawns born before august
cutoff = '07-31'
data$early.late <- ifelse(data$monthday > cutoff, 2, 1)

write.csv(data, './cleaned/birthdates.csv', row.names = F)

#load capture data including fawns
cap<-read.csv('./cleaned/capture_cleaned_fawns_long.csv', header = T)

#filter for individuals born in DMP
cap<- cap %>% filter(bs == 'dmp')
