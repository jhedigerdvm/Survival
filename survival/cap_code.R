#this code will take long format and put it into capture history in wide format
library(dplyr)
library(tidyr)
library(magrittr)
library(tidyverse)
library(here)

data<-read.csv("./raw/mortality_data_22.csv", header = T)

# Rename the variable
data$status[data$status == 'Captured'] <- 1
data$status[data$status == 'Nat Mortality'] <- 0
data$status[data$status == 'Alive-Cuddy'] <- 1
data$status[data$status == 'Harvested'] <- 2
data$status[data$status == 'Dead-Capture'] <- 2
unique(data$status)
data$check<-duplicated(data) #check for duplicates
#270-2010-1-019 
data[402,2] <- '2010'
#270-2011-1-215 
data<-data[-611,]
#270-2019-1-909
data<-data[-1917,]
#90-2020-0-2080
data<-data[-3538,]
data[314,2] <- '2010'
data[513,2] <- '2011'

#correct birthyear
for (i in c(1:2037, 3613:3619)){
  x <- substr(data$animal_id[i], 5, 8)
  data$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}

for (i in 2038:3612){
  x <- substr(data$animal_id[i], 4, 7)
  data$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}

write.csv(data,'cap_code_long.csv', row.names = F)

data<- read.csv('cap_code_long.csv', header = T)
data$check<-duplicated(data) #check for duplicates
data<- data[-c(404,1916,3543),-c(6)]
data2<-data %>% pivot_wider(names_from = cap_year, values_from = status, values_fill = 0, id_cols = animal_id)#
#lost birth_year with pivot wider
# birthyear
for (i in c(1:457)){
  x <- substr(data2$animal_id[i], 5, 8)
  data2$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}

for (i in 458:730){
  x <- substr(data2$animal_id[i], 4, 7)
  data2$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}


#add birth site

for (i in 1:nrow(data2)) {
  x <- substr(data2$animal_id[i], 1, 1)
  if (x == "2") {data2$bs[i] <- "wy"} else 
  { data2$bs[i] <- "ey"}
}
##make a new column with discernment of DMP and pasture born
for (i in 1:nrow(data2)) {
  x <- substr(data2$animal_id[i], 10, 10)
  if (x == "1") {data2$bs[i] <- "dmp"}
}

data3<-data2[,c("animal_id",'birth_year', 'bs', "2007","2008","2009","2010","2011","2012","2013","2014","2015",
        "2016","2017","2018","2019","2020","2021", '2022')]


write.csv(data3, './cleaned/caphx2022.csv', row.names = F) 


#remove 0.5 year age class from analysis
  data4<- subset(data, age != 0.5)
  
  
  data5<-data4 %>% pivot_wider(names_from = cap_year, values_from = status, values_fill = 0, id_cols = animal_id)#
  
   #add birth year
  for (i in c(1:285)){
    x <- substr(data5$animal_id[i], 5, 8)
    data5$birth_year[i]<-x
    # if (x == "1") {data$bs[i] <- "dmp"}
  }
  
  for (i in 286:512){
    x <- substr(data5$animal_id[i], 4, 7)
    data5$birth_year[i]<-x
    # if (x == "1") {data$bs[i] <- "dmp"}
  }
  
  #add birth site
  
  for (i in 1:nrow(data5)) {
    x <- substr(data5$animal_id[i], 1, 1)
    if (x == "2") {data5$bs[i] <- "wy"} else 
    { data5$bs[i] <- "ey"}
  }
  ##make a new column with discernment of DMP and pasture born
  for (i in 1:nrow(data5)) {
    x <- substr(data5$animal_id[i], 10, 10)
    if (x == "1") {data5$bs[i] <- "dmp"}
  }
  data6<-data5[,c("animal_id",'birth_year','bs',"2008","2009","2010","2011","2012","2013","2014","2015",
                  "2016","2017","2018","2019","2020","2021",'2022')]
  
 
  
  
  write.csv(data6, './cleaned/caphx2022_nofawns.csv', row.names = F)
  
