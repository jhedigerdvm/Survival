#code to clean master data
library(here)
library(tidyr)
library(dplyr)

#taking raw capture history data
data<- read.csv('./raw/master_caphist.csv', header = T)

#rename column names
names(data)<- c('animal_id', 'birth_year', 'age', 'status', 'cap_year')

# replace status with numbers
data$status[data$status == 'Captured'] <- 1
data$status[data$status == 'Nat Mortality'] <- 0 #found dead
data$status[data$status == 'Alive-Cuddy'] <- 1  #not capture physically but seen on camera
data$status[data$status == 'Harvested'] <- 2    #accidentally harvested
data$status[data$status == 'Dead-Capture'] <- 2 #died during capture event

#check for duplicated entries
data$check<-duplicated(data) #check for duplicates
#270-2010-1-019 #270-2011-1-215 #270-2019-1-909 #90-2020-0-2080
condition<- data$check=='TRUE'
data1<- data[!condition,]

#some birthyears are wrong, code to assign birthyear based upon animal id
#first we need to reset row names
rownames(data1) = seq(length=nrow(data1))

for (i in c(1:2034, 3609:3615)){
  x <- substr(data1$animal_id[i], 5, 8)
  data1$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}

for (i in 2035:3608){
  x <- substr(data1$animal_id[i], 4, 7)
  data1$birth_year[i]<-x
  # if (x == "1") {data$bs[i] <- "dmp"}
}

#ensure capture year is correct for all
data1$birth_year<-as.numeric(data1$birth_year)
data2<-data1

#make age whole numbers, 0.5 -> 0, 1.5 -> 1
data2$age1 <- data1$age - 0.5

#to confirm capture year is correct, add age to birth year and make new col capyear1
for(i in 1:nrow(data2)){
  x <- data2$birth_year[i] + data2$age1[i]
  data2$cap_year1[i]<-x
}

#check to make sure that cap_year1 and capyear are the same
check<-data2$cap_year - data2$cap_year1
any(check>0) #they are

#remove unnecessary columns
data3<-data2[,-c(3,8,6)]
names(data3)[names(data3) == 'age1'] <- 'age'

#add birth-site to data
for (i in 1:nrow(data3)) {
  x <- substr(data3$animal_id[i], 1, 1)
  if (x == "2") {data3$bs[i] <- "wy"} else 
  { data3$bs[i] <- "ey"}
}
##make a new column with discernment of DMP and pasture born
for (i in 1:nrow(data3)) {
  x <- substr(data3$animal_id[i], 10, 10)
  if (x == "1") {data3$bs[i] <- "dmp"}
}

write.csv(data3,'./cleaned/capture_cleaned.csv', row.names = F)

#remove 0.5 year age class from analysis
data4<- subset(data3, age != 0)

write.csv(data4,'./cleaned/capture_cleaned_nofawns.csv',row.names = F)

#prepare to put data4 into wide format
data5<-data4[,-c(2,5,6)]
data5$cap_year<-as.character(data5$cap_year)
data6<-data4 %>% pivot_wider(names_from = cap_year, values_from = status, values_fill = '0', id_cols = animal_id)#


#some birthyears are wrong, code to assign birthyear based upon animal id
for (i in c(1:285)){
  x <- substr(data6$animal_id[i], 5, 8)
  data6$birth_year[i]<-x
}

for (i in 286:512){
  x <- substr(data6$animal_id[i], 4, 7)
  data6$birth_year[i]<-x
}

#add birth-site to data
for (i in 1:nrow(data6)) {
  x <- substr(data6$animal_id[i], 1, 1)
  if (x == "2") {data6$bs[i] <- "wy"} else 
  { data6$bs[i] <- "ey"}
}
##make a new column with discernment of DMP and pasture born
for (i in 1:nrow(data6)) {
  x <- substr(data6$animal_id[i], 10, 10)
  if (x == "1") {data6$bs[i] <- "dmp"}
}

#reorder columns to move from left to right sequentially in time
data7<- data6[,c(1,17,18,2,3,4,5,6,7,8,10,9,11,13,12,14,15,16)]

data8<-data7[,-c(1:3)]
data8$check<-apply(data8, 1, function(r) any(r %in% c("2")))
# 6 individuals first and last capture occasion were found
# 131 141 151 169 204 483
data9<- data7[-c(131,141,151,169,204,483),]#removing individuals whose first and last capture occasions were the same

write.csv(data9,'./cleaned/capture_cleaned_nofawns_wide.csv', row.names = F)

#take wide format and put back into long format
data_wide<- read.csv('./cleaned/capture_cleaned_nofawns_wide.csv', header = T)
data_long <- data_wide %>%  pivot_longer(cols = -c(animal_id, birth_year, bs),
                                         names_to = 'year', values_to = 'status')

#remove prefix x from 'year'
data_long<- data_long %>% mutate(year = substr(year,2,5))
data_long$year <-as.integer(data_long$year)

#add age class 
data_long$ageclass<- 0

  for (i in 1:nrow(data_long)){
    data_long$ageclass[i]<- data_long$year[i] - data_long$birth_year[i]
  }

data_long$ageclass[data_long$ageclass<=0] <- NA #remove any ages that are less than 0

#add rainfall to data long
rainfall<- read.csv('./cleaned/rainfall_clean_nov_oct1.csv', header = T)
dmp <- rainfall[rainfall$site=='wy',] #need to add site 'dmp' to rainfall, same as wy
dmp$site<- 'dmp' #rename wy to dmp 
rainfall1<- rbind(rainfall,dmp)
write.csv(rainfall1, './cleaned/rainfall.csv')
data_long1<- left_join(data_long,rainfall1, by = c('year' = 'Year',
                                                  'bs' = 'site' )) 

data<- data_long1 
data<- data[data$birth_year != '2021',] #remove 2021 cohort because first capture is also last capture

#remove individuals with no capture occasions
data1<-
  data %>% 
  group_by(animal_id) %>%
  summarise(status_sum = sum(status)) %>%
  filter(status_sum < 1) -> to_remove

data2 <- data[!data$animal_id %in% to_remove$animal_id,]

data2$annual.sc <-scale(data2$annual) #scale and center data
data2$sum.march.apr.may.sc<- scale(data2$sum.march.apr.may)
data2$sum.jun.jul.aug.sc<- scale(data2$sum.jun.jul.aug)

write.csv(data2, './cleaned/caphx.rainfall.nov.oct1.csv', row.names = F)
                    

###########
#try to merge bcs and weight into caphx 
ch.rain <- read.csv('./cleaned/caphx.rainfall.nov.oct1.csv', header = T)
bucks <- read.csv("C:/Users/Joe/Documents/3-R Projects/Growth Curves/clean/nofawns22.csv")

bucks$age <- bucks$year_cap - bucks$year_birth 
unique(bucks$year_birth)

#remove 2021 born fawns because first capture is also last capture
bucks<-bucks[bucks$year_birth != '2021',] #remove 2021 cohort because first capture is also last capture

#rename column year cap to year to match capture history format
bucks <- bucks %>% rename(year = year_cap)
bucks <- bucks %>% rename(birth_year = year_birth)

bucks <- bucks %>% rename(bs = birthsite)
bucks <- bucks %>% rename(ageclass = age)

#add weights and bcs to capture history 
join<- ch.rain %>%  left_join(bucks)

join$weightkg <- join$weight/2.2
join$bcscm <- join$bcsin * 2.54

write.csv(join, './cleaned/final.ch.csv', row.names = F)
