
#create birth_sites with for loop
for (i in 1:nrow(data)) {
  x <- substr(data$animal_id[i], 1, 1)
  if (x == "2") {data$wyey[i] <- "west yana"} else 
  { data$wyey[i] <- "east yana"}
}
##make a new column with discernment of wy and pasture born
for (i in 1:nrow(data)) {
  x <- substr(data$animal_id[i], 10, 10)
  if (x == "1") {data$dmp_pasture[i] <- "dmp"} else 
  { data$dmp_pasture[i] <- "pasture"}
}

##1 vector with 3 class variable 
data$birth_site <- paste(data$wyey, data$dmp_pasture, sep = "_")

data$birth_site<-as.factor(data$birth_site)
data$animal_id<-as.factor(data$animal_id)

##make new column with birth year
#create birth_sites with for loop
for (i in 1:110) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2007") {data$birth_year[i] <- "2007"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 111:205) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2008") {data$birth_year[i] <- "2008"}
}
##make new column with birth year
#create birth_sites with for loop
for (i in 206:302) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2009") {data$birth_year[i] <- "2009"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 303:487) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2010") {data$birth_year[i] <- "2010"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 488:733) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2011") {data$birth_year[i] <- "2011"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 734:989) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2012") {data$birth_year[i] <- "2012"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 989:1188) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2013") {data$birth_year[i] <- "2013"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1189:1363) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2014") {data$birth_year[i] <- "2014"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1364:1517) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2015") {data$birth_year[i] <- "2015"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1518:1661) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2016") {data$birth_year[i] <- "2016"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1662:1727) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2017") {data$birth_year[i] <- "2017"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1728:1791) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2018") {data$birth_year[i] <- "2018"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1792:1847) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2019") {data$birth_year[i] <- "2019"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 1848:1881) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2020") {data$birth_year[i] <- "2020"}
}

for (i in 1882:1902) {
  x <- substr(data$animal_id[i], 5, 8)
  if (x == "2021") {data$birth_year[i] <- "2021"}
}


##make new column with birth year
#create birth_sites with for loop
for (i in 1903:2026) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2007") {data$birth_year[i] <- "2007"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2027:2163) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2008") {data$birth_year[i] <- "2008"}
}
##make new column with birth year
#create birth_sites with for loop
for (i in 2164:2337) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2009") {data$birth_year[i] <- "2009"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2338:2472) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2010") {data$birth_year[i] <- "2010"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2473:2606) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2011") {data$birth_year[i] <- "2011"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2607:2731) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2012") {data$birth_year[i] <- "2012"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2732:2864) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2013") {data$birth_year[i] <- "2013"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2865:2990) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2014") {data$birth_year[i] <- "2014"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 2991:3104) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2015") {data$birth_year[i] <- "2015"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 3105:3179) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2016") {data$birth_year[i] <- "2016"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 3180:3238) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2017") {data$birth_year[i] <- "2017"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 3239:3281) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2018") {data$birth_year[i] <- "2018"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 3282:3319) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2019") {data$birth_year[i] <- "2019"}
}

##make new column with birth year
#create birth_sites with for loop
for (i in 3320:3349) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2020") {data$birth_year[i] <- "2020"}
}

for (i in 3350:3364) {
  x <- substr(data$animal_id[i], 4, 7)
  if (x == "2021") {data$birth_year[i] <- "2021"}
}

#creat CSV from data frame 
write.csv(data, '.\\cap_hx.csv',row.names = FALSE)


