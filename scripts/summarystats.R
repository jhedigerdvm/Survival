#summary stats for survival
#average pmdi and SE

library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

data <- read.csv('./cleaned/ch.pmdi.csv', header = T)

pmdistats <- data %>% summarise(mean = mean(pmdi_spring),
                                sd = sd(pmdi_spring))
agecount <- data %>% filter(status == 1) %>%  distinct(animal_id, ageclass) %>%   # keep one row per individual per ageclass
  count(ageclass)

pasture <- c('control', "control", "control","control","control","control","control","control","control",
             "control","control","control","control","control","control","control",
             'treatment', "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
             "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
             "treatment", "treatment")
year <- c('2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
          '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021',
          '2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
          '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021')
totaldeer<- c(137, 156,	135,	165	,136,	200	,189,	165,	238,	268	,124,	82,	79,	118,	117,	124,
             #west yana 
             132, 77,	59,	91,	82,	173,	180,	208	,241,	287	,143,	118,	108,	135,	206,	147)

census <- data.frame(pasture, year, totaldeer)
write.csv(census, './cleaned/census.csv', row.names = F)

census$densitykm2 <- census$totaldeer/4.4 #4.4km2

densities <- census %>% group_by(pasture) %>% summarise(mean = mean(densitykm2),
                                                        sd = sd(densitykm2))
