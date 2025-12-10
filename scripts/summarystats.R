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
census <- read.csv('./cleaned/census.csv', header = T)

census$densitykm2 <- census$totaldeer/4.4 #4.4km2

densities <- census %>% group_by(pasture) %>% summarise(mean = mean(densitykm2),
                                                        sd = sd(densitykm2))

agecount$ageclass <- agecount$ageclass +.5 #add 0.5 to ages

#histogram for age distribution
hist <- ggplot(agecount, aes(ageclass, n)) +
  geom_col(fill = "black") +
  scale_x_continuous(breaks = unique(agecount$ageclass)) +
  geom_text(aes(label = n), vjust = -0.3, size = 6) +
  labs(x = "AGE CLASS", y = "NUMBER OF INDIVIDUALS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.1),          # x, y inside the plot area
        legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 20, hjust = 0.5),
        axis.text = element_text(face='bold',size = 20),
        #axis.text.x = element_text(margin = margin(t = 5)),
        panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
hist
ggsave('./figures/AGECOUNT.jpg', hist, width = 10, height = 10)


#summary stats for temperature
temp <- read.csv('./cleaned/temperaturedataNOAA.csv', header = T)

temp %>% group_by(month) %>% 
  summarise(mean = mean(maxtempC),      sd = sd(maxtempC))

#summary stats for rainfall

data %>% group_by(month)  %>% 
  summarise(mean = mean(cy.rain*2.54),
            sd = sd(cy.rain*2.54))


data <- read.csv('./raw/master_caphist.csv', header = T)
data <- data %>%  filter(Age > 0.5)
data %>%  group_by(Status) %>%  summarise(count = n())
