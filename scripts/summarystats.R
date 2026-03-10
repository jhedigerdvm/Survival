#summary stats for survival
#average pmdi and SE

library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(mcmcr) 
library(viridis)
library(here)

# data <- read.csv('./cleaned/ch.pmdi.csv', header = T) #data that does not include fawns, used for surv analysis
data <- read.csv('./cleaned/fawncaphx.csv', header = T) #dataset that includes fawns

#cleanup birthyear for 2022 cohort
# unique(data$birth_year)
# data1 <- data %>%
#   mutate(birth_year_from_id = as.integer(str_split(animal_id, "-", simplify = TRUE)[,2]))
# data1 <- data1 %>%
#   mutate(birth_year = birth_year_from_id)
# 
# data <- data %>% mutate(across(c('animal_id', 'birth_year','age', 'status', 'cap_year'), factor))
# 
# #bin dead, harvest, nat mortality into known death
# data1$status1 <- 
#   recode(data1$status, "Alive-Cuddy" = "Captured",
#          "Dead-Capture" = "Known Death",
#          "Harvested" = "Known Death",
#          "Nat Mortality" = "Known Death")
# 
# 
# 
# write.csv(data1, './cleaned/fawncaphx.csv', row.names = F)

unique(data$status)

# pmdistats <- data %>% summarise(mean = mean(pmdi_spring),
#                                 sd = sd(pmdi_spring))



# ---- Counts by ageclass and plots----



#create a datafram with all individuals entering the study, this includes fawns and older who were captured physically or on camera
agecount <- data %>% filter(status1%in% c("Captured")) %>%  
  distinct(animal_id, age) %>%   
  count(age) 

install.packages("ggpattern")
library(ggpattern)

agecount$my_pattern <- ifelse(agecount$age == "0.5", "stripe", "none")
# agecount$ageclass <- agecount$ageclass +.5 #add 0.5 to ages

agecount <- agecount %>% mutate(across(c(age, my_pattern),factor))

#histogram for age distribution
hist <- ggplot(agecount, aes(age, n, fill = my_pattern)) +
  geom_col(color = "black") +
  # scale_x_continuous(breaks = unique(agecount$age)) +
  scale_fill_manual(values = c("stripe" = "grey70", "none" = "white"))+
  geom_text(aes(label = n), vjust = -0.3, size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  labs(x = "AGE CLASS", y = "NUMBER OF INDIVIDUALS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        # legend.position.inside = c(0.1,0.1),          # x, y inside the plot area
        # legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        # legend.text = element_text(size = 28),
        # legend.title = element_blank(),
        # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 20, hjust = 0.5),
        axis.text = element_text(face='bold',size = 20),
        #axis.text.x = element_text(margin = margin(t = 5)),
        panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
hist
ggsave('./figures/AGECOUNT.jpg', hist, width = 10, height = 10)

#count how many individual occasions an individual was captured
# data1 <- data %>% filter(status%in% "Captured") #filter by individuals captured each year
capcounts<-
  data %>%filter(status%in% c("Captured", "Alive-Cuddy")) %>%  
  distinct(animal_id, cap_year, .keep_all = TRUE) %>%  # one row per animal per year
  count(animal_id) %>%        # total captures per individual
  #filter(n %in% c(13, 14))

  count(n, name = "n_individuals") %>% # how many individuals fall into each total
  rename(n_occasions = n)

hist <- ggplot(capcounts, aes(n_occasions, n_individuals)) +
  geom_col(fill = "white", color = "black") +
  scale_x_continuous(breaks = unique(capcounts$n_occasions)) +
  geom_text(aes(label = n_individuals), vjust = -0.3, size = 6) + #adds count on top of bars
  labs(x = "NUMBER OF CAPTURES PER INDIVIDUAL", y = "NUMBER OF INDIVIDUALS") +
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
ggsave('./figures/capcountwfawns.jpg', hist, width = 10, height = 10)


# ---- CAPTURE PER SITE PER YEAR ----

#captures per site per year
# data1 <- data %>% filter(status == "1") #filter by individuals captured each year
# data1 <- data1 %>%
#   mutate(bs = recode(bs,
#                      "dmp" = "wy"))
cap.site.year<-
  data %>%
  group_by(bs,cap_year) %>%        
  summarise(n_captures = n(), .groups = "drop")

cap.site.year$bs <- as.factor(cap.site.year$bs)
cap.site.year <- cap.site.year %>% mutate(bs = recode(bs, "wy" = "West Yana", "ey" = "East Yana"))

hist <- ggplot(cap.site.year, aes(cap_year, n_captures, fill = bs)) +
  geom_col(position = position_dodge(width=.9)) +
  scale_x_continuous(breaks = unique(cap.site.year$cap_year)) + #tick for each year
  geom_text(position = position_dodge(width = 0.9), aes(label = n_captures), vjust = -0.3, size = 4) + #adds count on top of bars
  labs(x = "CAPTURE YEAR", y = "NUMBER OF INDIVIDUALS") +
  scale_fill_grey(start = 0.4, end = 0.8, name = "bs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.8),          # x, y inside the plot area
        legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        
            
        # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 14, hjust = 0.5),
        axis.text = element_text(face='bold',size = 14),
        #axis.text.x = element_text(margin = margin(t = 5)),
        panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
        plot.background = element_rect(fill='transparent', color=NA)
         # light/dark shades
  ) #transparent plot bg)
hist
ggsave('./figures/capcountbyyearwfawns.jpg', hist, width = 10, height = 5)



# ---- CAPTURE HISTORY PLOTS ----
data <- data %>% mutate(across(c('animal_id', 'birth_year','age', 'status', 'cap_year'), factor))


status<- data %>% group_by(status) %>% 
            summarise(n = n())

#bin dead, harvest, nat mortality into known death
data1$status1 <- 
          recode(data1$status, "Alive-Cuddy" = "Captured",
                              "Dead-Capture" = "Known Death",
                              "Harvested" = "Known Death",
                              "Nat Mortality" = "Known Death")


plot_caphx<- 
ggplot(data1, aes(
  x = cap_year,
  y = fct_reorder(animal_id, birth_year, .desc = TRUE),
  group = animal_id,
  shape = status1,
  color = status1
)) +
  geom_point(size = 1) +
  labs(
    x = "Capture Year",
    y = "Individual Capture History"
  ) +
  scale_shape_manual(
    values = c(16, 17),
    guide = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values = c("grey20", "darkred"))+ # 
  scale_y_discrete(expand = expansion(add = 10))+
  geom_line(color = "grey80", alpha = 0.3) +
  # scale_color_viridis_d(option = "B") +
  
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.position = "inside",
      legend.position.inside = c(0.1,0.2),          # x, y inside the plot area
      legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
      legend.text = element_text(size = 32),
      legend.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
      axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
      axis.text = element_text(face='bold',size = 24),
      #axis.text.x = element_text(margin = margin(t = 5)),
      panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
      plot.background = element_rect(fill='transparent', color=NA)
      # light/dark shades
) #transparent plot bg)

plot_caphx
# ggsave('./figures/phi.AGE.1site.jpg', phi.plot, width = 10, height = 10)
# ggsave('./figures/phi.AGE.1site.jpg', phi.plot, width = 10, height = 10)
ggsave('./figures/caphx.plot.jpg', plot_caphx, width = 15, height = 20)

# ---- CENSUS DATA PLOTS, DENSITY ----

# 
# pasture <- c('control', "control", "control","control","control","control","control","control","control",
#              "control","control","control","control","control","control","control",
#              'treatment', "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
#              "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
#              "treatment", "treatment")
# year <- c('2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
#           '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021',
#           '2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
#           '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021')
# totaldeer<- c(137, 156,	135,	165	,136,	200	,189,	165,	238,	268	,124,	82,	79,	118,	117,	124,
#              #west yana 
#              132, 77,	59,	91,	82,	173,	180,	208	,241,	287	,143,	118,	108,	135,	206,	147)
# 
# census <- data.frame(pasture, year, totaldeer)
# write.csv(census, './cleaned/census.csv', row.names = F)
# census <- read.csv('./cleaned/census.csv', header = T)
# 
# census$densitykm2 <- census$totaldeer/4.4 #4.4km2
# 
# densities <- census %>% group_by(pasture) %>% summarise(mean = mean(densitykm2),
#                                                         sd = sd(densitykm2))


#create densities of adults deer, not fawns
# pasture <- c('control', "control", "control","control","control","control","control","control","control",
#              "control","control","control","control","control","control","control",
#              'treatment', "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
#              "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", 
#              "treatment", "treatment")
# year <- c('2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
#           '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021',
#           '2006', '2007',	'2008',	'2009',	'2010',	'2011',	'2012',	'2013',	'2014',
#           '2015',	'2016',	'2017',	'2018',	'2019',	'2020',	'2021')
# totaldeer<- c(137, 156,	135,	165	,136,	200	,189,	165,	238,	268	,124,	82,	79,	118,	117,	124,
#              #west yana 
#              132, 77,	59,	91,	82,	173,	180,	208	,241,	287	,143,	118,	108,	135,	206,	147)
# 
# census <- data.frame(pasture, year, totaldeer)
#plots for census data, total deer, bucks, does and fawns
census <- read.csv('./cleaned/census.csv', header = T)


census <- census %>% mutate(pasture = recode(pasture, "treatment" = "West Yana", "control" = "East Yana"))


# ---Density stats --- ##

dens <- read.csv('./cleaned/adultdensities.csv', header = T)

dens %>% group_by(pasture) %>% 
  summarise(mean = mean(densitykm2),
            sd = sd(densitykm2))


dens <- read.csv('./cleaned/adultdensities.csv', header = T)

dens %>% group_by(pasture) %>% 
  summarise(mean = mean(densitykm2),
            sd = sd(densitykm2))

hist <- ggplot(census, aes(year, totaldeer, fill = pasture)) +
  geom_col(position = position_dodge(width=.9)) +
  scale_x_continuous(breaks = unique(census$year)) + #tick for each year
  geom_text(position = position_dodge(width = 0.9), aes(label = totaldeer), vjust = -0.3, size = 4) + #adds count on top of bars
  labs(x = "YEAR", y = "TOTAL NUMBER OF DEER") +
  scale_fill_grey(start = 0.4, end = 0.8, name = "pasture") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.8),          # x, y inside the plot area
        legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        
        
        # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 14, hjust = 0.5),
        axis.text = element_text(face='bold',size = 14),
        #axis.text.x = element_text(margin = margin(t = 5)),
        panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
        plot.background = element_rect(fill='transparent', color=NA)
        # light/dark shades
  ) #transparent plot bg)
hist
ggsave('./figures/censusdata.jpg', hist, width = 10, height = 5)


#plots and calculations for adult density
census <- read.csv('./cleaned/adultdensity.csv', header = T)

census$densitykm2 <- census$value/4.4 #4.4km2
census$densitykm2<- round(census$densitykm2, 2)

densities <- census %>% group_by(pasture) %>% summarise(mean = mean(densitykm2),
                                                        sd = sd(densitykm2))
census$year <- as.factor(census$year)

write.csv(census, './cleaned/adultdensities.csv', row.names = F)

hist <- ggplot(census, aes(year, densitykm2, fill = pasture)) +
  geom_col(position = position_dodge(width=.9)) +
  # scale_x_continuous(breaks = unique(census$year)) + #tick for each year
  # geom_text(position = position_dodge(width = 0.9), aes(label = densitykm2), vjust = -0.3, size = 4) + #adds count on top of bars
  labs(x = "YEAR", y = "ADULT DEER DENSITY (KM2)") +
  scale_fill_grey(start = 0.4, end = 0.8, name = "pasture") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.8),          # x, y inside the plot area
        legend.justification = c("left", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        # plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 14, hjust = 0.5),
        axis.text = element_text(face='bold',size = 14),
        #axis.text.x = element_text(margin = margin(t = 5)),
        panel.background = element_rect(fill='transparent'), #transparenhttp://127.0.0.1:46083/graphics/815b1ae8-dcf1-4f7c-921f-7bb4b3b81021.pngt panel bg
        plot.background = element_rect(fill='transparent', color=NA)
        # light/dark shades
  ) #transparent plot bg)
hist
ggsave('./figures/densitiesadults.jpg', hist, width = 10, height = 5)


# ---- TEMPERATURE STATS FROM NOAA ----

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

