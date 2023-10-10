count<-data %>%  group_by(bs, birth_year) %>% summarise(n=n())

data2<- read.csv('./cleaned/caphx2022.csv', header = T)
count2<-data2 %>%  group_by(bs, birth_year) %>% summarise(n=n())
caphx<- read.csv('./cleaned/capture_log22.csv', header = T)
bucks<- read.csv('C:/Users/Joe/Documents/R Projects/Growth Curves/raw/bucks_nofawns.csv', header = T)
count3<-bucks %>%  group_by(birthsite, year_birth) %>% summarise(n=n())
