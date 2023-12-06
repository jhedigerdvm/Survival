
# write.csv(gather, './output/posterior_long.csv', row.names = F)
# posterior<- tidy_draws(cjs.rain.site.age)
# posterior<- posterior[,-c(1:26)]
# posterior <- posterior[,-33001]
# posterior<- posterior[,1:3000]
# 
# #create dataframe with posteriors of just survival age1 across the three sites
# #pivot longer puts them in a tibble format
# posterior_long <- posterior %>% pivot_longer(everything())

# #string search
# posterior_long %>%
#   as_tibble() %>%
#   mutate(rain.sim = rep(rain.sim, 9000),
#          site = ifelse(str_detect("\\d+\\,1\\,1\\]$"),"one",
#                        ifelse(str_detect("\\d+\\,2\\,1\\]$"),"two","three")),
#          age = ifelse(str_detect("1\\]$"),"one",
#                       ifelse(str_detect("2\\]$"),"two","three"))
#                               
# str_detect(pattern = "\\d+\\,1\\,1\\]$", posterior_long$name[1:100])

#make column for rainfall data
# low<- cjs.rain.site.age$q2.5$survival
# 
# mean1 <- cjs.rain.site.age$mean$survival %>% as_tibble %>% pivot_longer(everything(), cols_vary = 'slowest')
# low<- cjs.rain.site.age$q2.5$survival %>% as_tibble() %>%  pivot_longer(everything())
# high<- cjs.rain.site.age$q97.5$survival %>% as_tibble() %>%  pivot_longer(everything())
# 
# posterior2 <- cbind(mean, low, high)
# 
# 
# posterior_long$rain <- rep(rain.sim, nrow(posterior_long)/1000)
# posterior_long$p2.5 <- rep(p2.5$p2.5, nrow(posterior_long)/3000)
# posterior_long$p97.5 <- rep(p97.5$p97.5, nrow(posterior_long)/3000)

# # Set total rows
# num_rows <- 9000000
# 
# # Create empty values vector
# site <- vector(length = num_rows)
# 
# # Index for filling
# index <- 1
# 
# # Repeat till last row
# while(index <= num_rows){
# 
#   # Insert 1
#   site[index:(index+999)] <- 1
#   index <- index + 1000
# 
#   # Insert 2
#   site[index:(index+999)] <- 2
#   index <- index + 1000
# 
#   # Insert 3
#   site[index:(index+999)] <- 3
#   index <- index + 1000
# 
# }
# 
# 
# posterior_long<- cbind(posterior_long,site)
# 
# posterior_long$site<- as.factor(posterior_long$site)

# data<- read.csv('./output/posterior_long.csv', header = T)