
#census cleaning
df <- read_csv("./raw/faithcensus.csv")

# Rename first column to something useful
df <- df %>%
  rename(category = `26-Jan-26`)

# Remove the first row (it's just repeated year labels)
df <- df %>%
  slice(-1)

# Pivot from wide to long
df_long <- df %>%
  pivot_longer(
    cols = -category,
    names_to = "year",
    values_to = "value"
  )

df_long$pasture <- NA
df_long$pasture[1:705] <- "East Yana"
df_long <- df_long %>% slice(-1411:-5190)
df_long$pasture[706:1410] <- "West Yana"
df_long <- df_long %>% slice(-1:-135,-151:-210, -256:-360, -376:-600, -676:-840, -856:-915, -961:-1065,
                             -1081:-1305, -1381:-1410)
df_long <- df_long[-c(1:135, 151:210,256:360,376:600,)]

#add 2022 data
add22ey<-
  tibble(
  category = c("TOTAL BUCKS", "DOES", "FAWNS", "TOTAL DEER"),  # etc
  year = 2022,
  value = c(37, 35, 18, 90),   # your 2022 numbers
  pasture = "East Yana")

add22wy<-
  tibble(
    category = c("TOTAL BUCKS", "DOES", "FAWNS", "TOTAL DEER"),  # etc
    year = 2022,
    value = c(52, 22, 6, 80),   # your 2022 numbers
    pasture = "West Yana")

df_long$year <- as.numeric(df_long$year)
df_long$value <- as.numeric(df_long$value)

df_long <- bind_rows(df_long, add22ey,add22wy)


write.csv(df_long, './cleaned/censusupto2025.csv', row.names = F)

library(dplyr)

total_adults <- df_long %>%
  filter(category %in% c("TOTAL BUCKS", "DOES")) %>%
  group_by(pasture, year) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(category = "total adults")

write.csv(total_adults, './cleaned/adultdensity.csv', row.names = F)
