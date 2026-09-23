#code to clean master data
library(here)
library(tidyr)
library(dplyr)
library(stringr)

#taking raw capture history data
# data<- read.csv('./raw/master_caphist.csv', header = T)


data <- read.csv('./cleaned/fawncaphx.csv', header = T)

data <- data %>%
  mutate(
    birth_year = as.numeric(word(animal_id, 2, sep = fixed("-"))) # fill year with animal id
  )

data <- data %>%                                  # fill bs with animal ID
  mutate(
    part1 = word(animal_id, 1, sep = fixed("-")),
    part3 = word(animal_id, 3, sep = fixed("-")),
    
    bs = case_when(
      nchar(part1) == 2 ~ "ey",
      nchar(part1) == 3 & part3 == "0" ~ "wy",
      nchar(part1) == 3 & part3 == "1" ~ "dmp",
      TRUE ~ NA_character_
    )
  )

#QC check that all IDs follow correct format
data %>%
  filter(
    !str_detect(animal_id,
                "^([0-9]{2}|[0-9]{3})-[0-9]{4}-[01]-[0-9]+$")
  )
#check for any NA values or abnormal years
table(data$bs, useNA = "ifany")
table(data$birth_year, useNA = "ifany")
table(data$cap_year, useNA = "ifany")


#read in PMDI data
    library(lubridate)
    pmdi <- read.csv("./raw/pmdi.csv", header = F)
    pmdi <- pmdi[-c(1:2),]
    pmdi <- pmdi %>% rename(
                Date = V1,
                Value = V2
              )
    
    pmdi$Value <- as.numeric(pmdi$Value)
    
    pmdi <- pmdi %>%
      mutate(
        year = as.integer(substr(Date, 1, 4)),
        month = as.integer(substr(Date, 5, 6)),
        winter_year = ifelse(month == 12, year + 1, year) #biological winter, December 2015, January 2016, February 2016
      )
    
    winter_vals <- pmdi %>%   #
      filter(month %in% c(12, 1, 2)) %>%
      group_by(year = winter_year) %>%
      summarise(
        pmdi_winter = mean(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    seasonal_vals <- pmdi %>%
      group_by(year) %>%
      summarise(
        pmdi_annual = mean(Value, na.rm = TRUE),
        pmdi_spring = mean(Value[month %in% c(3,4,5)], na.rm = TRUE),
        pmdi_summer = mean(Value[month %in% c(6,7,8)], na.rm = TRUE),
        pmdi_fall = mean(Value[month %in% c(9,10,11)], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(winter_vals, by = "year")
    
    fawncaphx <- data %>%
      left_join(
        seasonal_vals,
        by = c("cap_year" = "year")
      )

#add densities
    library(dplyr)
    
    dens <- read.csv("./cleaned/adultdensities.csv") %>%
      mutate(
        bs = case_when(
          pasture == "East Yana" ~ "ey",
          pasture == "West Yana" ~ "wy"
        )
      ) %>%
      select(year, bs, densitykm2)
    
    # duplicate West Yana values for DMP
    dens_dmp <- dens %>%
      filter(bs == "wy") %>%
      mutate(bs = "dmp")
    
    dens <- bind_rows(dens, dens_dmp)

# join to capture history
    fawncaphx <- fawncaphx %>%
      left_join(
        dens,
        by = c("cap_year" = "year", "bs" = "bs")
      ) %>%
      rename(adult_density = densitykm2)
    
    table(is.na(fawncaphx$adult_density))
    
    fawncaphx %>%
      distinct(cap_year, bs, adult_density) %>%
      arrange(cap_year, bs)

#keep columns of interest
    names(fawncaphx)
    fawncaphx <- fawncaphx %>%
      select(
        animal_id,
        birth_year,
        cap_year,
        age,
        bs,
        status,
        adult_density,
        pmdi_annual,
        pmdi_spring,
        pmdi_summer,
        pmdi_fall,
        pmdi_winter
      )
    names(fawncaphx)
    unique(fawncaphx$status)

#update status to 0, 1, 2 based upon variable  
    fawncaphx <- fawncaphx %>%
      mutate(
        status_cam = case_when( #these include data obtained from cameras               
          status %in% c("Captured", "Alive-Cuddy") ~ 1,
          TRUE ~ 2
        ),
        
        status_nocam = case_when( #this exclude camera data
          status == "Captured" ~ 1,
          status == "Alive-Cuddy" ~ 0,
          TRUE ~ 2
        )
      )
    table(fawncaphx$status, fawncaphx$status_cam)
    
    table(fawncaphx$status, fawncaphx$status_nocam)
    
#create non-detection data
    study_end <- max(fawncaphx$cap_year)
    
    fawncaphx_full <- fawncaphx %>%
      group_by(animal_id) %>%
      complete(
        cap_year = seq(first(birth_year), study_end)
      ) %>%
      ungroup()

    fawncaphx_full <- fawncaphx_full %>%
      group_by(animal_id) %>%
      fill(
        birth_year,
        bs,
        .direction = "downup"
      ) %>%
      ungroup()  
    
    fawncaphx_full <- fawncaphx_full %>%
      mutate(
        status = ifelse(is.na(status), "not_detected", status),
        status_cam = ifelse(is.na(status_cam), 0, status_cam),
        status_nocam = ifelse(is.na(status_nocam), 0, status_nocam),
        age = cap_year - birth_year + 0.5
      )
    
# Remove old versions first to avoid duplicate columns
    fawncaphx_full <- fawncaphx_full %>%
      select(
        -adult_density,
        -pmdi_annual,
        -pmdi_spring,
        -pmdi_summer,
        -pmdi_fall,
        -pmdi_winter
      )
    
# Join annual density and PMDI information back by year
    fawncaphx_full <- fawncaphx_full %>%
      left_join(
        dens,
        by = c("cap_year" = "year", "bs" = "bs")
      ) %>%
      rename(adult_density = densitykm2) %>%
      left_join(
        seasonal_vals,
        by = c("cap_year" = "year")
      )

#QC
    fawncaphx_full %>%
      filter(animal_id == "270-2007-1-001") %>%
      arrange(cap_year)
    
#identify any individuals that have a detection in the same year for both alive-cuddy and captured
    fawncaphx_full %>%
      group_by(animal_id, cap_year) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1)
    
    fawncaphx_full %>%
      group_by(animal_id, cap_year) %>%
      filter(n() > 1) %>%
      arrange(animal_id, cap_year)

#only keep duplicate with "capture"
    fawncaphx_full <- fawncaphx_full %>%
      group_by(animal_id, cap_year) %>%
      filter(
        if(any(status == "Captured")) {
          status == "Captured"
        } else {
          TRUE
        }
      ) %>%
      slice(1) %>%
      ungroup()

#confirm there are no duplicates
    fawncaphx_full %>%
      count(animal_id, cap_year) %>%
      filter(n > 1)
    
#save file
    write.csv(fawncaphx_full, 'cleaned/fawncaphx.csv', row.names = F)
    
    
#how many detections with cam and without
    # Total detections
    cam_detections <- fawncaphx_full %>%
      filter(status_cam != 0) %>%
      nrow()
    
    # Unique animals detected
    cam_individuals <- fawncaphx_full %>%
      filter(status_cam != 0) %>%
      distinct(animal_id) %>%
      nrow()
    
    cam_detections
    cam_individuals
    
    # Total detections
    nocam_detections <- fawncaphx_full %>%
      filter(status_nocam != 0) %>%
      nrow()
    
    # Unique animals detected
    nocam_individuals <- fawncaphx_full %>%
      filter(status_nocam != 0) %>%
      distinct(animal_id) %>%
      nrow()
    
    nocam_detections
    nocam_individuals

#detection summary table
    summary_counts <- data.frame(
      dataset = c("CAM", "NOCAM"),
      detections = c(
        sum(fawncaphx_full$status_cam != 0),
        sum(fawncaphx_full$status_nocam != 0)
      ),
      unique_individuals = c(
        n_distinct(
          fawncaphx_full$animal_id[fawncaphx_full$status_cam != 0]
        ),
        n_distinct(
          fawncaphx_full$animal_id[fawncaphx_full$status_nocam != 0]
        )
      )
    )
    
    summary_counts
    