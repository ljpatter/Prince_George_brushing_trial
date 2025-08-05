# ---
# title: "2_Clean data"
# author: "Leonard Patterson"
# created: "2025-12-12"
# description: "Code to clean data brought in from wildtrax"
# ---

# Clear environment
rm(list = ls())  # Removes all objects from the environment

##Load packages----
library(tidyverse)
library(wildrtrax)
library(xfun)
library(expss)
library(lubridate)

##////////////////////////////////////////////////////////////////
# Login to WildTrax----
# NOTE: Edit the 'loginexample.R' script to include your WildTrax 
# login details and rename to 'login.R'. 
# DO NOT PUSH YOUR LOGIN TO GITHUB
config <- "Scripts/login.R"
source(config)
wt_auth()

## Import data----
#load("Output/R Data/PGBT_raw_aru_2025-07-03.rData")
#load("Output/R Data/PGBT_raw_aru_2025-07-03.rData")

## Interim load
PGBT_2023 <- read.csv("Input/Tag reports/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2023_tag_report.csv")
PGBT_2024 <- read.csv("Input/Tag reports/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2024_tag_report.csv")
PGBT_2025 <- read.csv("Input/Tag reports/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2025_Limited_Amplitude_Processing_tag_report.csv")



# Filter to only LP as observer
PGBT_2023 <- PGBT_2023 %>%
  filter(observer =="Leonard Patterson") 
PGBT_2024 <- PGBT_2024 %>%
  filter(observer =="Leonard Patterson") 
PGBT_2025 <- PGBT_2025 %>%
  filter(observer =="Leonard Patterson") 

# Abundance column was character for two reports; change to integer
PGBT_2023$abundance <- as.integer(PGBT_2023$abundance)
PGBT_2024$abundance <- as.integer(PGBT_2024$abundance)
PGBT_2025$abundance <- as.integer(PGBT_2025$abundance)

# Bind 2023 and 2024 data into a single dataframe
dat1 <- bind_rows(PGBT_2023, PGBT_2024, PGBT_2025, id = NULL)


### Remove water-associated species from tag report

# Get AVES species, excluding water/shore/game birds, etc.
aves <- wt_get_species() %>%
  filter(
    species_class == "AVES",
    !species_order %in% c("Podicipediformes","Pelecaniformes","Anseriformes",
                          "Gaviiformes","Gruiformes","Galliformes",
                          "Piciformes","Charadriiformes")
  )

# Define specific species codes to drop
drop_codes <- c("NOWA","TRES","COYE","UNTR","UNPA","MODO",
                "BEKI","UNCV","UNFL","UNWA","SWSP", "UNBI")

# Remove unwanted species
dat1_filtered <- dat1 %>%
  # Keep only target AVES (by species_code)
  semi_join(aves, by = "species_code") %>%
  # Drop specific unwanted codes
  filter(!species_code %in% drop_codes)

# Save
write.csv(dat1_filtered, "Input/Tag reports/tags_all_years.csv")










## ============================================================
## MEAN COUNT — UNLIMITED DISTANCE (UD) 
## ============================================================

# Load data
dat1_UD <- read.csv("Input/Tag reports/tags_all_years.csv")

# Ensure individual_order is numeric and extract year

dat2_UD <- dat1_UD %>%
  mutate(
    individual_order = as.numeric(individual_order),
    year = year(ymd_hms(recording_date_time))
  ) %>%
  
  # Aggregate: maximum individual_order per site-year-species
  group_by(location, year, recording_date_time, species_code) %>%
  summarise(
    max_count = max(individual_order, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Wide format: one row per site-year, one column per species
  pivot_wider(
    id_cols    = c(location, year, recording_date_time),
    names_from = species_code,
    values_from = max_count,
    values_fill = 0
  ) %>%
  arrange(location, year, recording_date_time)

# Save WITH SINGLETONS
write.csv(dat2_UD, "Output/Tabular Data/mean_count_all_years_UD_with_singletons.csv", row.names = FALSE)

# Remove singletons (species occurring at <= 3 site-years)
dat3_UD <- dat2_UD[ , sapply(dat2_UD, function(x) !is.numeric(x) || sum(x != 0, na.rm = TRUE) > 3) ]

# Save WITHOUT SINGLETONS
write.csv(dat3_UD, "Output/Tabular Data/mean_count_all_years_UD.csv", row.names = FALSE)















## ============================================================
## MEAN COUNT — LIMITED AMPLITUDE (LA) 
## ============================================================

# Load data
#dat_LA <- read.csv("Output/Tabular Data/truncated_count_150m_with_distances.csv")
dat_LA <- read.csv("Output/Tabular Data/truncated_count_150m.csv")

# Create year column
dat2_LA <- dat_LA %>%
  mutate(
    year = year(ymd_hms(recording_date_time))
  )

# Save UD mean count table
write.csv(dat2_LA,
          "Output/Tabular Data/mean_count_all_years_LA_with_singletons.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# 3. Remove rare species (e.g., species present at <= 3 site-years)
# ------------------------------------------------------------

# Remove singletons (species occurring at <= 3 site-years)
dat3_LA <- dat2_LA[ , sapply(dat2_LA, function(x) !is.numeric(x) || sum(x != 0, na.rm = TRUE) > 3) ]

# Save WITHOUT SINGLETONS
write.csv(dat3_LA,
          "Output/Tabular Data/mean_count_all_years_LA.csv",
          row.names = FALSE)
