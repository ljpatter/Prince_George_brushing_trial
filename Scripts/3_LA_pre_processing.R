# ---
# title: "LA Pre-processing"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: "Clean and wrangle tag data before limited amplitude processing"
# ---

# Clear environment
rm(list = ls())  # Removes all objects from the environment

# Load packages
library(tidyverse)

# Load tag data
dat1 <- read.csv("Input/Tag reports/tags_all_years.csv")

### Tag report does not include weather/malfunction notes - BOO. 
### Load in task report and join comments with tag data

task_2023 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2023_Tasks_2025-12-08.csv")
task_2024 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2024_Tasks_2025-12-08.csv")
task_2025 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2025 Limited Amplitude Processing_Tasks_2025-12-08.csv")

# Filter to just LP tasks for 2023 and 2024 (only LP in 2025)
task_2023_filt <- task_2023  %>%
  filter(observer == "Leonard Patterson")

task_2024_filt <- task_2024  %>%
  filter(observer == "Leonard Patterson")

task_2025_filt <- task_2025  %>%
  filter(observer == "Leonard Patterson")

# Combined all recordings
all_recordings <- bind_rows(task_2023_filt, task_2024_filt, task_2025_filt)

### For some reason, there are recordings that were not transcribed that are coming up as completed
### Get completed sites from tag data

# Get unique survey times from tag data
dat2 <- dat1 %>%
  distinct(location, recording_date_time)

# Trim all recordings so can join with surveys
dat2_clean <- dat2 %>%
  mutate(
    location            = str_trim(as.character(location)),
    recording_date_time = str_trim(as.character(recording_date_time))
  )

all_recordings_clean <- all_recordings %>%
  mutate(
    location            = str_trim(as.character(location)),
    recording_date_time = str_trim(as.character(recording_date_time))
  )

# Remove recordings that are not present in tag data
all_recordings_filtered <- all_recordings_clean %>%
  semi_join(dat2_clean, by = c("location", "recording_date_time"))

# Save cleaned recordings
write.csv(all_recordings_filtered, "Output/Tabular Data/LA inputs/PGBT_recordings_all_years.csv")






### A hand full of tags did not extract amplitude values; code below extracts
### tags WITH amplitude values; had to transcribe individual tags in a new project

# Load new tag data
fixed_2025 <- read.csv("Input/Fixed tags/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2024_tag_report.csv")
fixed_2023_2024 <- read.csv("Input/Fixed tags/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2025_Limited_Amplitude_Processing_tag_report.csv")

# Filter to only new 2025 tags
fixed_2025_1 <- fixed_2025 %>%
  mutate(year = year(recording_date_time)) %>% 
  filter(year == "2025")

# Filter to only new 2023/2024 tags
fixed_2023_2024_1 <- fixed_2023_2024 %>%
  mutate(year = year(recording_date_time)) %>% 
  filter(year == "2023" | year == "2024")

## Merge fixed 2023, 2024, and 2025 tags
# Abundance is a character in one on the dataframes for some reason. Change to integer.
fixed_2025_1$abundance <- as.integer(fixed_2025_1$abundance)

# Now merge fixed 2024 and 2025 dfs and remove year
fixed_all_tags <- bind_rows(fixed_2023_2024_1, fixed_2025_1)

# Remove faulty tags from main tag df
main_tag_df <- dat1 %>%
  filter(!is.na(right_full_freq_tag_rms_peak_dbfs)) %>%
  select(-X)

# Fix column structure so can merge w/ main tag df
fixed_all_tags <- fixed_all_tags %>%
  mutate(
    task_duration = as.integer(task_duration),
    min_tag_freq  = as.integer(round(as.numeric(gsub("kHz", "", min_tag_freq)) * 1000)),
    max_tag_freq  = as.integer(round(as.numeric(gsub("kHz", "", max_tag_freq)) * 1000))
  ) %>%
  rename(task_is_complete = is_complete) %>%
  select(-year)

# Finally, merge main df with new fixed tags
main_tag_df_2 <- bind_rows(main_tag_df, fixed_all_tags)



########## For recordings with only one good mic channel, replace
########## with amplitude values from the working mic

# Subset task reports to create comment df to just to tag report
comment_2023 <- task_2023  %>%
  select(location, recording_date_time, task_comments)
comment_2024 <- task_2024  %>%
  select(location, recording_date_time, task_comments)
comment_2025 <- task_2025  %>%
  select(location, recording_date_time, task_comments)

# Merge into one df
comment_combined <- bind_rows(comment_2023, comment_2024, comment_2025)

# Filter to only include sites with relevant comments
comment_combined <- comment_combined %>%
  filter(task_comments == "USE RM" | task_comments == "USE LM")



### Now, join comments and replace amplitude values in non-working mic
### with those from working mic

# Trim recording_date_time to remove extra space before values
comment_combined <- comment_combined %>%
  mutate(recording_date_time = str_trim(as.character(recording_date_time)))

# Join dfs and replace amp values
dat1_fixed <- main_tag_df_2 %>%
  # Join on the flagged recordings
  left_join(comment_combined,
            by = c("location", "recording_date_time")) %>%
  mutate(
    # LEFT mic: only change it when we are told "USE RM" (right mic is good)
    left_full_freq_tag_rms_peak_dbfs = case_when(
      task_comments == "USE RM" ~ right_full_freq_tag_rms_peak_dbfs,
      TRUE                      ~ left_full_freq_tag_rms_peak_dbfs
    ),
    # RIGHT mic: only change it when we are told "USE LM" (left mic is good)
    right_full_freq_tag_rms_peak_dbfs = case_when(
      task_comments == "USE LM" ~ left_full_freq_tag_rms_peak_dbfs,
      TRUE                      ~ right_full_freq_tag_rms_peak_dbfs
    )
  ) %>%
  # Drop the comments and year
  select(-task_comments)

# Save cleaned tags
write.csv(dat1, "Output/Tabular Data/LA inputs/PGBT_tags_all_years.csv", row.names = FALSE)







######### Create a location, ARU, canopy openness df for LA processing

# Select location column from previous task_2023 df
LA_metadata <- task_2023 %>%
  select(location) %>%
  distinct(location) %>%
  mutate("SM2" = "0") %>%
  mutate("canopy" = "1") 


# Save
write.csv(LA_metadata, "Output/Tabular Data/LA inputs/metadata.csv")
