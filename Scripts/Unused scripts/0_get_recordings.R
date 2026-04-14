# Load packages

library(dplyr)
library(stringr)
library(lubridate)
library(readr)
library(tidyr)

# Load data

recordings_2023 <- read.csv("Input/CANFOR_2023_recordings.csv")
recordings_2024 <- read.csv("Input/CANFOR_2024_recordings.csv")

# Remove unneeded columns

recordings_2023 <- recordings_2023 %>%
  select(Name)

recordings_2024 <- recordings_2024 %>%
  select(Name)






######################################## 2023
########################################
########################################

recordings_2023_processed <- recordings_2023 %>%
  mutate(
    # Use str_match to extract all components at once for .wav files
    # This avoids the lookbehind issue by capturing the relevant parts as groups.
    # The regex pattern is:
    # ^                            - Start of the string
    # ([A-Za-z0-9-]+)             - Group 1: Site (e.g., MUS061A-HR-C)
    # _                            - Literal underscore
    # (\\d+\\+?\\d*)              - Group 2: The variable part (e.g., 0+1). This was the problematic part in lookbehind.
    # _                            - Literal underscore
    # (\\d{8})                     - Group 3: Date (8 digits, e.g., 20230601)
    # _                            - Literal underscore
    # (\\d{6})                     - Group 4: Time (6 digits, e.g., 061500)
    # \\.wav$                      - Literal ".wav" at the end of the string
    matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$")
  ) %>%
  mutate(
    # Assign extracted components to new columns from the 'matches' matrix
    # If the regex doesn't match (e.g., for 'YIBGW2.csv'), 'matches' will contain NAs,
    # and these assignments will correctly result in NA for the respective columns.
    Site = matches[, 2], # Group 1 is the Site
    Date = matches[, 4], # Group 3 is the Date
    Time = matches[, 5]  # Group 4 is the Time
  ) %>%
  # Now create Year, Month, Day, Hour from the extracted Date and Time columns
  mutate(
    # Step 5: Create 'Year', 'Month', 'Day' from the 'Date' column
    # str_sub extracts substrings based on start and end positions
    Year = if_else(!is.na(Date), as.integer(str_sub(Date, 1, 4)), NA_integer_),
    Month = if_else(!is.na(Date), str_sub(Date, 5, 6), NA_character_),
    Day = if_else(!is.na(Date), str_sub(Date, 7, 8), NA_character_),
    
    # Step 6: Create 'Hour' column by converting the first two digits of 'Time'
    Hour = if_else(!is.na(Time), as.integer(str_sub(Time, 1, 2)), NA_integer_)
  ) %>%
  select(-matches) # Remove the temporary 'matches' column, as it's no longer needed

recordings_2023_processed <- recordings_2023_processed[-1,]




################## Get 2023 date ranges

recordings_2023_processed <- recordings_2023 %>%
  mutate(
    # Perform str_match once and store the result in a temporary column '.matches'
    # This pattern only applies to .wav files; non-matching rows will have NA for .matches.
    # The regex captures Site, a variable part, Date, and Time.
    .matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$"),
    
    # Assign extracted components to new columns from the '.matches' matrix
    Site = .matches[, 2], # Group 1 is the Site
    Date = .matches[, 4], # Group 3 is the Date (YYYYMMDD)
    Time = .matches[, 5], # Group 4 is the Time (HHMMSS)
    
    # Create 'Year', 'Month', 'Day' from the 'Date' column
    Year = if_else(!is.na(Date), as.integer(str_sub(Date, 1, 4)), NA_integer_),
    Month = if_else(!is.na(Date), str_sub(Date, 5, 6), NA_character_),
    Day = if_else(!is.na(Date), str_sub(Date, 7, 8), NA_character_),
    
    # Create 'Hour' column by converting the first two digits of 'Time'
    Hour = if_else(!is.na(Time), as.integer(str_sub(Time, 1, 2)), NA_integer_)
  ) %>%
  select(-.matches) # Remove the temporary '.matches' column

# 4. Calculate min and max date range for each unique Site
# First, ensure 'Date' is in a proper Date format for correct min/max calculation
# Then, group by Site and summarize to find min and max dates.
site_date_ranges_2023 <- recordings_2023_processed %>%
  # Filter out rows where Site or Date might be NA (e.g., the .csv file)
  filter(!is.na(Site) & !is.na(Date)) %>%
  # Convert the Date column to a proper Date object for accurate min/max
  mutate(Date_Parsed = as.Date(Date, format = "%Y%m%d")) %>%
  group_by(Site) %>%
  summarise(
    MinDate = min(Date_Parsed, na.rm = TRUE),
    MaxDate = max(Date_Parsed, na.rm = TRUE),
    .groups = 'drop' # Drop the grouping structure after summarizing
  )

## Calculate date range
site_date_ranges_2023 <- site_date_ranges_2023 %>%
  mutate(
    MinDate = ymd(MinDate), # Ensure MinDate is in date format
    MaxDate = ymd(MaxDate), # Ensure MaxDate is in date format
    recording_date_range = as.numeric(MaxDate - MinDate)
  )





######## Find overlap dates between sites
########

# Calculate the common recording range for all sites
common_start_date_val <- max(site_date_ranges_2023$MinDate)
common_end_date_val <- min(site_date_ranges_2023$MaxDate)

# Add the common start and end dates as new columns to the dataframe
site_date_ranges_2023 <- site_date_ranges_2023 %>%
  mutate(
    CommonRecordingStartDate = common_start_date_val,
    CommonRecordingEndDate = common_end_date_val
  )

write.csv(site_date_ranges_2023, "Output/Tabular Data/CANFOR_2023_recording_ranges.csv")










###### Get 2023 recordings (3 per site)

# Step 1: Keep rows with relevant dates and hour
valid_dates <- c(20230602, 20230603, 20230605, 20230608)

filtered <- recordings_2023_processed %>%
  filter(Date %in% valid_dates & Hour == 5)

# Step 2: Create preferred date groups
filtered <- filtered %>%
  mutate(
    date_group = case_when(
      Date %in% c(20230602, 20230603) ~ "A",
      Date == 20230605 ~ "B",
      Date == 20230608 ~ "C",
      TRUE ~ NA_character_
    )
  )

# Step 3: Prioritize earlier date within group A (02 preferred over 03)
filtered_grouped <- filtered %>%
  group_by(Site, date_group) %>%
  arrange(Date) %>%  # ensures 20230602 is selected over 20230603
  slice(1) %>%       # keep only the first (preferred) row
  ungroup()

# Step 4: Ensure all 3 date groups exist per site
sites_with_all_three <- filtered_grouped %>%
  count(Site) %>%
  filter(n == 3) %>%
  pull(Site)

# Final result: only sites that have all 3 required recordings
recordings_2023_upload <- filtered_grouped %>%
  filter(Site %in% sites_with_all_three)






######## Pull recordings from Cirrus and save on desktop

# Define soure folder
recordings_2023_upload <- recordings_2023_upload %>%
  dplyr::mutate(
    source_path = file.path("Z:/CANFOR/ARU/IBWG/2023/V1", Site, Name)
  )

# Define destination folder
dest_folder <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"

# Create folder if it doesn't exist
if (!dir.exists(dest_folder)) {
  dir.create(dest_folder, recursive = TRUE)
}

# Copy files
file.copy(
  from = recordings_2023_upload$source_path,
  to = file.path(dest_folder, recordings_2023_upload$Name),
  overwrite = FALSE
)




### Check if overlap between recordings pulled and those already transcribed

# Define paths
desktop_folder <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"
recordings_on_desktop <- list.files(desktop_folder)

# STEP 1: Identify and remove recordings that match transcribed names
matching_names <- recordings_on_desktop[
  recordings_on_desktop %in% transcribed_recordings_2023$source_file_name
]

# Remove matched files from Desktop
file.remove(file.path(desktop_folder, matching_names))

# STEP 2: Get site/date info for matched recordings
matched_info <- recordings_2023_processed %>%
  filter(Name %in% matching_names) %>%
  select(Site, Date)

# STEP 3: Find Hour 6 replacements for same Site and Date
replacement_candidates <- recordings_2023_processed %>%
  filter(Hour == 6) %>%
  inner_join(matched_info, by = c("Site", "Date"))

# Keep one replacement per Site-Date pair
replacement_files <- replacement_candidates %>%
  group_by(Site, Date) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    source_path = file.path("Z:/CANFOR/ARU/IBWG/2023/V1", Site, Name),
    dest_path = file.path(desktop_folder, Name)
  )

# STEP 4: Copy replacement files to Desktop folder
file.copy(from = replacement_files$source_path, to = replacement_files$dest_path, overwrite = FALSE)

# STEP 5: Confirm no remaining matches
final_files <- list.files(desktop_folder)
final_matches <- final_files %in% transcribed_recordings_2023$source_file_name
cat("Final number of matches after replacement:", sum(final_matches), "\n")

# Optional: list remaining matches (should be zero)
final_matching_files <- final_files[final_matches]











######################################## 2024
########################################
########################################

recordings_2024_processed <- recordings_2024 %>%
  dplyr::mutate(
    # Use str_match to extract all components at once for .wav files
    # This avoids the lookbehind issue by capturing the relevant parts as groups.
    # The regex pattern is:
    # ^                            - Start of the string
    # ([A-Za-z0-9-]+)             - Group 1: Site (e.g., MUS061A-HR-C)
    # _                            - Literal underscore
    # (\\d+\\+?\\d*)              - Group 2: The variable part (e.g., 0+1). This was the problematic part in lookbehind.
    # _                            - Literal underscore
    # (\\d{8})                     - Group 3: Date (8 digits, e.g., 20230601)
    # _                            - Literal underscore
    # (\\d{6})                     - Group 4: Time (6 digits, e.g., 061500)
    # \\.wav$                      - Literal ".wav" at the end of the string
    matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$")
  ) %>%
  mutate(
    # Assign extracted components to new columns from the 'matches' matrix
    # If the regex doesn't match (e.g., for 'YIBGW2.csv'), 'matches' will contain NAs,
    # and these assignments will correctly result in NA for the respective columns.
    Site = matches[, 2], # Group 1 is the Site
    Date = matches[, 4], # Group 3 is the Date
    Time = matches[, 5]  # Group 4 is the Time
  ) %>%
  # Now create Year, Month, Day, Hour from the extracted Date and Time columns
  mutate(
    # Step 5: Create 'Year', 'Month', 'Day' from the 'Date' column
    # str_sub extracts substrings based on start and end positions
    Year = if_else(!is.na(Date), as.integer(str_sub(Date, 1, 4)), NA_integer_),
    Month = if_else(!is.na(Date), str_sub(Date, 5, 6), NA_character_),
    Day = if_else(!is.na(Date), str_sub(Date, 7, 8), NA_character_),
    
    # Step 6: Create 'Hour' column by converting the first two digits of 'Time'
    Hour = if_else(!is.na(Time), as.integer(str_sub(Time, 1, 2)), NA_integer_)
  ) %>%
  select(-matches) # Remove the temporary 'matches' column, as it's no longer needed

recordings_2024_processed <- recordings_2024_processed[-1,]




################## Get 2024 date ranges

recordings_2024_processed <- recordings_2024 %>%
  mutate(
    # Perform str_match once and store the result in a temporary column '.matches'
    # This pattern only applies to .wav files; non-matching rows will have NA for .matches.
    # The regex captures Site, a variable part, Date, and Time.
    .matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$"),
    
    # Assign extracted components to new columns from the '.matches' matrix
    Site = .matches[, 2], # Group 1 is the Site
    Date = .matches[, 4], # Group 3 is the Date (YYYYMMDD)
    Time = .matches[, 5], # Group 4 is the Time (HHMMSS)
    
    # Create 'Year', 'Month', 'Day' from the 'Date' column
    Year = if_else(!is.na(Date), as.integer(str_sub(Date, 1, 4)), NA_integer_),
    Month = if_else(!is.na(Date), str_sub(Date, 5, 6), NA_character_),
    Day = if_else(!is.na(Date), str_sub(Date, 7, 8), NA_character_),
    
    # Create 'Hour' column by converting the first two digits of 'Time'
    Hour = if_else(!is.na(Time), as.integer(str_sub(Time, 1, 2)), NA_integer_)
  ) %>%
  select(-.matches) # Remove the temporary '.matches' column

# 3. Print the resulting dataframe to see the new columns (optional, for debugging)
# print(recordings_2023_processed)

# 4. Calculate min and max date range for each unique Site
# First, ensure 'Date' is in a proper Date format for correct min/max calculation
# Then, group by Site and summarize to find min and max dates.
site_date_ranges_2024 <- recordings_2024_processed %>%
  # Filter out rows where Site or Date might be NA (e.g., the .csv file)
  filter(!is.na(Site) & !is.na(Date)) %>%
  # Convert the Date column to a proper Date object for accurate min/max
  mutate(Date_Parsed = as.Date(Date, format = "%Y%m%d")) %>%
  group_by(Site) %>%
  summarise(
    MinDate = min(Date_Parsed, na.rm = TRUE),
    MaxDate = max(Date_Parsed, na.rm = TRUE),
    .groups = 'drop' # Drop the grouping structure after summarizing
  )












### Select 3 recordings per site for 2024 (MUS124-FR-C is missing many recordings, hence 
### this site having different dates and time than the other sites)

# Step 1: Define custom selection for MUS124-FR-C
custom_mus124 <- recordings_2024_processed %>%
  filter(
    Site == "MUS124-FR-C" &
      ((Date == 20240620 & Hour == 4) |
         (Date == 20240623 & Hour == 5) |
         (Date == 20240627 & Hour == 5))
  )

# Step 2: Define default selection for all other sites
default_dates <- c(20240604, 20240607, 20240610)

default_selection <- recordings_2024_processed %>%
  filter(Site != "MUS124-FR-C", Date %in% default_dates, Hour == 5)

# Step 3: Combine selections
recordings_2024_upload <- bind_rows(custom_mus124, default_selection) %>%
  mutate(
    source_path = file.path("Z:/CANFOR/ARU/IBWG/2024/V1", Site, Name)
  )

# Step 4: Define destination folder
dest_folder <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"

# Create folder if it doesn't exist
if (!dir.exists(dest_folder)) {
  dir.create(dest_folder, recursive = TRUE)
}

# Step 5: Copy files
file.copy(
  from = recordings_2024_upload$source_path,
  to = file.path(dest_folder, recordings_2024_upload$Name),
  overwrite = FALSE
)

# Step 6: Check which files were successfully copied
successful_copies <- recordings_2024_upload$Name[
  file.exists(file.path(dest_folder, recordings_2024_upload$Name))
]

cat("Copied", length(successful_copies), "files to Desktop folder.\n")









####### Check for overlap between new recordings and those already transcribed

# Load names of already transcribed recordings
transcribed_recordings_2024 <- read.csv("Input/CANFOR_Interior_Broadleaf_Working_Group_ARU_Monitoring_2024_recording_report.csv")

# Step 1: List files in the CANFOR_2024 folder
desktop_folder <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"
recordings_on_desktop <- list.files(desktop_folder)

# Step 2: Check for matches with transcribed_recordings_2023$source_file_name
matched_recordings <- recordings_on_desktop %in% transcribed_recordings_2024$source_file_name

# Step 3: Create summary
num_matches <- sum(matched_recordings)
total_files <- length(recordings_on_desktop)

cat("Number of matches found:", num_matches, "of", total_files, "files.\n")

# Optional: list which recordings matched
matching_names <- recordings_on_desktop[matched_recordings]











######### Identify missing recordings

library(dplyr)
library(tidyr)

# Define the required dates
required_dates <- c(20240604, 20240607, 20240610)

# Step 1: Create a complete grid of expected recordings (24 sites × 3 dates)
# Get the list of all sites from the full recordings data
all_sites <- recordings_2024_processed %>%
  filter(Hour == 5 & Date %in% required_dates) %>%
  distinct(Site) %>%
  pull(Site)

expected_grid <- expand.grid(
  Site = all_sites,
  Date = required_dates
)

# Step 2: Get actual recordings that were filtered and copied
recordings_2024_upload$Date <- as.numeric(recordings_2024_upload$Date)
actual_grid <- recordings_2024_upload %>%
  filter(Hour == 5 & Date %in% required_dates) %>%
  select(Site, Date)

# Step 3: Identify missing recordings by anti_join
missing_recordings <- expected_grid %>%
  anti_join(actual_grid, by = c("Site", "Date"))

# Step 4: Count how many recordings each site is missing
missing_summary <- missing_recordings %>%
  group_by(Site) %>%
  summarise(MissingDates = paste(Date, collapse = ", "),
            MissingCount = n(), .groups = "drop")

# Output
print(missing_summary)

cat("Number of sites missing at least one recording:", nrow(missing_summary), "\n")
cat("Total missing recordings:", nrow(missing_recordings), "\n")





# No recordings missing; therefore, an entire site is missing. This code determines which site 
# is missing recordings

library(dplyr)
library(stringr)

# Step 1: Get site list from CANFOR_2023 filenames
canfor_2023_files <- list.files("C:/Users/leona/OneDrive/Desktop/CANFOR_2023")

# Extract site names from file names (assumes format: SITE_...wav)
canfor_2023_sites <- str_match(canfor_2023_files, "^([A-Za-z0-9-]+)_")[, 2] %>%
  unique() %>%
  sort()

# Step 2: Get site list from 2024 upload
canfor_2024_files <- list.files("C:/Users/leona/OneDrive/Desktop/CANFOR_2024")
canfor_2024_sites <- str_match(canfor_2024_files, "^([A-Za-z0-9-]+)_")[, 2] %>%
  unique() %>%
  sort()

# Step 3: Find the site that is in 2023 but missing in 2024
missing_2024_site <- setdiff(canfor_2023_sites, canfor_2024_sites)

# Output result
cat("Site(s) present in CANFOR_2023 but missing in CANFOR_2024:\n")
print(missing_2024_site)













############ Get more recordings (b/c of bad weather or malfunctions)

raw_dat2023 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2023_Tasks_202521.csv")
raw_dat2024 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2024_Tasks_202521.csv")

# Filter to only LP tasks

dat2023 <- raw_dat2023 %>%
  filter(observer == "Leonard Patterson")
dat2024 <- raw_dat2024 %>%
  filter(observer == "Leonard Patterson")


# --- Catalog inputs (file name lists with column Name) ---
in_2023 <- "Input/CANFOR_2023_recordings.csv"
in_2024 <- "Input/CANFOR_2024_recordings.csv"

# --- Storage roots and destinations ---
src_root_2023 <- "Z:/CANFOR/ARU/IBWG/2023/V1"
src_root_2024 <- "Z:/CANFOR/ARU/IBWG/2024/V1"
dest_2023 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"
dest_2024 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"

dir.create(dest_2023, recursive = TRUE, showWarnings = FALSE)
dir.create(dest_2024, recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)


# --- Parse catalog filenames into joinable keys ---
parse_recordings <- function(df, year_root) {
  df %>%
    mutate(matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$")) %>%
    transmute(
      Name,
      Site  = matches[, 2],
      Date8 = matches[, 4],
      Time6 = matches[, 5],
      Hour  = suppressWarnings(as.integer(str_sub(matches[, 5], 1, 2))),
      Date_only = as.Date(suppressWarnings(ymd(matches[, 4]))),
      DateTime = suppressWarnings(ymd_hms(
        paste0(matches[, 4], " ",
               str_sub(matches[, 5], 1, 2), ":",
               str_sub(matches[, 5], 3, 4), ":",
               str_sub(matches[, 5], 5, 6)))),
      source_path = ifelse(!is.na(Name) & !is.na(Site),
                           file.path(year_root, Site, Name),
                           NA_character_)
    )
}

stopifnot(file.exists(in_2023))
cat_2023 <- read_csv(in_2023, show_col_types = FALSE) %>%
  select(Name) %>% distinct() %>% parse_recordings(src_root_2023)

cat_2024 <- if (file.exists(in_2024)) {
  read_csv(in_2024, show_col_types = FALSE) %>%
    select(Name) %>% distinct() %>% parse_recordings(src_root_2024)
} else tibble()

# --- Blocklist: anything already transcribed by NON-LP (avoid duplicates) ---
already_transcribed_by_others <- function(raw_df) {
  raw_df %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      Hour = hour(dt),
      Date_only = as.Date(dt)
    ) %>%
    filter(status_clean == "transcribed",
           !is.na(transcriber),
           transcriber != "Leonard Patterson") %>%
    filter(!is.na(location), !is.na(Date_only), !is.na(Hour)) %>%
    distinct(Site = location, Date_only, Hour)
}

# --- Core builder: next day same hour; if blocked/missing ⇒ next day hour+1 ---
build_replacements_nextday_then_plus1 <- function(dat_lp, raw_all, catalog_df, year_label) {
  if (nrow(dat_lp) == 0) return(tibble())
  
  blocked <- already_transcribed_by_others(raw_all)
  
  failed <- dat_lp %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      is_fail = status_clean %in% c("bad weather", "malfunction"),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      BaseDate = as.Date(dt),
      BaseHour = hour(dt),
      NextDate = BaseDate + days(1L)
    ) %>%
    filter(is_fail, !is.na(location), !is.na(NextDate), !is.na(BaseHour)) %>%
    mutate(req_id = row_number(),
           Site = location)
  
  if (nrow(failed) == 0) return(tibble())
  
  # Candidate 1: next day, same hour (skip if blocked by others)
  cand1 <- failed %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_same_hour",
           CandidateDate = NextDate,
           CandidateHour = BaseHour)
  
  # Which req_ids are satisfied by cand1 (found a Name)?
  ids1 <- cand1 %>% filter(!is.na(Name)) %>% pull(req_id)
  
  # Candidate 2: next day, hour+1 (0..23), also skip if blocked by others
  fallback <- failed %>%
    filter(!(req_id %in% ids1)) %>%
    mutate(FallbackHour = BaseHour + 1L) %>%
    filter(FallbackHour >= 0L, FallbackHour <= 23L)
  
  cand2 <- fallback %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_hour_plus1",
           CandidateDate = NextDate,
           CandidateHour = FallbackHour)
  
  # Combine best available; if neither found, emit a placeholder row for logging
  best <- bind_rows(
    cand1 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice),
    cand2 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice)
  ) %>%
    group_by(req_id) %>%
    arrange(match(choice, c("next_day_same_hour", "next_day_hour_plus1")),
            is.na(Name), DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  missing_ids <- setdiff(failed$req_id, best$req_id)
  if (length(missing_ids)) {
    best <- bind_rows(
      best,
      failed %>%
        filter(req_id %in% missing_ids) %>%
        transmute(req_id, Site,
                  CandidateDate = NextDate,
                  CandidateHour = BaseHour,
                  Name = NA_character_,
                  source_path = NA_character_,
                  DateTime = as.POSIXct(NA),
                  choice = "none")
    )
  }
  
  best %>%
    mutate(year = year_label,
           found = !is.na(Name),
           dest_path = NA_character_)
}

# --- Copier with safety checks ---
do_copy <- function(df, dest_folder) {
  if (nrow(df) == 0) return(invisible(df))
  df <- df %>%
    mutate(
      dest_path = ifelse(!is.na(Name), file.path(dest_folder, Name), NA_character_),
      source_exists = ifelse(!is.na(source_path), file.exists(source_path), FALSE),
      already_on_desktop = ifelse(!is.na(dest_path), file.exists(dest_path), FALSE)
    )
  mask <- df$source_exists & !df$already_on_desktop & !is.na(df$dest_path)
  copied <- rep(FALSE, nrow(df))
  if (any(mask)) {
    copied[mask] <- file.copy(from = df$source_path[mask],
                              to   = df$dest_path[mask],
                              overwrite = FALSE)
  }
  df$copied <- copied
  df
}


# =========================
# Run for 2023 (needs: raw_dat2023, dat2023)
# =========================

stopifnot(exists("raw_dat2023"), exists("dat2023"))

repl_2023 <- build_replacements_nextday_then_plus1(
  dat_lp     = dat2023,
  raw_all    = raw_dat2023,
  catalog_df = cat_2023,
  year_label = "2023"
)
repl_2023_done <- do_copy(repl_2023, dest_2023)

cat("\n--- 2023 Replacement Summary (next day same hour; else +1h) ---\n")
cat("LP failures considered:     ",
    nrow(dat2023 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
cat("Matched (any choice):       ", sum(repl_2023_done$found, na.rm = TRUE), "\n")
cat("Source exists on disk:      ", sum(repl_2023_done$source_exists, na.rm = TRUE), "\n")
cat("Already on Desktop:         ", sum(repl_2023_done$already_on_desktop, na.rm = TRUE), "\n")
cat("Copied now:                 ", sum(repl_2023_done$copied, na.rm = TRUE), "\n")

write_csv(repl_2023_done %>%
            transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                      Name, source_path, dest_path,
                      source_exists, already_on_desktop, copied),
          "Output/Tabular Data/replacement_copy_log_2023.csv")

if (any(!repl_2023_done$found | !repl_2023_done$source_exists)) {
  cat("\nUnfilled or missing examples (up to 10):\n")
  print(repl_2023_done %>%
          filter(!found | !source_exists) %>%
          select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
          head(10))
}

# =========================
# Run for 2024 (if available)
# =========================
if (exists("raw_dat2024") && exists("dat2024") && nrow(cat_2024)) {
  repl_2024 <- build_replacements_nextday_then_plus1(
    dat_lp     = dat2024,
    raw_all    = raw_dat2024,
    catalog_df = cat_2024,
    year_label = "2024"
  )
  repl_2024_done <- do_copy(repl_2024, dest_2024)
  
  cat("\n--- 2024 Replacement Summary (next day same hour; else +1h) ---\n")
  cat("LP failures considered:     ",
      nrow(dat2024 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
  cat("Matched (any choice):       ", sum(repl_2024_done$found, na.rm = TRUE), "\n")
  cat("Source exists on disk:      ", sum(repl_2024_done$source_exists, na.rm = TRUE), "\n")
  cat("Already on Desktop:         ", sum(repl_2024_done$already_on_desktop, na.rm = TRUE), "\n")
  cat("Copied now:                 ", sum(repl_2024_done$copied, na.rm = TRUE), "\n")
  
  write_csv(repl_2024_done %>%
              transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                        Name, source_path, dest_path,
                        source_exists, already_on_desktop, copied),
            "Output/Tabular Data/replacement_copy_log_2024.csv")
  
  if (any(!repl_2024_done$found | !repl_2024_done$source_exists)) {
    cat("\nUnfilled or missing examples (up to 10):\n")
    print(repl_2024_done %>%
            filter(!found | !source_exists) %>%
            select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
            head(10))
  }
} else {
  message("Skipped 2024 (missing raw_dat2024/dat2024 or 2024 catalog).")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

############ Get more recordings (b/c of bad weather or malfunctions)

raw_dat2023 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2023_Tasks_20250819.csv")
raw_dat2024 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2024_Tasks_20250819.csv")

# Filter to only LP tasks

dat2023 <- raw_dat2023 %>%
  filter(transcriber == "Leonard Patterson")
dat2024 <- raw_dat2024 %>%
  filter(transcriber == "Leonard Patterson")


# --- Catalog inputs (file name lists with column Name) ---
in_2023 <- "Input/CANFOR_2023_recordings.csv"
in_2024 <- "Input/CANFOR_2024_recordings.csv"

# --- Storage roots and destinations ---
src_root_2023 <- "Z:/CANFOR/ARU/IBWG/2023/V1"
src_root_2024 <- "Z:/CANFOR/ARU/IBWG/2024/V1"
dest_2023 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"
dest_2024 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"

dir.create(dest_2023, recursive = TRUE, showWarnings = FALSE)
dir.create(dest_2024, recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)


# --- Parse catalog filenames into joinable keys ---
parse_recordings <- function(df, year_root) {
  df %>%
    mutate(matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$")) %>%
    transmute(
      Name,
      Site  = matches[, 2],
      Date8 = matches[, 4],
      Time6 = matches[, 5],
      Hour  = suppressWarnings(as.integer(str_sub(matches[, 5], 1, 2))),
      Date_only = as.Date(suppressWarnings(ymd(matches[, 4]))),
      DateTime = suppressWarnings(ymd_hms(
        paste0(matches[, 4], " ",
               str_sub(matches[, 5], 1, 2), ":",
               str_sub(matches[, 5], 3, 4), ":",
               str_sub(matches[, 5], 5, 6)))),
      source_path = ifelse(!is.na(Name) & !is.na(Site),
                           file.path(year_root, Site, Name),
                           NA_character_)
    )
}

stopifnot(file.exists(in_2023))
cat_2023 <- read_csv(in_2023, show_col_types = FALSE) %>%
  select(Name) %>% distinct() %>% parse_recordings(src_root_2023)

cat_2024 <- if (file.exists(in_2024)) {
  read_csv(in_2024, show_col_types = FALSE) %>%
    select(Name) %>% distinct() %>% parse_recordings(src_root_2024)
} else tibble()

# --- Blocklist: anything already transcribed by NON-LP (avoid duplicates) ---
already_transcribed_by_others <- function(raw_df) {
  raw_df %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      Hour = hour(dt),
      Date_only = as.Date(dt)
    ) %>%
    filter(status_clean == "transcribed",
           !is.na(transcriber),
           transcriber != "Leonard Patterson") %>%
    filter(!is.na(location), !is.na(Date_only), !is.na(Hour)) %>%
    distinct(Site = location, Date_only, Hour)
}

# --- Core builder: next day same hour; if blocked/missing ⇒ next day hour+1 ---
build_replacements_nextday_then_plus1 <- function(dat_lp, raw_all, catalog_df, year_label) {
  if (nrow(dat_lp) == 0) return(tibble())
  
  blocked <- already_transcribed_by_others(raw_all)
  
  failed <- dat_lp %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      is_fail = status_clean %in% c("bad weather", "malfunction"),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      BaseDate = as.Date(dt),
      BaseHour = hour(dt),
      NextDate = BaseDate + days(1L)
    ) %>%
    filter(is_fail, !is.na(location), !is.na(NextDate), !is.na(BaseHour)) %>%
    mutate(req_id = row_number(),
           Site = location)
  
  if (nrow(failed) == 0) return(tibble())
  
  # Candidate 1: next day, same hour (skip if blocked by others)
  cand1 <- failed %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_same_hour",
           CandidateDate = NextDate,
           CandidateHour = BaseHour)
  
  # Which req_ids are satisfied by cand1 (found a Name)?
  ids1 <- cand1 %>% filter(!is.na(Name)) %>% pull(req_id)
  
  # Candidate 2: next day, hour+1 (0..23), also skip if blocked by others
  fallback <- failed %>%
    filter(!(req_id %in% ids1)) %>%
    mutate(FallbackHour = BaseHour + 1L) %>%
    filter(FallbackHour >= 0L, FallbackHour <= 23L)
  
  cand2 <- fallback %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_hour_plus1",
           CandidateDate = NextDate,
           CandidateHour = FallbackHour)
  
  # Combine best available; if neither found, emit a placeholder row for logging
  best <- bind_rows(
    cand1 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice),
    cand2 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice)
  ) %>%
    group_by(req_id) %>%
    arrange(match(choice, c("next_day_same_hour", "next_day_hour_plus1")),
            is.na(Name), DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  missing_ids <- setdiff(failed$req_id, best$req_id)
  if (length(missing_ids)) {
    best <- bind_rows(
      best,
      failed %>%
        filter(req_id %in% missing_ids) %>%
        transmute(req_id, Site,
                  CandidateDate = NextDate,
                  CandidateHour = BaseHour,
                  Name = NA_character_,
                  source_path = NA_character_,
                  DateTime = as.POSIXct(NA),
                  choice = "none")
    )
  }
  
  best %>%
    mutate(year = year_label,
           found = !is.na(Name),
           dest_path = NA_character_)
}

# --- Copier with safety checks ---
do_copy <- function(df, dest_folder) {
  if (nrow(df) == 0) return(invisible(df))
  df <- df %>%
    mutate(
      dest_path = ifelse(!is.na(Name), file.path(dest_folder, Name), NA_character_),
      source_exists = ifelse(!is.na(source_path), file.exists(source_path), FALSE),
      already_on_desktop = ifelse(!is.na(dest_path), file.exists(dest_path), FALSE)
    )
  mask <- df$source_exists & !df$already_on_desktop & !is.na(df$dest_path)
  copied <- rep(FALSE, nrow(df))
  if (any(mask)) {
    copied[mask] <- file.copy(from = df$source_path[mask],
                              to   = df$dest_path[mask],
                              overwrite = FALSE)
  }
  df$copied <- copied
  df
}


# =========================
# Run for 2023 (needs: raw_dat2023, dat2023)
# =========================

stopifnot(exists("raw_dat2023"), exists("dat2023"))

repl_2023 <- build_replacements_nextday_then_plus1(
  dat_lp     = dat2023,
  raw_all    = raw_dat2023,
  catalog_df = cat_2023,
  year_label = "2023"
)
repl_2023_done <- do_copy(repl_2023, dest_2023)

cat("\n--- 2023 Replacement Summary (next day same hour; else +1h) ---\n")
cat("LP failures considered:     ",
    nrow(dat2023 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
cat("Matched (any choice):       ", sum(repl_2023_done$found, na.rm = TRUE), "\n")
cat("Source exists on disk:      ", sum(repl_2023_done$source_exists, na.rm = TRUE), "\n")
cat("Already on Desktop:         ", sum(repl_2023_done$already_on_desktop, na.rm = TRUE), "\n")
cat("Copied now:                 ", sum(repl_2023_done$copied, na.rm = TRUE), "\n")

write_csv(repl_2023_done %>%
            transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                      Name, source_path, dest_path,
                      source_exists, already_on_desktop, copied),
          "Output/Tabular Data/replacement_copy_log_2023.csv")

if (any(!repl_2023_done$found | !repl_2023_done$source_exists)) {
  cat("\nUnfilled or missing examples (up to 10):\n")
  print(repl_2023_done %>%
          filter(!found | !source_exists) %>%
          select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
          head(10))
}

# =========================
# Run for 2024 (if available)
# =========================
if (exists("raw_dat2024") && exists("dat2024") && nrow(cat_2024)) {
  repl_2024 <- build_replacements_nextday_then_plus1(
    dat_lp     = dat2024,
    raw_all    = raw_dat2024,
    catalog_df = cat_2024,
    year_label = "2024"
  )
  repl_2024_done <- do_copy(repl_2024, dest_2024)
  
  cat("\n--- 2024 Replacement Summary (next day same hour; else +1h) ---\n")
  cat("LP failures considered:     ",
      nrow(dat2024 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
  cat("Matched (any choice):       ", sum(repl_2024_done$found, na.rm = TRUE), "\n")
  cat("Source exists on disk:      ", sum(repl_2024_done$source_exists, na.rm = TRUE), "\n")
  cat("Already on Desktop:         ", sum(repl_2024_done$already_on_desktop, na.rm = TRUE), "\n")
  cat("Copied now:                 ", sum(repl_2024_done$copied, na.rm = TRUE), "\n")
  
  write_csv(repl_2024_done %>%
              transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                        Name, source_path, dest_path,
                        source_exists, already_on_desktop, copied),
            "Output/Tabular Data/replacement_copy_log_2024.csv")
  
  if (any(!repl_2024_done$found | !repl_2024_done$source_exists)) {
    cat("\nUnfilled or missing examples (up to 10):\n")
    print(repl_2024_done %>%
            filter(!found | !source_exists) %>%
            select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
            head(10))
  }
} else {
  message("Skipped 2024 (missing raw_dat2024/dat2024 or 2024 catalog).")
}

















library(readr)



############ Get more recordings (b/c of bad weather or malfunctions)

raw_dat2023 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2023_Tasks_20250824.csv")
raw_dat2024 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2024_Tasks_20250824.csv")

# Filter to only LP tasks

dat2023 <- raw_dat2023 %>%
  filter(transcriber == "Leonard Patterson")
dat2024 <- raw_dat2024 %>%
  filter(transcriber == "Leonard Patterson")


# --- Catalog inputs (file name lists with column Name) ---
in_2023 <- "Input/CANFOR_2023_recordings.csv"
in_2024 <- "Input/CANFOR_2024_recordings.csv"

# --- Storage roots and destinations ---
src_root_2023 <- "Z:/CANFOR/ARU/IBWG/2023/V1"
src_root_2024 <- "Z:/CANFOR/ARU/IBWG/2024/V1"
dest_2023 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"
dest_2024 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"

dir.create(dest_2023, recursive = TRUE, showWarnings = FALSE)
dir.create(dest_2024, recursive = TRUE, showWarnings = FALSE)
dir.create("Output/Tabular Data", recursive = TRUE, showWarnings = FALSE)


# --- Parse catalog filenames into joinable keys ---
parse_recordings <- function(df, year_root) {
  df %>%
    mutate(matches = str_match(Name, "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.wav$")) %>%
    transmute(
      Name,
      Site  = matches[, 2],
      Date8 = matches[, 4],
      Time6 = matches[, 5],
      Hour  = suppressWarnings(as.integer(str_sub(matches[, 5], 1, 2))),
      Date_only = as.Date(suppressWarnings(ymd(matches[, 4]))),
      DateTime = suppressWarnings(ymd_hms(
        paste0(matches[, 4], " ",
               str_sub(matches[, 5], 1, 2), ":",
               str_sub(matches[, 5], 3, 4), ":",
               str_sub(matches[, 5], 5, 6)))),
      source_path = ifelse(!is.na(Name) & !is.na(Site),
                           file.path(year_root, Site, Name),
                           NA_character_)
    )
}

stopifnot(file.exists(in_2023))
cat_2023 <- read_csv(in_2023, show_col_types = FALSE) %>%
  select(Name) %>% distinct() %>% parse_recordings(src_root_2023)

cat_2024 <- if (file.exists(in_2024)) {
  read_csv(in_2024, show_col_types = FALSE) %>%
    select(Name) %>% distinct() %>% parse_recordings(src_root_2024)
} else tibble()

# --- Blocklist: anything already transcribed by NON-LP (avoid duplicates) ---
already_transcribed_by_others <- function(raw_df) {
  raw_df %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      Hour = hour(dt),
      Date_only = as.Date(dt)
    ) %>%
    filter(status_clean == "transcribed",
           !is.na(transcriber),
           transcriber != "Leonard Patterson") %>%
    filter(!is.na(location), !is.na(Date_only), !is.na(Hour)) %>%
    distinct(Site = location, Date_only, Hour)
}

# --- Core builder: next day same hour; if blocked/missing ⇒ next day hour+1 ---
build_replacements_nextday_then_plus1 <- function(dat_lp, raw_all, catalog_df, year_label) {
  if (nrow(dat_lp) == 0) return(tibble())
  
  blocked <- already_transcribed_by_others(raw_all)
  
  failed <- dat_lp %>%
    mutate(
      status_clean = str_to_lower(str_trim(as.character(status))),
      is_fail = status_clean %in% c("bad weather", "malfunction"),
      dt = suppressWarnings(ymd_hms(recordingDate)),
      BaseDate = as.Date(dt),
      BaseHour = hour(dt),
      NextDate = BaseDate + days(1L)
    ) %>%
    filter(is_fail, !is.na(location), !is.na(NextDate), !is.na(BaseHour)) %>%
    mutate(req_id = row_number(),
           Site = location)
  
  if (nrow(failed) == 0) return(tibble())
  
  # Candidate 1: next day, same hour (skip if blocked by others)
  cand1 <- failed %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "BaseHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_same_hour",
           CandidateDate = NextDate,
           CandidateHour = BaseHour)
  
  # Which req_ids are satisfied by cand1 (found a Name)?
  ids1 <- cand1 %>% filter(!is.na(Name)) %>% pull(req_id)
  
  # Candidate 2: next day, hour+1 (0..23), also skip if blocked by others
  fallback <- failed %>%
    filter(!(req_id %in% ids1)) %>%
    mutate(FallbackHour = BaseHour + 1L) %>%
    filter(FallbackHour >= 0L, FallbackHour <= 23L)
  
  cand2 <- fallback %>%
    anti_join(blocked, by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    left_join(catalog_df,
              by = c("Site" = "Site", "NextDate" = "Date_only", "FallbackHour" = "Hour")) %>%
    group_by(req_id) %>%
    arrange(DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(choice = "next_day_hour_plus1",
           CandidateDate = NextDate,
           CandidateHour = FallbackHour)
  
  # Combine best available; if neither found, emit a placeholder row for logging
  best <- bind_rows(
    cand1 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice),
    cand2 %>% select(req_id, Site, CandidateDate, CandidateHour, Name, source_path, DateTime, choice)
  ) %>%
    group_by(req_id) %>%
    arrange(match(choice, c("next_day_same_hour", "next_day_hour_plus1")),
            is.na(Name), DateTime, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  missing_ids <- setdiff(failed$req_id, best$req_id)
  if (length(missing_ids)) {
    best <- bind_rows(
      best,
      failed %>%
        filter(req_id %in% missing_ids) %>%
        transmute(req_id, Site,
                  CandidateDate = NextDate,
                  CandidateHour = BaseHour,
                  Name = NA_character_,
                  source_path = NA_character_,
                  DateTime = as.POSIXct(NA),
                  choice = "none")
    )
  }
  
  best %>%
    mutate(year = year_label,
           found = !is.na(Name),
           dest_path = NA_character_)
}

# --- Copier with safety checks ---
do_copy <- function(df, dest_folder) {
  if (nrow(df) == 0) return(invisible(df))
  df <- df %>%
    mutate(
      dest_path = ifelse(!is.na(Name), file.path(dest_folder, Name), NA_character_),
      source_exists = ifelse(!is.na(source_path), file.exists(source_path), FALSE),
      already_on_desktop = ifelse(!is.na(dest_path), file.exists(dest_path), FALSE)
    )
  mask <- df$source_exists & !df$already_on_desktop & !is.na(df$dest_path)
  copied <- rep(FALSE, nrow(df))
  if (any(mask)) {
    copied[mask] <- file.copy(from = df$source_path[mask],
                              to   = df$dest_path[mask],
                              overwrite = FALSE)
  }
  df$copied <- copied
  df
}


# =========================
# Run for 2023 (needs: raw_dat2023, dat2023)
# =========================

stopifnot(exists("raw_dat2023"), exists("dat2023"))

repl_2023 <- build_replacements_nextday_then_plus1(
  dat_lp     = dat2023,
  raw_all    = raw_dat2023,
  catalog_df = cat_2023,
  year_label = "2023"
)
repl_2023_done <- do_copy(repl_2023, dest_2023)

cat("\n--- 2023 Replacement Summary (next day same hour; else +1h) ---\n")
cat("LP failures considered:     ",
    nrow(dat2023 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
cat("Matched (any choice):       ", sum(repl_2023_done$found, na.rm = TRUE), "\n")
cat("Source exists on disk:      ", sum(repl_2023_done$source_exists, na.rm = TRUE), "\n")
cat("Already on Desktop:         ", sum(repl_2023_done$already_on_desktop, na.rm = TRUE), "\n")
cat("Copied now:                 ", sum(repl_2023_done$copied, na.rm = TRUE), "\n")

write_csv(repl_2023_done %>%
            transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                      Name, source_path, dest_path,
                      source_exists, already_on_desktop, copied),
          "Output/Tabular Data/replacement_copy_log_2023.csv")

if (any(!repl_2023_done$found | !repl_2023_done$source_exists)) {
  cat("\nUnfilled or missing examples (up to 10):\n")
  print(repl_2023_done %>%
          filter(!found | !source_exists) %>%
          select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
          head(10))
}

# =========================
# Run for 2024 (if available)
# =========================
if (exists("raw_dat2024") && exists("dat2024") && nrow(cat_2024)) {
  repl_2024 <- build_replacements_nextday_then_plus1(
    dat_lp     = dat2024,
    raw_all    = raw_dat2024,
    catalog_df = cat_2024,
    year_label = "2024"
  )
  repl_2024_done <- do_copy(repl_2024, dest_2024)
  
  cat("\n--- 2024 Replacement Summary (next day same hour; else +1h) ---\n")
  cat("LP failures considered:     ",
      nrow(dat2024 %>% filter(str_to_lower(status) %in% c('bad weather','malfunction'))), "\n")
  cat("Matched (any choice):       ", sum(repl_2024_done$found, na.rm = TRUE), "\n")
  cat("Source exists on disk:      ", sum(repl_2024_done$source_exists, na.rm = TRUE), "\n")
  cat("Already on Desktop:         ", sum(repl_2024_done$already_on_desktop, na.rm = TRUE), "\n")
  cat("Copied now:                 ", sum(repl_2024_done$copied, na.rm = TRUE), "\n")
  
  write_csv(repl_2024_done %>%
              transmute(year, req_id, Site, choice, CandidateDate, CandidateHour,
                        Name, source_path, dest_path,
                        source_exists, already_on_desktop, copied),
            "Output/Tabular Data/replacement_copy_log_2024.csv")
  
  if (any(!repl_2024_done$found | !repl_2024_done$source_exists)) {
    cat("\nUnfilled or missing examples (up to 10):\n")
    print(repl_2024_done %>%
            filter(!found | !source_exists) %>%
            select(Site, choice, CandidateDate, CandidateHour, Name, source_path) %>%
            head(10))
  }
} else {
  message("Skipped 2024 (missing raw_dat2024/dat2024 or 2024 catalog).")
}













################# 2025




# ---- 2025: Gather file names from folders, filter, copy ----

library(dplyr)
library(stringr)
library(lubridate)
library(readr)
library(tidyr)

# Root that contains site subfolders with WAV files
src_root_2025 <- "Z:/CANFOR/ARU/IBWG/2025/V1"

# 1) List all .wav files (case-insensitive), recursively
paths_2025 <- list.files(
  path = src_root_2025,
  pattern = "\\.(?i:wav)$",     # case-insensitive .wav
  full.names = TRUE,
  recursive = TRUE
)

# If there are non-audio extras (e.g., hidden/system), keep only files
paths_2025 <- paths_2025[file.exists(paths_2025)]

# 2) Build a data frame from paths
#    - Name = filename (e.g., "MUS061A-HR-C_0+1_20250604_050000.wav")
#    - SiteFolder = last dir in the path (for sanity check)
#    - Parse Site/Date/Time/Hour from Name via regex
recordings_2025_processed <- tibble(source_path = paths_2025) %>%
  mutate(
    Name       = basename(source_path),
    SiteFolder = basename(dirname(source_path)),
    .m = str_match(
      Name,
      "^([A-Za-z0-9-]+)_(\\d+\\+?\\d*)_(\\d{8})_(\\d{6})\\.(?i:wav)$"
    ),
    Site  = .m[, 2],
    Date  = .m[, 4],   # "YYYYMMDD" (character)
    Time  = .m[, 5],   # "HHMMSS"   (character)
    Hour  = suppressWarnings(as.integer(str_sub(Time, 1, 2))),
    # Optional sanity check: Site from filename should match folder name
    site_match = if_else(!is.na(Site), Site == SiteFolder, NA)
  ) %>%
  select(-.m)

# 3) Keep the three dates at 05:00
wanted_dates <- c("20250604", "20250607", "20250610")

recordings_2025_upload <- recordings_2025_processed %>%
  filter(!is.na(Site), !is.na(Date), !is.na(Hour)) %>%
  filter(Date %in% wanted_dates, Hour == 5)

# 4) Destination and copy
dest_folder_2025 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2025"
if (!dir.exists(dest_folder_2025)) dir.create(dest_folder_2025, recursive = TRUE)

copied_ok <- file.copy(
  from = recordings_2025_upload$source_path,
  to   = file.path(dest_folder_2025, recordings_2025_upload$Name),
  overwrite = FALSE
)

cat("2025: attempted to copy", nrow(recordings_2025_upload),
    "files; successfully copied", sum(copied_ok, na.rm = TRUE), "files.\n")

# 5) Helpful diagnostics

# Any files that matched the filter but don't exist anymore?
missing_sources <- recordings_2025_upload %>%
  filter(!file.exists(source_path)) %>%
  select(Site, Date, Time, Name, source_path)

if (nrow(missing_sources)) {
  cat("Warning: Some matched files were not found on disk. Showing up to 10:\n")
  print(head(missing_sources, 10))
}

# Site-by-site counts
site_counts <- recordings_2025_upload %>%
  count(Site, name = "n_files") %>%
  arrange(Site)

print(site_counts)

# Optional: save a manifest of what you copied (or tried to)
write_csv(
  recordings_2025_upload %>%
    mutate(dest_path = file.path(dest_folder_2025, Name),
           copied = copied_ok),
  "Output/Tabular Data/CANFOR_2025_copy_manifest.csv"
)

# Optional: flag any mismatches where folder name != filename site token
mismatch <- recordings_2025_upload %>%
  filter(!is.na(site_match) & site_match == FALSE) %>%
  select(SiteFolder, Site, Name, source_path)

if (nrow(mismatch)) {
  cat("Note: Some files have Site in filename that differs from folder name. Showing up to 10:\n")
  print(head(mismatch, 10))
}











# Load the necessary library for data manipulation
# If you don't have it installed, run: install.packages("dplyr")
library(dplyr)

dat1 <- read.csv("Input/Add_more_recordings_08-19-2025/Interior Broadleaf Working Group ARU Monitoring 2025 Limited Amplitude Processing_Tags_202527.csv")

dat1 <- dat1 %>%
  filter(observer =="Leonard Patterson") 

# Assuming your data frame is named 'dat1'

# Calculate the number of unique recording_date_time values per location
unique_recording_summary <- dat1 %>%
  # 1. Group the data by the 'location' column
  group_by(location) %>%
  
  # 2. Summarize the data by creating a new column called 'unique_recording_count'
  #    n_distinct() counts the number of unique, non-missing values in the specified column
  summarise(
    unique_recording_count = n_distinct(recording_date_time)
  ) %>%
  
  # 3. Ungroup the data (good practice after using summarise)
  ungroup()

# Print the resulting summary data frame
View(unique_recording_summary)

dat2 <- unique_recording_summary %>%
  filter(unique_recording_count < 3)



##### GET FOUR RECORDINGS

## Interim load
PGBT_2023 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2023_Tasks_202536.csv")
PGBT_2024 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2024_Tasks_202536.csv")
PGBT_2025 <- read.csv("Input/Task reports/Interior Broadleaf Working Group ARU Monitoring 2025 Limited Amplitude Processing_Tasks_202536.csv")



# Filter to only LP as observer
PGBT_2023 <- PGBT_2023 %>%
  filter(observer =="Leonard Patterson") 
PGBT_2024 <- PGBT_2024 %>%
  filter(observer =="Leonard Patterson") 
PGBT_2025 <- PGBT_2025 %>%
  filter(observer =="Leonard Patterson") 




PGBT_2023 <- PGBT_2023 %>%
  select(location, recording_date_time) 
PGBT_2024 <- PGBT_2024 %>%
  select(location, recording_date_time) 
PGBT_2025 <- PGBT_2025 %>%
  select(location, recording_date_time) 






# Load necessary libraries
# We use 'tidyverse' which includes 'dplyr' for data manipulation and 'lubridate' for date/time handling.
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)
library(lubridate)

# --- Define Source and Destination Paths ---
# IMPORTANT: These paths are based on your request. Please ensure they are correct
# and accessible from your R environment before running the code.
src_root_2023 <- "Z:/CANFOR/ARU/IBWG/2023/V1"
src_root_2024 <- "Z:/CANFOR/ARU/IBWG/2024/V1"
src_root_2025 <- "Z:/CANFOR/ARU/IBWG/2025/V1"

dest_2023 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2023"
dest_2024 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2024"
dest_2025 <- "C:/Users/leona/OneDrive/Desktop/CANFOR_2025"

# --- Simulate Data Frames (PGBT_2023, PGBT_2024, PGBT_2025) ---
# NOTE: The dates in your example were 2025, but the dataframe was PGBT_2023.
# We will adjust the simulated data to match the year in the dataframe name (2023, 2024, 2025)
# to ensure the extracted filenames are consistent.








# --- Core Processing Function ---

#' Processes the ARU data to find the target 3-day offset recording and copies the file.
#' The file search now only requires a match on the date (3 days later) and the hour of the original recording.
#' If multiple files match, one is randomly selected.
#'
#' @param df The data frame (PGBT_2023, PGBT_2024, or PGBT_2025).
#' @param src_dir The source directory path (e.g., "Z:/.../V1").
#' @param dest_dir The destination directory path (e.g., "C:/.../CANFOR_2023").
#' @return A data frame containing the list of target recordings and the status of the file copy.
process_aru_data <- function(df, src_dir, dest_dir) {
  # Defensive: make sure there is at least one row
  if (nrow(df) == 0) {
    cat("Data frame is empty, nothing to do.\n")
    return(tibble())
  }
  
  df <- df %>%
    mutate(recording_date_time = lubridate::as_datetime(recording_date_time))
  
  current_year <- unique(lubridate::year(df$recording_date_time[1]))
  cat(paste0("--- Processing data for year: ", current_year, " ---\n"))
  
  # 1. For each location, take the latest recording and compute target date/hour
  target_recordings <- df %>%
    group_by(location) %>%
    slice_max(recording_date_time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      target_date_time = recording_date_time + lubridate::days(4),
      target_date_str  = format(target_date_time, "%Y%m%d"),
      target_hour_str  = format(target_date_time, "%H")
    )
  
  # 2. Ensure destination directory exists
  if (!dir.exists(dest_dir)) {
    cat(paste("Creating destination directory:", dest_dir, "\n"))
    dir.create(dest_dir, recursive = TRUE)
  }
  
  # 3. Loop over rows and do the file search + copy
  out_list <- vector("list", nrow(target_recordings))
  
  for (i in seq_len(nrow(target_recordings))) {
    loc            <- target_recordings$location[i]
    latest_dt      <- target_recordings$recording_date_time[i]
    target_dt      <- target_recordings$target_date_time[i]
    target_dateStr <- target_recordings$target_date_str[i]
    target_hourStr <- target_recordings$target_hour_str[i]
    
    search_dir <- file.path(src_dir, loc)
    
    # Pattern: e.g. MUS111-NT-C_0+1_20230612_05*.wav
    search_pattern <- paste0(
      "^", loc, "_0\\+1_", 
      target_dateStr, "_", 
      target_hourStr,
      ".*\\.wav$"
    )
    
    if (dir.exists(search_dir)) {
      candidates <- list.files(
        path       = search_dir,
        pattern    = search_pattern,
        full.names = FALSE
      )
    } else {
      candidates <- character(0)
    }
    
    if (length(candidates) > 0) {
      # Randomly select one 5am recording if more than one
      target_filename <- sample(candidates, 1)
      src_path  <- file.path(search_dir, target_filename)
      dest_path <- file.path(dest_dir, target_filename)
      
      cat("Attempting to copy:", target_filename, "from", search_dir, "to", dest_dir, "\n")
      copy_success <- file.copy(from = src_path, to = dest_path, overwrite = FALSE)
      file_exists  <- TRUE
    } else {
      target_filename <- NA_character_
      copy_success    <- FALSE
      file_exists     <- FALSE
    }
    
    out_list[[i]] <- tibble(
      site_id          = loc,
      latest_recording = latest_dt,
      target_date_time = target_dt,
      target_filename  = target_filename,
      file_exists      = file_exists,
      copy_success     = copy_success
    )
  }
  
  copy_results <- bind_rows(out_list)
  
  cat("--- Processing Complete ---\n\n")
  return(copy_results)
}


# --- Execute the Function for each Year ---

# 1. Process 2023 Data
results_2023 <- process_aru_data(PGBT_2023, src_root_2023, dest_2023)
print(results_2023)

# 2. Process 2024 Data
results_2024 <- process_aru_data(PGBT_2024, src_root_2024, dest_2024)
print(results_2024)

# 3. Process 2025 Data
results_2025 <- process_aru_data(PGBT_2025, src_root_2025, dest_2025)
print(results_2025)

# --- Summary of results ---
# You can combine all results into a single table:
final_results <- bind_rows(results_2023, results_2024) %>% bind_rows(results_2025)
cat("\n\n--- FINAL SUMMARY OF TARGET FILES AND STATUS ---\n")
print(final_results)

# A final message about files not found
files_not_found <- final_results %>% filter(!file_exists)
if (nrow(files_not_found) > 0) {
  cat("\n\n!! ATTENTION !! The following", nrow(files_not_found), "target files were not found in the source directories:\n")
  print(files_not_found %>% select(site_id, target_filename, file_exists))
} else {
  cat("\n\nAll calculated target files were found in the source directories.\n")
}


PGBT_2023 <- PGBT_2023 %>%
  summarise(unique(location))
