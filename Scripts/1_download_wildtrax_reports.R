# ---
# title: "Download WT reports"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: 
# ---


## Load packages----
# install.packages("remotes")
#library(remotes)
#remotes::install_github("ABbiodiversity/wildRtrax")
# remotes::install_github("ABbiodiversity/wildRtrax@development")
library(wildrtrax)
library(tidyverse)

##////////////////////////////////////////////////////////////////
# Login to WildTrax----
# NOTE: Edit the 'loginexample.R' script to include your WildTrax 
# login details and rename to 'login.R'. 
# DO NOT PUSH YOUR LOGIN TO GITHUB
config <- "Scripts/login.R"
source(config)
wt_auth()

# Get WildTrax project summary
wildtrax_projects <- wt_get_download_summary(
  sensor_id = "ARU"
)


# Download PGBT reports
PGBT_report_2023 <- wt_download_report(
  project_id = 2346,
  sensor_id = "ARU",
  report = "main",
  weather_cols = FALSE
)
  
PGBT_report_2024 <- wt_download_report(
  project_id = 2967,
  sensor_id = "ARU",
  report = "main",
  weather_cols = FALSE
)

PGBT_report_2025 <- wt_download_report(
  project_id = ####,
  sensor_id = "ARU",
  report = "main",
  weather_cols = FALSE
)

# Save----
save(PGBT_report, file=paste0("Output/R Data/PGBT_raw_aru_", Sys.Date(), ".rData"))
