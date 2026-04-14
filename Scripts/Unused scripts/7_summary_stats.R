# ---
# title: "Summary stats"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: 
# ---



############# Summary stats


### UD
# Data 
dat1_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

# Identify species columns in the filtered dataset
species_cols <- setdiff(names(dat1_UD), c("location", "year", "recording_date_time"))

dat1_UD_2 <- dat1_UD %>%
  mutate(
    across(
      all_of(species_cols),
      ~ if_else(is.na(.x), 0L, as.integer(.x > 0))
    )
  )

dat1_UD_3 <- dat1_UD_2 %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)))

dat1_UD_4 <- dat1_UD_3 %>%
  pivot_longer(species_cols)

dat1_UD_5 <- dat1_UD_4 %>%
  mutate(prop_obs = value / 216)



### LA 
# Data 
dat1_LA <- read.csv("Output/Tabular Data/mean_count_all_years_LA.csv")

# Identify species columns in the filtered dataset
species_cols <- setdiff(names(dat1_LA), c("location", "year", "recording_date_time"))

dat1_LA_2 <- dat1_LA %>%
  mutate(
    across(
      all_of(species_cols),
      ~ if_else(is.na(.x), 0L, as.integer(.x > 0))
    )
  )

dat1_LA_3 <- dat1_LA_2 %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)))

dat1_LA_4 <- dat1_LA_3 %>%
  pivot_longer(species_cols)

dat1_LA_5 <- dat1_LA_4 %>%
  mutate(prop_obs = value / 216)

