# ---
# title: "Exploratory plots"
# author: "Leonard Patterson"
# created: "2025-12-12"
# description: 
# ---

# Load packages
library(tidyverse)



##### LIMITED AMPLITUDE

# Load data
LA_data <- read.csv("Output/Tabular Data/mean_count_all_years_LA.csv")

# --- 1. Data Preparation and Reshaping ---

# The original 'aggregated_counts' columns were: Total_Count, treatment, year, Species.
# We need to transform LA_data to match this 'long' structure.
LA_data_long <- LA_data %>%
  
  # A. Extract the Treatment type from the 'location' column.
  # Assuming the treatment is always the 4th and 5th characters after the first hyphen (e.g., DOC213-FR-C -> FR)
  # A safer assumption is often to use regex, but based on the sample, 'FR', 'HR', 'LR', 'NT' are the treatments.
  # Let's assume the treatment is the two characters immediately before the second hyphen.
  # Example: DOC213-FR-C -> "FR"
  mutate(treatment = sub(".*-(\\w{2})-.*", "\\1", location)) %>%
  
  # B. Reshape the data from 'wide' to 'long' format.
  # We gather all species columns (excluding the first three: location, year, treatment)
  # into two new columns: 'Species' (key) and 'Total_Count' (value).
  pivot_longer(
    cols = c(ALFL:YRWA),        # Selects all columns from ALFL up to YRWA
    names_to = "Species",       # New column for the species names (e.g., ALFL, AMRE)
    values_to = "Total_Count"   # New column for the count values
  ) %>%
  
  # C. Filter out zero counts for cleaner plotting (optional, but often helpful)
  filter(Total_Count > 0)


# --- 2. Create the Bar Plot using ggplot2 ---

ggplot(LA_data_long, 
       # Treatment is on the X axis, Year is the FILL color
       aes(x = treatment, y = Total_Count, fill = factor(year))) + # 'year' must be a factor for fill
  
  geom_bar(
    stat = "identity",       
    position = "dodge",      # Use position="dodge" to put the years side-by-side
    width = 0.8              
  ) +
  
  # Create Facets (One plot per species)
  facet_wrap(
    ~ Species,               
    scales = "free_y",       # Allows Y-axis scale to vary per species panel
    ncol = 4                 
  ) +
  
  labs(
    title = "Total Species Abundance by Treatment and Year",
    x = "Treatment Type",
    y = "Total Aggregated Count",
    fill = "Year" # Legend title for the fill color
  ) +
  
  # CRITICAL FIX: Use a distinct color palette for the years (fill)
  scale_fill_brewer(palette = "Set1") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # Keep X-axis text centered under the treatment groups
    axis.text.x = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )







##### UNLIMITED DISTANCE 

# Load data
UD_data <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

# --- 1. Data Preparation and Reshaping ---

# The original 'aggregated_counts' columns were: Total_Count, treatment, year, Species.
# We need to transform LA_data to match this 'long' structure.
UD_data_long <- UD_data %>%
  
  # A. Extract the Treatment type from the 'location' column.
  # Assuming the treatment is always the 4th and 5th characters after the first hyphen (e.g., DOC213-FR-C -> FR)
  # A safer assumption is often to use regex, but based on the sample, 'FR', 'HR', 'LR', 'NT' are the treatments.
  # Let's assume the treatment is the two characters immediately before the second hyphen.
  # Example: DOC213-FR-C -> "FR"
  mutate(treatment = sub(".*-(\\w{2})-.*", "\\1", location)) %>%
  
  # B. Reshape the data from 'wide' to 'long' format.
  # We gather all species columns (excluding the first three: location, year, treatment)
  # into two new columns: 'Species' (key) and 'Total_Count' (value).
  pivot_longer(
    cols = c(ALFL:MAWA),        # Selects all columns from ALFL up to YRWA
    names_to = "Species",       # New column for the species names (e.g., ALFL, AMRE)
    values_to = "Total_Count"   # New column for the count values
  ) %>%
  
  # C. Filter out zero counts for cleaner plotting (optional, but often helpful)
  filter(Total_Count > 0)


# --- 2. Create the Bar Plot using ggplot2 ---

ggplot(UD_data_long, 
       # Treatment is on the X axis, Year is the FILL color
       aes(x = treatment, y = Total_Count, fill = factor(year))) + # 'year' must be a factor for fill
  
  geom_bar(
    stat = "identity",       
    position = "dodge",      # Use position="dodge" to put the years side-by-side
    width = 0.8              
  ) +
  
  # Create Facets (One plot per species)
  facet_wrap(
    ~ Species,               
    scales = "free_y",       # Allows Y-axis scale to vary per species panel
    ncol = 4                 
  ) +
  
  labs(
    title = "Total Species Abundance by Treatment and Year",
    x = "Treatment Type",
    y = "Total Aggregated Count",
    fill = "Year" # Legend title for the fill color
  ) +
  
  # CRITICAL FIX: Use a distinct color palette for the years (fill)
  scale_fill_brewer(palette = "Set1") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # Keep X-axis text centered under the treatment groups
    axis.text.x = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
