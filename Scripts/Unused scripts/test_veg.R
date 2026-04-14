# ---
# title: "Singe species analysis - LA - MEAN COUNT"
# author: "Leonard Patterson"
# created: "2025-12-12"
# description: 
# ---

# Clear environment
rm(list = ls())  # Removes all objects from the environment

# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(tibble)
  library(viridis)
  library(DHARMa)
  library(MuMIn)
  library(performance)
  library(sf)
  })


# Load data
dat1_LA <- read.csv("Output/Tabular Data/mean_count_all_years_LA.csv")

# Extract site, treatment, plot from location
dat2_LA <- dat1_LA %>%
  separate_wider_delim(location, delim = "-", names = c("block", "treatment", "plot")) %>%
  mutate(
    # Make Year a factor for BACI analysis
    year = factor(year, levels = sort(unique(year)))
  )

# Create Site_ID and treatment_applied
pre_year <- min(as.numeric(as.character(dat2_LA$year)), na.rm = TRUE)

dat3_LA <- dat2_LA %>%
  mutate(
    # Unique experimental unit (site × assigned treatment)
    site = paste(block, treatment, sep = "-"),
    
    # BACI-ready treatment_applied: set pre-treatment year to NT
    treatment_applied = ifelse(as.numeric(as.character(year)) == pre_year, "NT", treatment),
    treatment_applied = factor(treatment_applied, levels = c("NT","LR","HR","FR"))
  )

# Reorganize dat3
spp <- names(dat3_LA)[grepl("^[A-Z]{4}$", names(dat3_LA))]
dat4_LA <- dat3_LA %>%
  select(block, site, treatment, year, all_of(spp))

# Create the new 'time_period' factor 
dat4_pooled_LA <- dat4_LA %>%
  mutate(
    time_period = case_when(
      year == 2023 ~ "Before",
      year %in% c(2024, 2025) ~ "After",
      TRUE ~ NA_character_
    ),
    # Ensure the factor order is correct: Before (reference) then After
    time_period = factor(time_period, levels = c("Before", "After")),
    
    # Ensure other factors remain correctly defined (NT must be the reference level)
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")), 
    # Define site and block as factors for the dual random effects
    site = factor(site), 
    block = factor(block), 
    year = factor(year)
  )

# Some sites are spread over different two blocks so block-site structure doesn't currently reflect
# nested structure. Hard code correct block-site structure

dat4_pooled_LA <- dat4_pooled_LA %>%
  mutate(
    block = case_when(
      # still update all MUS113 blocks
      block == "MUS113" ~ "MUS124",
      # only this specific site gets moved to MUS061A
      block == "MUS124" & site == "MUS124-NT" ~ "MUS061A",
      TRUE ~ block
    )
  )


# Save
write.csv(dat4_pooled_LA, "Output/Tabular Data/mean_count_all_years_LA_SSM.csv")








### VEG DATA

veg_2023 <- read.csv("Input/Veg data/IBWG Brushing Trial Pre-Treatment Data.csv") %>%
  rename(block = Block) %>%
  rename(treatment = Treatment) %>%
  mutate(site = paste(block, treatment, sep = "-")) 

# Rename
veg_2023_1 <- veg_2023 %>%
  select(block,treatment,site, Conifer.Density,Deciduous.Density,Slope,Moisture,Alder.Cover....,Willow.Cover....,Saskatoon.Cover....,Mountain.Ash.Cover....,High.Bush.Cranberry.Cover....,Red.Osier.Dog.Cover....) %>%
  rename(
    Alder_Cover = Alder.Cover....,
    Willow_Cover = Willow.Cover....,
    Saskatoon_Cover = Saskatoon.Cover....,
    MountainAsh_Cover = Mountain.Ash.Cover....,
    Cranberry_Cover = High.Bush.Cranberry.Cover....,
    Dogwood_Cover = Red.Osier.Dog.Cover....
  )

# Remove NAs
veg_2023_1 <- veg_2023_1 %>%
  mutate(across(ends_with("_Cover"), ~ {
    # 1. Force the column to character so it can "see" the "-"
    val <- as.character(.x)
    # 2. Replace "-" with NA
    val <- na_if(val, "-")
    # 3. Convert everything to numeric (NA stays NA, numbers become numbers)
    val <- as.numeric(val)
    # 4. Finally, turn all NAs into 0
    replace_na(val, 0)
  }))

# Summarize by site
veg_2023_2 <- veg_2023_1 %>%
  group_by(site) %>%
  summarise(
    Conifer.Density   = mean(Conifer.Density, na.rm = TRUE),
    Deciduous.Density = mean(Deciduous.Density, na.rm = TRUE),
    Slope = mean(Slope, na.rm = TRUE),
    Alder_Cover = sum(Alder_Cover, na.rm = TRUE),
    Willow_Cover = sum(Willow_Cover, na.rm = TRUE),
    Saskatoon_Cover = sum(Saskatoon_Cover, na.rm = TRUE),
    MountainAsh_Cover = sum(MountainAsh_Cover, na.rm = TRUE),
    Cranberry_Cover = sum(Cranberry_Cover, na.rm = TRUE),
    Dogwood_Cover = sum(Dogwood_Cover, na.rm = TRUE),
    .groups = "drop"
  )

# Sum shrubs
veg_2023_3 <- veg_2023_2 %>%
  mutate(Total_Shrub_Cover = rowSums(pick(ends_with("_Cover")), na.rm = TRUE))

# Extract site vars
microsite <- veg_2023 %>%
  filter(Plot == "C") %>%
  select(site,Moisture, Aspect,Position)

# merge veg and microsite
veg_2023_4 <- veg_2023_3 %>%
  left_join(microsite, by = "site")

# Merge w/ full site data
dat4_pooled_veg_1 <- dat4_pooled_LA %>%
  left_join(veg_2023_4, by = "site")


# Ensure levels are set and scale decid
dat4_pooled_veg <- dat4_pooled_veg_1 %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After")),
    moisture = factor(Moisture, levels = c("Mesic", "Hydric", "Subxeric", "SubHygric")),
    Decid_pre_sc = as.numeric(scale(Deciduous.Density)),
    Conif_pre_sc = as.numeric(scale(Conifer.Density)),
    Shrub_Cover_Scaled = as.numeric(scale(Total_Shrub_Cover)),
    Willow_Cover_Scaled = as.numeric(scale(Willow_Cover)),
    Alder_Cover_Scaled = as.numeric(scale(Alder_Cover)),
    # This must stay inside mutate to work!
    Slope = as.numeric(Slope)
  ) %>%
  select(-c(Deciduous.Density, Conifer.Density, Willow_Cover, Alder_Cover, Saskatoon_Cover, 
            MountainAsh_Cover, Cranberry_Cover, Dogwood_Cover))






### Calculate distance to wetlands
# Read in shapefiles
wetlands <- st_read("Input/Site shapefile/wetlands.shp")
sites <- st_read("Input/Site shapefile/sites.shp")

# Calculate distance from each site to the nearest wetland
# st_distance returns a matrix; we take the minimum distance for each row
dist_matrix <- st_distance(sites, wetlands)
sites$dist_to_wetland <- apply(dist_matrix, 1, min)
sites$dist_to_wetland <- as.numeric(sites$dist_to_wetland)
sites_df <- st_drop_geometry(sites) %>%
  rename(site = Plot_ID) %>%
  # Remove "-C" from the end of the site names
  mutate(site = str_remove(site, "-C$"))

# Merge the distance data into dat4_pooled_veg
# We use left_join to keep all records in your vegetation data
dat4_pooled_veg <- dat4_pooled_veg %>%
  left_join(sites_df, by = "site")

# Center and scale distance to wetland
dat4_pooled_veg$dist_scaled <- scale(dat4_pooled_veg$dist_to_wetland)




### Calculate dist to forest
dist_to_forest <- read.csv("Output/dist_to_forest.csv") %>%
  select(-X) %>%
  rename(distance_to_forest = X.1) 
dist_to_forest <- dist_to_forest[1:24,]
dist_to_forest$dist_for_scaled <- scale(dist_to_forest$distance_to_forest)
#Join
dat4_pooled_veg <- dat4_pooled_veg %>%
  left_join(dist_to_forest, by = "site")





### Plot treatments by covariates ###


## Pre-treatment subalpine fir
p1 <- ggplot(veg_2023, aes(x = treatment, y = BB.sph)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment subalpine fir (stems/ha)") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1)
  )
out_fig <- "Figures/Appendix/SubAlpineFir.tiff"
ggsave(
  filename    = out_fig,
  plot        = p1,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")



## Pre-treatment pine
p2 <- ggplot(veg_2023, aes(x = treatment, y = Pli.sph)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment lodgepole pine (stems/ha)") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1))
    
out_fig <- "Figures/Appendix/LodgepolePine.tiff"
ggsave(
  filename    = out_fig,
  plot        = p2,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")



## Pre-treatment spruce
p3 <- ggplot(veg_2023, aes(x = treatment, y = Sx.sph)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment hydrid spruce (stems/ha)") + 
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1))
out_fig <- "Figures/Appendix/Spruce.tiff"
ggsave(
  filename    = out_fig,
  plot        = p3,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")
                                                            


## Pre-treatment aspen
p4 <- ggplot(veg_2023, aes(x = treatment, y = At.sph)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 30000, by = 5000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment trembling aspen (stems/ha)") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1))
out_fig <- "Figures/Appendix/Aspen.tiff"
ggsave(
  filename    = out_fig,
  plot        = p4,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")


## Pre-treatment decid density
p5 <- ggplot(veg_2023, aes(x = treatment, y = Deciduous.Density)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 45000, by = 3000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment broadleaves (stems/ha)") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1))
out_fig <- "Figures/Appendix/Broadleaf.tiff"
ggsave(
  filename    = out_fig,
  plot        = p5,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")



## Pre-treatment conifer density
p6 <- ggplot(veg_2023, aes(x = treatment, y = Conifer.Density)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  labs(
    x = "Treatment",
    y = "Pre-treatment conifers (stems/ha)") + 
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1))
out_fig <- "Figures/Appendix/Conifers.tiff"
ggsave(
  filename    = out_fig,
  plot        = p6,
  device      = "tiff",
  width       = 18,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")




# Calculate range of aspen densities, summarizing by block
aspen_den <- veg_2023 %>%
  group_by(block) %>%
  summarise(aspen_den = mean(At.sph))
median(aspen_den$aspen_den)

fir_den <- veg_2023 %>%
  group_by(block) %>%
  summarise(fir_den = mean(BB.sph))
median(fir_den$fir_den)



# Dist to wetland ~ treatment
ggplot(dat4_pooled_veg, aes(x = treatment, y = dist_scaled)) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Dist to wetland"
  )


# Dist to forest ~ treatment
ggplot(dat4_pooled_veg, aes(x = treatment, y = dist_for_scaled)) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Distance to forest"
  )


boxplot(veg_2023$At.sph ~ treatment)

?boxplot

### Plot showing treatments by moisture

# Plot treatments by moisture regime
# 1. Filter and Prepare Data
moisture_plot <- veg_2023 %>%
  filter(Plot == "C") %>%
  select(treatment, Shrub_Cover_Scaled)

# 2. Define the Order (The "Factor" step)
# This tells R exactly how to sequence the data and the legend
moisture_plot$Moisture <- factor(moisture_plot$Moisture, 
                                 levels = c("Subxeric", "Mesic", "SubHygric", "Hydric"))

# 3. Re-run the table (it will now follow the factor levels)
counts <- table(moisture_plot$Moisture, moisture_plot$treatment)

# 4. Plot with specific colors
# Note: The colors match the order of the levels defined above
barplot(counts, 
        main = "Moisture by Treatment",
        xlab = "Treatment", 
        col = c("yellow", "darkgreen", "lightblue", "blue"),
        legend = rownames(counts),
        beside = TRUE,
        args.legend = list(x = "topright", bty = "n")) # Optional: cleans up legend box







### AMRE sensitivity

datafter <- dat4_pooled_veg_5 %>%
  filter(time_period == "After")

after_WTSP <- glmmTMB(WTSP ~ treatment + WTSP_baseline_centered + (1|block/site), family = compois(link = "log"), data = datafter)
summary(after_WTSP)
after_ALFL <- glmmTMB(ALFL ~ treatment + ALFL_baseline_centered + (1|block/site), family = poisson(link = "log"), data = datafter)
summary(after_ALFL)






### Create baseline covariates
# Create and Center baseline covariates
species_vars <- c("AMRE", "ALFL", "DUFL", "OCWA", "WAVI", 
                  "WTSP", "SWTH", "DEJU", "YRWA")

baseline_all <- dat4_pooled_veg %>%
  filter(time_period == "Before") %>%
  group_by(site) %>%
  summarise(
    across(
      all_of(species_vars),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_baseline"
    ),
    .groups = "drop"
  ) %>%
  # Now center the baseline values by subtracting the mean of all sites
  mutate(
    across(
      ends_with("_baseline"),
      ~ .x - mean(.x, na.rm = TRUE),
      .names = "{.col}_centered"
    )
  )

# Join back to full dataset
dat4_pooled_veg_5 <- dat4_pooled_veg %>%
  left_join(baseline_all, by = "site")


AMRE1 <- glmmTMB(AMRE ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE2 <- glmmTMB(AMRE ~ treatment * time_period  + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE3 <- glmmTMB(AMRE ~ treatment * time_period + AMRE_baseline_centered + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE4 <- glmmTMB(AMRE ~ treatment * time_period + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AICc(AMRE1)
AICc(AMRE2)
AICc(AMRE3)
AICc(AMRE4)
check_collinearity(AMRE1)
check_collinearity(AMRE2)
check_collinearity(AMRE3)
check_collinearity(AMRE4)

### AMRE MODEL SELECTION

AMRE <- glmmTMB(AMRE ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE1 <- glmmTMB(AMRE ~ treatment * time_period + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE2 <- glmmTMB(AMRE ~ treatment * time_period + moisture + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE3 <- glmmTMB(AMRE ~ treatment * time_period + moisture  + AMRE_baseline + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE4 <- glmmTMB(AMRE ~ treatment * time_period + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE5 <- glmmTMB(AMRE ~ treatment * time_period + Aspect + moisture + Shrub_Cover_Scaled + I(Shrub_Cover_Scaled^2) + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE6 <- glmmTMB(AMRE ~ treatment * time_period + moisture + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
AMRE7 <- glmmTMB(AMRE ~ treatment * time_period + AMRE_baseline + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)


AICc(AMRE)
AICc(AMRE1)
AICc(AMRE2)
AICc(AMRE3) # best
AICc(AMRE4) 
AICc(AMRE5)
AICc(AMRE6)
AICc(AMRE7)


summary(AMRE)
summary(AMRE1)
summary(AMRE2)
summary(AMRE3)
summary(AMRE4)
summary(AMRE5)
summary(AMRE7)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_AMRE1 <- simulateResiduals(
  fittedModel = AMRE3,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_AMRE1)


### ALFL MODEL SELECTION

ALFL <- glmmTMB(ALFL ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg) # base
ALFL1 <- glmmTMB(ALFL ~ treatment * time_period + dist_scaled + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # moisture
ALFL2 <- glmmTMB(ALFL ~ treatment * time_period + Decid_pre_sc + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # pre decid
ALFL3 <- glmmTMB(ALFL ~ treatment * time_period + Decid_pre_sc + I(Decid_pre_sc^2) + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # pre decid ^2
ALFL4 <- glmmTMB(ALFL ~ treatment * time_period + moisture + Decid_pre_sc + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # moist and pre decid
ALFL5 <- glmmTMB(ALFL ~ treatment * time_period + ALFL_baseline_centered + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # moist and pre decid ^2
ALFL6 <- glmmTMB(ALFL ~ treatment * time_period + Aspect + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # aspect only
ALFL7 <- glmmTMB(ALFL ~ treatment * time_period + Shrub_Cover_Scaled + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # shrub cover only
ALFL8 <- glmmTMB(ALFL ~ treatment * time_period + moisture + Shrub_Cover_Scaled + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # moisture and shrub cover
ALFL9 <- glmmTMB(ALFL ~ treatment * time_period + moisture + Shrub_Cover_Scaled + I(Shrub_Cover_Scaled^2) + (1|block/site), family = poisson(link = "log"), data = dat4_pooled_veg) # moisture and shrub cover

AICc(ALFL) # Best
AICc(ALFL1)
AICc(ALFL2)
AICc(ALFL3)
AICc(ALFL4) 
AICc(ALFL6)
AICc(ALFL7)
AICc(ALFL8)
AICc(ALFL9)

summary(ALFL)



### DUFL

DUFL <- glmmTMB(DUFL ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL1 <- glmmTMB(DUFL ~ treatment * time_period + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL2 <- glmmTMB(DUFL ~ treatment * time_period + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL3 <- glmmTMB(DUFL ~ treatment * time_period + moisture + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL4 <- glmmTMB(DUFL ~ treatment * time_period + moisture + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL5 <- glmmTMB(DUFL ~ treatment * time_period + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL6 <- glmmTMB(DUFL ~ treatment * time_period + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL7 <- glmmTMB(DUFL ~ treatment * time_period + moisture + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
DUFL8 <- glmmTMB(DUFL ~ treatment * time_period + DUFL_baseline_centered + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)

AICc(DUFL)
AICc(DUFL1)
AICc(DUFL2) # best
AICc(DUFL3) 
AICc(DUFL4) 
AICc(DUFL5)
AICc(DUFL6)
AICc(DUFL7)
AICc(DUFL8)

summary(DUFL2)



### OCWA

OCWA <- glmmTMB(OCWA ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA1 <- glmmTMB(OCWA ~ treatment * time_period + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA2 <- glmmTMB(OCWA ~ treatment * time_period + Decid_pre_sc + I(Decid_pre_sc^2) + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA3 <- glmmTMB(OCWA ~ treatment * time_period + moisture + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA4 <- glmmTMB(OCWA ~ treatment * time_period + moisture + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA5 <- glmmTMB(OCWA ~ treatment * time_period + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA6 <- glmmTMB(OCWA ~ treatment * time_period + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA7 <- glmmTMB(OCWA ~ treatment * time_period + moisture + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
OCWA8 <- glmmTMB(DUFL ~ treatment * time_period + OCWA_baseline_centered + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)

AICc(OCWA)
AICc(OCWA1)
AICc(OCWA2) 
AICc(OCWA3) # Best
AICc(OCWA4) 
AICc(OCWA5)
AICc(OCWA6)
AICc(OCWA7)
AICc(OCWA8)





### WAVI

WAVI <- glmmTMB(WAVI ~ treatment * time_period + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI1 <- glmmTMB(WAVI ~ treatment * time_period + moisture + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI2 <- glmmTMB(WAVI ~ treatment * time_period + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI3 <- glmmTMB(WAVI ~ treatment * time_period + moisture + Decid_pre_sc + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI4 <- glmmTMB(WAVI ~ treatment * time_period + moisture + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI5 <- glmmTMB(WAVI ~ treatment * time_period + Aspect + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI6 <- glmmTMB(WAVI ~ treatment * time_period + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)
WAVI7 <- glmmTMB(WAVI ~ treatment * time_period + moisture + Shrub_Cover_Scaled + (1|block/site), family = compois(link = "log"), data = dat4_pooled_veg)

AICc(WAVI)
AICc(WAVI1)
AICc(WAVI2) 
AICc(WAVI3) # Best
AICc(WAVI4) 
AICc(WAVI5)
AICc(WAVI6)
AICc(WAVI7)

summary(OCWA3)





### AMRE 

##EMMs

# EMMs averaged over moisture classes
emm_AMRE3 <- emmeans(
  AMRE1,
  ~ treatment * time_period,
  type = "link"
)

# View table
emm_AMRE3_df <- as.data.frame(emm_AMRE3)
print(emm_AMRE3_df)

# Plot
p_AMRE3_emm <- ggplot(
  emm_AMRE3_df,
  aes(x = treatment, y = response, color = time_period, group = time_period)
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "AMRE estimated marginal means",
    x = "Treatment",
    y = "Estimated mean abundance",
    color = "Time period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_AMRE3_emm)



## AMRE BACI

# -----------------------------
# 1. EMMs on the link scale
# -----------------------------
emm_AMRE3_link <- emmeans(
  AMRE3,
  ~ treatment * time_period,
  type = "link"
)

# View link-scale EMMs if wanted
emm_AMRE3_link_df <- as.data.frame(emm_AMRE3_link)
print(emm_AMRE3_link_df)

# -----------------------------
# 2. BACI contrasts:
#    (After - Before)_treatment - (After - Before)_NT
# -----------------------------
baci_AMRE3 <- contrast(
  emm_AMRE3_link,
  method = list(
    LR_vs_NT = c( 1, -1,  0,  0, -1,  1,  0,  0),
    HR_vs_NT = c( 1,  0, -1,  0, -1,  0,  1,  0),
    FR_vs_NT = c( 1,  0,  0, -1, -1,  0,  0,  1)
  )
)

# Summary on link scale
baci_AMRE3_link_df <- as.data.frame(summary(baci_AMRE3, infer = c(TRUE, TRUE)))
print(baci_AMRE3_link_df)

# -----------------------------
# 3. Back-transform to response scale
#    Because model uses log link, exp(contrast) gives
#    ratio of change relative to NT
# -----------------------------
baci_AMRE3_resp_df <- baci_AMRE3_link_df %>%
  mutate(
    ratio = exp(estimate),
    lower = exp(asymp.LCL),
    upper = exp(asymp.UCL)
  )

print(baci_AMRE3_resp_df)

# -----------------------------
# 4. Plot BACI contrasts
#    Reference line at 1 = no BACI effect relative to NT
# -----------------------------
p_baci_AMRE3 <- ggplot(
  baci_AMRE3_resp_df,
  aes(x = contrast, y = ratio)
) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.8
  ) +
  labs(
    title = "AMRE BACI contrasts",
    subtitle = "Ratio of temporal change relative to NT",
    x = "Treatment contrast",
    y = "BACI ratio (After/Before relative to NT)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_baci_AMRE3)

# -----------------------------
# 5. Optional: nicer labels
# -----------------------------
baci_AMRE3_resp_df <- baci_AMRE3_resp_df %>%
  mutate(
    contrast = recode(
      contrast,
      LR_vs_NT = "LR vs NT",
      HR_vs_NT = "HR vs NT",
      FR_vs_NT = "FR vs NT"
    )
  )

p_baci_AMRE3 <- ggplot(
  baci_AMRE3_resp_df,
  aes(x = contrast, y = ratio)
) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.8
  ) +
  labs(
    title = "AMRE BACI contrasts",
    subtitle = "Ratio of temporal change relative to NT",
    x = NULL,
    y = "BACI ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_baci_AMRE3)











### ALFL EMM

# EMMs averaged over moisture classes
emm_ALFL <- emmeans(
  ALFL,
  ~ treatment * time_period,
  type = "link"
)

# View table
emm_ALFL_df <- as.data.frame(emm_ALFL)
print(emm_ALFL_df)

# Plot
p_ALFL_emm <- ggplot(
  emm_ALFL_df,
  aes(x = treatment, y = rate, color = time_period, group = time_period)
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "ALFL estimated marginal means",
    x = "Treatment",
    y = "Estimated mean abundance",
    color = "Time period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_ALFL_emm)





### DUFL EMM

# EMMs averaged over moisture classes
emm_DUFL2 <- emmeans(
  DUFL,
  ~ treatment * time_period,
  type = "response"
)

# View table
emm_DUFL2_df <- as.data.frame(emm_DUFL2)
print(emm_DUFL2_df)

# Plot
p_DUFL2_emm <- ggplot(
  emm_DUFL2_df,
  aes(x = treatment, y = response, color = time_period, group = time_period)
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "DUFL estimated marginal means",
    x = "Treatment",
    y = "Estimated mean abundance",
    color = "Time period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_DUFL2_emm)






### OCWA EMM

# EMMs averaged over moisture classes
emm_OCWA3 <- emmeans(
  WTSP_model_LA,
  ~ treatment * time_period,
  type = "response"
)

# View table
emm_OCWA3_df <- as.data.frame(emm_OCWA3)
print(emm_OCWA3_df)

# Plot
p_OCWA3_emm <- ggplot(
  emm_OCWA3_df,
  aes(x = treatment, y = response, color = time_period, group = time_period)
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "OCWA estimated marginal means",
    subtitle = "Averaged over moisture classes and pre- deciduous stems",
    x = "Treatment",
    y = "Estimated mean abundance",
    color = "Time period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_OCWA3_emm)






### WAVI EMM

# EMMs averaged over moisture classes
emm_WAVI3 <- emmeans(
  after_WTSP,
  ~ treatment,
  type = "response"
)

# View table
emm_WAVI3_df <- as.data.frame(emm_WAVI3)
print(emm_WAVI3_df)

# Plot
p_WAVI3_emm <- ggplot(
  emm_WAVI3_df,
  aes(x = treatment, y = response)
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "WAVI estimated marginal means",
    subtitle = "Averaged over moisture classes and pre- decid stems",
    x = "Treatment",
    y = "Estimated mean abundance",
    color = "Time period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_WAVI3_emm)

















####################### NEW





##############################
## BACI + PRETREATMENT STEMS
## Limited Amplitude dataset
##############################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(glmmTMB)
  library(emmeans)
  library(DHARMa)
  library(MuMIn)
  library(cowplot)
  library(stringr)
})

## =========================================================
## 1. DATA PREP
## =========================================================

dat4_pooled_veg <- dat4_pooled_veg %>%
  mutate(
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After")),
    Decid_pre_sc = as.numeric(scale(Deciduous.Density))
  )

# Optional: inspect scaling
decid_mean <- mean(dat4_pooled_veg$Deciduous.Density, na.rm = TRUE)
decid_sd   <- sd(dat4_pooled_veg$Deciduous.Density, na.rm = TRUE)

cat("\nPre-treatment deciduous density summary:\n")
cat("Mean =", round(decid_mean, 2), "\n")
cat("SD   =", round(decid_sd, 2), "\n")
cat("Low  = Mean - 1 SD =", round(decid_mean - decid_sd, 2), "\n")
cat("High = Mean + 1 SD =", round(decid_mean + decid_sd, 2), "\n")

# Define stem-density levels for conditional BACI / EMM estimation
stem_levels <- c(Low = -1, Mean = 0, High = 1)

## =========================================================
## 2. SPECIES / FAMILY SETUP
## =========================================================

species_families <- tribble(
  ~Species, ~Family,
  "AMRE", "compois",
  "ALFL", "poisson",
  "DUFL", "compois",
  "OCWA", "compois",
  "SWTH", "poisson",
  "WAVI", "compois",
  "WTSP", "poisson"
)

species_lookup <- c(
  AMRE = "American Redstart",
  ALFL = "Alder Flycatcher",
  DUFL = "Dusky Flycatcher",
  OCWA = "Orange-crowned Warbler",
  SWTH = "Swainson's Thrush",
  WAVI = "Warbling Vireo",
  WTSP = "White-throated Sparrow"
)

## =========================================================
## 3. HELPER FUNCTIONS
## =========================================================

get_family_object <- function(fam_name) {
  switch(
    fam_name,
    poisson = poisson(link = "log"),
    compois = compois(link = "log"),
    stop("Unsupported family: ", fam_name)
  )
}

fit_species_model <- function(species_name, fam_name, data) {
  fm <- as.formula(
    paste0(species_name, " ~ treatment * time_period * Decid_pre_sc + (1|block/site)")
  )
  
  glmmTMB(
    formula = fm,
    family  = get_family_object(fam_name),
    data    = data
  )
}

standardize_ci_cols <- function(tbl) {
  if (all(c("lower.CL", "upper.CL") %in% names(tbl))) {
    tbl <- tbl %>% rename(LCL = lower.CL, UCL = upper.CL)
  } else if (all(c(".lower", ".upper") %in% names(tbl))) {
    tbl <- tbl %>% rename(LCL = .lower, UCL = .upper)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(tbl))) {
    tbl <- tbl %>% rename(LCL = asymp.LCL, UCL = asymp.UCL)
  } else if (!all(c("LCL", "UCL") %in% names(tbl))) {
    stop("Could not find CI columns. Names were: ", paste(names(tbl), collapse = ", "))
  }
  tbl
}

calculate_baci_log_contrast <- function(model, species_name,
                                        cov_value,
                                        cov_name = "Decid_pre_sc",
                                        cov_label = "Mean") {
  
  emm <- emmeans(
    model,
    ~ treatment * time_period,
    type = "link",
    at = setNames(list(cov_value), cov_name)
  )
  
  # Row order:
  # 1 NT Before
  # 2 LR Before
  # 3 HR Before
  # 4 FR Before
  # 5 NT After
  # 6 LR After
  # 7 HR After
  # 8 FR After
  
  baci <- contrast(
    emm,
    method = list(
      LR = c( 1, -1,  0,  0, -1,  1,  0,  0),
      HR = c( 1,  0, -1,  0, -1,  0,  1,  0),
      FR = c( 1,  0,  0, -1, -1,  0,  0,  1)
    ),
    adjust = "none"
  )
  
  sum_tbl <- as.data.frame(summary(baci))
  ci_tbl  <- as.data.frame(confint(baci)) %>% standardize_ci_cols()
  
  out <- left_join(
    sum_tbl[, c("contrast", "estimate", "p.value")],
    ci_tbl[, c("contrast", "LCL", "UCL")],
    by = "contrast"
  ) %>%
    transmute(
      Species     = species_name,
      StemLevel   = cov_label,
      Treatment   = contrast,
      LogContrast = estimate,
      LCL         = LCL,
      UCL         = UCL,
      PValue      = p.value
    ) %>%
    arrange(Species, StemLevel, Treatment)
  
  out
}

get_log_emm_table <- function(model, species_name,
                              cov_value,
                              cov_name = "Decid_pre_sc",
                              cov_label = "Mean") {
  
  emm <- emmeans(
    model,
    ~ treatment * time_period,
    type = "link",
    at = setNames(list(cov_value), cov_name)
  )
  
  sum_tbl <- summary(emm) %>%
    as_tibble() %>%
    select(treatment, time_period, emmean, SE)
  
  ci_tbl <- confint(emm) %>%
    as_tibble() %>%
    standardize_ci_cols() %>%
    select(treatment, time_period, LCL, UCL)
  
  out <- sum_tbl %>%
    left_join(ci_tbl, by = c("treatment", "time_period")) %>%
    transmute(
      Species     = species_name,
      StemLevel   = cov_label,
      treatment   = treatment,
      time_period = time_period,
      LogMean     = emmean,
      SE          = SE,
      LCL         = LCL,
      UCL         = UCL
    ) %>%
    arrange(Species, StemLevel, treatment, time_period)
  
  out
}

get_model_summary_row <- function(model, species_name, fam_name) {
  fe <- summary(model)$coefficients$cond
  fe_names <- rownames(fe)
  
  three_way_term <- fe_names[str_detect(fe_names, "treatment.*:time_period.*:Decid_pre_sc|treatment.*:Decid_pre_sc.*:time_period")]
  
  if (length(three_way_term) == 0) {
    three_way_term <- NA_character_
    three_way_est  <- NA_real_
    three_way_p    <- NA_real_
  } else {
    # keep all 3-way terms, not just one
    three_way_est <- paste(round(fe[three_way_term, "Estimate"], 3), collapse = "; ")
    three_way_p   <- paste(signif(fe[three_way_term, "Pr(>|z|)"], 3), collapse = "; ")
    three_way_term <- paste(three_way_term, collapse = "; ")
  }
  
  tibble(
    Species = species_name,
    Family  = fam_name,
    AICc    = AICc(model),
    ThreeWayTerms = three_way_term,
    ThreeWayEstimates = three_way_est,
    ThreeWayPValues   = three_way_p
  )
}

## =========================================================
## 4. FIT ALL MODELS
## =========================================================

model_list_LA <- list()
model_summary_tbl_LA <- tibble()
residual_list_LA <- list()

for (i in seq_len(nrow(species_families))) {
  sp  <- species_families$Species[i]
  fam <- species_families$Family[i]
  
  cat("\n=================================================\n")
  cat("Fitting model for:", sp, "| Family:", fam, "\n")
  cat("=================================================\n")
  
  mod <- fit_species_model(sp, fam, dat4_pooled_veg)
  
  model_list_LA[[sp]] <- mod
  model_summary_tbl_LA <- bind_rows(
    model_summary_tbl_LA,
    get_model_summary_row(mod, sp, fam)
  )
  
  set.seed(123)
  residual_list_LA[[sp]] <- simulateResiduals(
    fittedModel = mod,
    n = 1000
  )
  
  print(summary(mod))
  print(AICc(mod))
  plot(residual_list_LA[[sp]])
}

## Save model summary
write.csv(
  model_summary_tbl_LA,
  "Output/Tables/SSM_Model_Summary_with_DecidDensity_LA.csv",
  row.names = FALSE
)

## =========================================================
## 5. BACI CONTRASTS AT LOW / MEAN / HIGH STEM DENSITY
## =========================================================

all_contrasts_df_LA <- purrr::imap_dfr(
  model_list_LA,
  ~ bind_rows(
    calculate_baci_log_contrast(.x, .y, cov_value = stem_levels["Low"],  cov_label = "Low"),
    calculate_baci_log_contrast(.x, .y, cov_value = stem_levels["Mean"], cov_label = "Mean"),
    calculate_baci_log_contrast(.x, .y, cov_value = stem_levels["High"], cov_label = "High")
  )
)

print(
  all_contrasts_df_LA %>%
    mutate(across(c(LogContrast, LCL, UCL), ~ round(.x, 3))) %>%
    arrange(Species, StemLevel, Treatment)
)

write.csv(
  all_contrasts_df_LA,
  "Output/Tables/SSM_BACI_LogContrasts_LA_byStemLevel.csv",
  row.names = FALSE
)

## Convert log contrasts to IRRs
all_contrasts_IRR_LA <- all_contrasts_df_LA %>%
  mutate(
    IRR       = exp(LogContrast),
    IRR_LCL   = exp(LCL),
    IRR_UCL   = exp(UCL),
    Delta     = IRR - 1,
    Delta_pct = 100 * (IRR - 1)
  )

species_IRR_table_LA <- all_contrasts_IRR_LA %>%
  mutate(
    LogCI_half = (UCL - LCL) / 2,
    LogContrast = round(LogContrast, 3),
    LogCI_half  = round(LogCI_half, 3),
    LCL         = round(LCL, 3),
    UCL         = round(UCL, 3),
    IRR         = round(IRR, 2),
    IRR_LCL     = round(IRR_LCL, 2),
    IRR_UCL     = round(IRR_UCL, 2),
    Delta_pct   = round(Delta_pct, 1),
    Delta_pct_text = case_when(
      Delta_pct > 0 ~ paste0(Delta_pct, "% increase"),
      Delta_pct < 0 ~ paste0(abs(Delta_pct), "% decrease"),
      TRUE          ~ "0.0% (no change)"
    ),
    Interpretation = case_when(
      IRR_LCL > 1 ~ "Increase (CI > 1)",
      IRR_UCL < 1 ~ "Decrease (CI < 1)",
      TRUE        ~ "NS (includes 1)"
    )
  ) %>%
  arrange(Species, StemLevel, Treatment) %>%
  select(
    Species, StemLevel, Treatment,
    LogContrast, LogCI_half, LCL, UCL,
    IRR, IRR_LCL, IRR_UCL,
    Delta_pct, Delta_pct_text,
    PValue, Interpretation
  )

print(species_IRR_table_LA)

write.csv(
  species_IRR_table_LA,
  "Output/Tables/SSM_BACI_IRR_LA_byStemLevel.csv",
  row.names = FALSE
)

## =========================================================
## 6. EMMs AT LOW / MEAN / HIGH STEM DENSITY
## =========================================================

all_emmeans_log_LA <- purrr::imap_dfr(
  model_list_LA,
  ~ bind_rows(
    get_log_emm_table(.x, .y, cov_value = stem_levels["Low"],  cov_label = "Low"),
    get_log_emm_table(.x, .y, cov_value = stem_levels["Mean"], cov_label = "Mean"),
    get_log_emm_table(.x, .y, cov_value = stem_levels["High"], cov_label = "High")
  )
)

print(
  all_emmeans_log_LA %>%
    mutate(across(c(LogMean, SE, LCL, UCL), ~ round(.x, 3))) %>%
    arrange(Species, StemLevel, treatment, time_period)
)

write.csv(
  all_emmeans_log_LA,
  "Output/Tables/SSM_EMMeans_LogScale_LA_byStemLevel.csv",
  row.names = FALSE
)

## Response-scale EMMs
all_emmeans_resp_LA <- all_emmeans_log_LA %>%
  mutate(
    Species     = unname(species_lookup[as.character(Species)]),
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After")),
    StemLevel   = factor(StemLevel, levels = c("Low", "Mean", "High")),
    Mean        = exp(LogMean),
    LCL_resp    = exp(LCL),
    UCL_resp    = exp(UCL)
  )

all_emmeans_resp_LA <- all_emmeans_resp_LA %>%
  mutate(
    Species = factor(Species),
    Species = forcats::fct_rev(forcats::fct_relevel(Species, sort(levels(Species))))
  )

## =========================================================
## 7. PLOT CONDITIONAL BACI CONTRASTS
## =========================================================

# Exclusions carried over from your old script
excluded_species <- c("OSFL", "RCKI", "LISP", "AMRO")

all_contrasts_df_LA_filtered <- all_contrasts_df_LA %>%
  filter(!Species %in% excluded_species) %>%
  filter(!is.na(LogContrast)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR")),
    StemLevel = factor(StemLevel, levels = c("Low", "Mean", "High")),
    Species   = factor(Species),
    Species   = forcats::fct_rev(forcats::fct_relevel(Species, sort(levels(Species))))
  )

contrast_plot_log_treat_LA <- ggplot(all_contrasts_df_LA_filtered,
                                     aes(x = LogContrast, y = Species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 height = 0.2, linewidth = 0.9, color = "black") +
  geom_point(aes(color = Treatment), size = 2.8) +
  facet_grid(StemLevel ~ Treatment, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +
  labs(
    title = "Limited Amplitude",
    subtitle = "Conditional BACI contrasts at low, mean, and high pre-treatment deciduous density",
    x = "BACI contrast (log scale)",
    y = "Species"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y      = element_text(face = "italic"),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(contrast_plot_log_treat_LA)

cowplot::save_plot(
  "Figures/Species-level models/BACI_SSM_LA_logcoef_byStemLevel.png",
  contrast_plot_log_treat_LA,
  base_width = 10,
  base_height = 9
)

## =========================================================
## 8. PLOT RESPONSE-SCALE EMMs BY STEM LEVEL
## =========================================================

all_emmeans_resp_plot_LA <- all_emmeans_resp_LA %>%
  filter(!Species %in% unname(species_lookup[excluded_species])) %>%
  filter(!is.na(Mean))

p_ssm_emm_LA <- ggplot(
  all_emmeans_resp_plot_LA,
  aes(x = treatment, y = Mean, color = time_period, group = time_period)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.3
  ) +
  geom_errorbar(
    aes(ymin = LCL_resp, ymax = UCL_resp),
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  facet_grid(StemLevel ~ Species, scales = "free_y") +
  labs(
    title = "Limited Amplitude",
    subtitle = "Estimated marginal mean abundance at low, mean, and high pre-treatment deciduous density",
    x     = "Treatment",
    y     = "Estimated marginal mean abundance",
    color = "Time period"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5),
    strip.text       = element_text(face = "italic", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "right",
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(p_ssm_emm_LA)

cowplot::save_plot(
  "Figures/Species-level models/SSM_LA_marginal_means_byStemLevel.png",
  p_ssm_emm_LA,
  base_width = 18,
  base_height = 10
)

## =========================================================
## 9. OPTIONAL: PLOT LOG-SCALE EMMs BY STEM LEVEL
## =========================================================

all_emmeans_log_plot_LA <- all_emmeans_log_LA %>%
  filter(!Species %in% excluded_species) %>%
  mutate(
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After")),
    StemLevel   = factor(StemLevel, levels = c("Low", "Mean", "High"))
  )

emm_forest_plot_LA <- ggplot(
  all_emmeans_log_plot_LA,
  aes(x = treatment, y = LogMean, color = time_period)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_errorbar(
    aes(ymin = LCL, ymax = UCL),
    position = position_dodge(width = 0.5),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2.3
  ) +
  facet_grid(StemLevel ~ Species, scales = "free_y") +
  labs(
    title = "Limited Amplitude",
    subtitle = "Estimated log mean abundance at low, mean, and high pre-treatment deciduous density",
    x = "Treatment",
    y = "Estimated log mean count (± 95% CI)",
    color = "Time period"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "italic", size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border     = element_rect(color = "gray90", fill = NA),
    legend.position  = "bottom"
  )

print(emm_forest_plot_LA)

cowplot::save_plot(
  "Figures/Species-level models/EMM_ForestPlot_LA_byStemLevel.png",
  emm_forest_plot_LA,
  base_width = 18,
  base_height = 10
)

## =========================================================
## 10. QUICK INTERPRETATION TABLE FOR MEAN STEM DENSITY ONLY
## =========================================================

species_IRR_table_LA_mean <- species_IRR_table_LA %>%
  filter(StemLevel == "Mean")

write.csv(
  species_IRR_table_LA_mean,
  "Output/Tables/SSM_BACI_IRR_LA_MeanStemDensityOnly.csv",
  row.names = FALSE
)

cat("\n=================================================\n")
cat("DONE\n")
cat("Outputs saved:\n")
cat("- SSM_Model_Summary_with_DecidDensity_LA.csv\n")
cat("- SSM_BACI_LogContrasts_LA_byStemLevel.csv\n")
cat("- SSM_BACI_IRR_LA_byStemLevel.csv\n")
cat("- SSM_BACI_IRR_LA_MeanStemDensityOnly.csv\n")
cat("- SSM_EMMeans_LogScale_LA_byStemLevel.csv\n")
cat("- BACI_SSM_LA_logcoef_byStemLevel.png\n")
cat("- SSM_LA_marginal_means_byStemLevel.png\n")
cat("- EMM_ForestPlot_LA_byStemLevel.png\n")
cat("=================================================\n")












#Make shapefile
library(sf)

sites <- st_read("Input/Site shapefile/sites.shp")

dat4_pooled_veg_2023 <- dat4_pooled_veg %>%
  filter(year == 2023) %>%
  rename(Plot_ID = site) %>%
  mutate(Plot_ID = paste0(Plot_ID, "-C"))

dat4_pooled_veg_2024 <- dat4_pooled_veg %>%
  filter(year == 2024) %>%
  rename(Plot_ID = site) %>%
  mutate(Plot_ID = paste0(Plot_ID, "-C"))

dat4_pooled_veg_2025 <- dat4_pooled_veg %>%
  filter(year == 2025) %>%
  rename(Plot_ID = site) %>%
  mutate(Plot_ID = paste0(Plot_ID, "-C"))

# Join
sites_2023 <- sites %>%
  left_join(dat4_pooled_veg_2023, by = "Plot_ID")
sites_2024 <- sites %>%
  left_join(dat4_pooled_veg_2024, by = "Plot_ID")
sites_2025 <- sites %>%
  left_join(dat4_pooled_veg_2025, by = "Plot_ID")


st_write(sites_2023, "Input/Site shapefile/sites_2023.shp", append = FALSE)
st_write(sites_2024, "Input/Site shapefile/sites_2024.shp", append = FALSE)
st_write(sites_2025, "Input/Site shapefile/sites_2025.shp", append = FALSE)










########### Dredge

suppressPackageStartupMessages({
  library(glmmTMB)
  library(MuMIn)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# -----------------------------
# 1. Species-specific families
# -----------------------------
species_families <- tribble(
  ~Species, ~Family,
  "AMRE", "compois",
  "ALFL", "poisson",
  "DUFL", "compois",
  "OCWA", "compois",
  "SWTH", "poisson",
  "WAVI", "compois",
  "WTSP", "poisson"
)

# keep only species you want here
species_vec <- c("AMRE", "ALFL", "DUFL", "OCWA", "WTSP", "SWTH")

# -----------------------------
# 2. Prepare data
# -----------------------------
dat_dredge <- dat4_pooled_veg %>%
  mutate(
    aspect = factor(Aspect),
    slope  = Slope
  )

options(na.action = "na.fail")

cand_terms <- c(
  "Shrub_Cover_Scaled",
  "Conif_pre_sc",
  "Decid_pre_sc",
  "moisture",
  "aspect",
  "slope"
)

# -----------------------------
# 3. Helper: map family names
# -----------------------------
get_family_object <- function(fam_name) {
  switch(
    fam_name,
    "poisson" = poisson(link = "log"),
    "compois" = compois(link = "log"),
    "nbinom2" = nbinom2(link = "log"),
    stop("Unsupported family: ", fam_name)
  )
}

# -----------------------------
# 4. Containers
# -----------------------------
global_model_list <- list()
dredge_list       <- list()
best_model_list   <- list()
top_models_list   <- list()

# -----------------------------
# 5. Fit + dredge loop
# -----------------------------
for (sp in species_vec) {
  
  fam_name <- species_families %>%
    filter(Species == sp) %>%
    pull(Family)
  
  if (length(fam_name) == 0) {
    stop("No family found for species: ", sp)
  }
  
  fam_obj <- get_family_object(fam_name)
  
  rhs <- paste(c("treatment * time_period", cand_terms), collapse = " + ")
  fmla <- as.formula(paste0(sp, " ~ ", rhs, " + (1 | block/site)"))
  
  cat("\n==============================\n")
  cat("Species:", sp, "\n")
  cat("Family :", fam_name, "\n")
  cat("==============================\n")
  
  global_mod <- glmmTMB(
    formula = fmla,
    family  = fam_obj,
    data    = dat_dredge
  )
  
  global_model_list[[sp]] <- global_mod
  
  dd <- dredge(
    global_mod,
    fixed = c("treatment", "time_period", "treatment:time_period"),
    rank = "AICc",
    trace = FALSE
  )
  
  dredge_list[[sp]] <- dd
  best_model_list[[sp]] <- get.models(dd, subset = 1)[[1]]
  top_models_list[[sp]] <- subset(dd, delta < 2)
  
  cat("\nTop models (ΔAICc < 2) for", sp, ":\n")
  print(top_models_list[[sp]])
}

# -----------------------------
# 6. Extract model terms helper
# -----------------------------
extract_terms <- function(model_formula) {
  tl <- attr(terms(model_formula), "term.labels")
  tl <- tl[!grepl("\\|", tl)]
  tl
}

# -----------------------------
# 7. Build compact top-model table
# -----------------------------
compact_top_table <- bind_rows(
  lapply(names(top_models_list), function(sp) {
    
    dd_sub <- as.data.frame(top_models_list[[sp]])
    
    mods <- get.models(
      dredge_list[[sp]],
      subset = delta < 2
    )
    
    bind_rows(lapply(seq_along(mods), function(i) {
      mod <- mods[[i]]
      fml <- formula(mod)
      trms <- extract_terms(fml)
      
      tibble(
        Species = sp,
        Family = species_families %>% filter(Species == sp) %>% pull(Family),
        Model_rank = i,
        AICc = dd_sub$AICc[i],
        Delta_AICc = dd_sub$delta[i],
        Weight = dd_sub$weight[i],
        K = dd_sub$df[i],
        Formula = deparse(fml),
        Includes_Shrub = "Shrub_Cover_Scaled" %in% trms,
        Includes_Conif = "Conif_pre_sc" %in% trms,
        Includes_Decid = "Decid_pre_sc" %in% trms,
        Includes_moisture = "moisture" %in% trms,
        Includes_aspect = "aspect" %in% trms,
        Includes_slope = "slope" %in% trms
      )
    }))
  })
)

print(compact_top_table, n = Inf)

# -----------------------------
# 8. Add compact optional-term column
# -----------------------------
compact_top_table_clean <- compact_top_table %>%
  rowwise() %>%
  mutate(
    Optional_terms = paste(
      c(
        if (Includes_Shrub) "Shrub_Cover_Scaled",
        if (Includes_Conif) "Conif_pre_sc",
        if (Includes_Decid) "Decid_pre_sc",
        if (Includes_moisture) "moisture",
        if (Includes_aspect) "aspect",
        if (Includes_slope) "slope"
      ),
      collapse = ", "
    )
  ) %>%
  ungroup()

print(
  compact_top_table_clean %>%
    select(Species, Family, Model_rank, AICc, Delta_AICc, Weight, K, Optional_terms, Formula),
  n = Inf
)

# -----------------------------
# 9. Optional: best model table
# -----------------------------
best_model_summary <- compact_top_table_clean %>%
  group_by(Species) %>%
  slice_min(order_by = AICc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Species, Family, AICc, Delta_AICc, Weight, K, Optional_terms, Formula)

print(best_model_summary, n = Inf)

# -----------------------------
# 10. Optional save
# -----------------------------
# saveRDS(global_model_list, "global_model_list_species_family_dredge.rds")
# saveRDS(dredge_list, "dredge_results_species_family.rds")
# saveRDS(best_model_list, "best_model_list_species_family_dredge.rds")
# saveRDS(top_models_list, "top_models_delta_lt2_species_family.rds")
# write.csv(compact_top_table_clean, "top_models_delta_lt2_species_family.csv", row.names = FALSE)
# write.csv(best_model_summary, "best_model_summary_species_family.csv", row.names = FALSE)
