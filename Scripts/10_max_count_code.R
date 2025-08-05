############# 2. CLEAN DATA



### CALCULATE MAX COUNT FOR UNLIMITED DISTANCE PC

# Load data
dat1_UD <- read.csv("Input/Tag reports/tags_all_years.csv")

# Ensure individual_order is numeric and extract year

dat2_UD <- dat1_UD %>%
  mutate(
    individual_order = as.numeric(individual_order),
    year = year(ymd_hms(recording_date_time))
  ) %>%
  
  # Aggregate: maximum individual_order per site-year-species
  group_by(location, year, species_code) %>%
  summarise(
    max_count = max(individual_order, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Wide format: one row per site-year, one column per species
  pivot_wider(
    id_cols    = c(location, year),
    names_from = species_code,
    values_from = max_count,
    values_fill = 0
  ) %>%
  arrange(location, year)

# Save WITH SINGLETONS
write.csv(dat2_UD, "Output/Tabular Data/max_count_all_years_UDPC_with_singletons.csv", row.names = FALSE)

# Remove singletons (species occurring at <= 5 site-years)
dat3_UD <- dat2_UD[ , sapply(dat2_UD, function(x) !is.numeric(x) || sum(x != 0, na.rm = TRUE) > 3) ]

# Save WITHOUT SINGLETONS
write.csv(dat3_UD, "Output/Tabular Data/max_count_all_years_UDPC.csv", row.names = FALSE)



### CALCULATE MAX COUNT FOR TRUNCATED PC

# Load truncated tag data
dat1_LA <- read.csv("Output/Tabular Data/truncated_count_150m.csv")

# Identify species columns (everything except location & datetime)
species_cols <- dat1_LA %>%
  select(-location, -recording_date_time) %>%
  names()

# Summarize to MAX count per location per year for each species
dat1_LA_max <- dat1_LA %>%
  mutate(
    # Extract year from recording_date_time
    year = year(ymd_hms(recording_date_time))
  ) %>%
  group_by(location, year) %>%
  summarise(
    across(
      all_of(species_cols),
      ~ max(.x, na.rm = TRUE)
    ),
    .groups = "drop"   # <-- this goes here, on summarise()
  )

# Check code worked
head(dat1_LA_max)

# Save WITH SINGLETONS
write.csv(dat1_LA_max, "Output/Tabular Data/max_count_all_years_LA_with_singletons.csv", row.names = FALSE)


# Remove singletons (species occurring at <= 3 site-years)
dat2_LA_max <- dat1_LA_max[ , sapply(dat1_LA_max, function(x) !is.numeric(x) || sum(x != 0, na.rm = TRUE) > 3) ]

# Save WITHOUT SINGLETONS
write.csv(dat2_LA_max, "Output/Tabular Data/max_count_all_years_LA.csv", row.names = FALSE)












##################### 5a Single Species Analysis UD MAX COUNT


# ---
# title: "Singe species analysis - UNLIMITED DISTANCE"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: 
# ---


# Load libraries

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(tibble)
  library(viridis)
  library(DHARMa)
})


# Load data
dat1_UD <- read.csv("Output/Tabular Data/max_count_all_years_UDPC.csv")

# Extract site, treatment, plot from location
dat2_UD <- dat1_UD %>%
  separate_wider_delim(location, delim = "-", names = c("block", "treatment", "plot")) %>%
  mutate(
    # Make Year a factor for BACI analysis
    year = factor(year, levels = sort(unique(year)))
  )

# Create Site_ID and treatment_applied
pre_year <- min(as.numeric(as.character(dat2_UD$year)), na.rm = TRUE)

dat3_UD <- dat2_UD %>%
  mutate(
    # Unique experimental unit (site × assigned treatment)
    site = paste(block, treatment, sep = "-"),
    
    # BACI-ready treatment_applied: set pre-treatment year to NT
    treatment_applied = ifelse(as.numeric(as.character(year)) == pre_year, "NT", treatment),
    treatment_applied = factor(treatment_applied, levels = c("NT","LR","HR","FR"))
  )

# Reorganize dat3
spp <- c("AMRE", "ALFL", "AMRO", "CHSP", "DUFL","DEJU", "HETH", "LISP", "MGWA","OCWA", "OSFL", "PISI", "RCKI", "SWTH", "WAVI", "WETA", "WTSP", "YRWA")
spp <- names(dat3_UD)[grepl("^[A-Z]{4}$", names(dat3_UD))]
dat4_UD <- dat3_UD %>%
  select(block, site, treatment, treatment_applied, year, all_of(spp))

# Create the new 'time_period' factor 
dat4_pooled_UD <- dat4_UD %>%
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


### Create RANDOM EFFECT for geographic grouping

# Create empty df
site_group <- data.frame(
  site  = factor(dat4_pooled_UD$site),
  block = factor(dat4_pooled_UD$block)
)

# Adding grouping by sites clusters
grp1 <- c(
  "MUS124-FR","MUS124-HR","MUS124-NT",
  "MUS113-LR","MUS113-NT",
  "MUS111-HR","MUS111-NT","MUS111-LR","MUS111-FR",
  "MUS061A-FR","MUS061A-HR","MUS061A-LR"
)

grp2 <- c(
  "MUS155-NT","MUS155-LR","MUS155-FR","MUS155-HR"
)

grp3 <- c(
  "DOC266-NT","DOC266-LR","DOC266-FR","DOC266-HR",
  "DOC213-HR","DOC213-LR","DOC213-NT","DOC213-FR"
)

site_group <- site_group %>%
  dplyr::mutate(
    site_group =
      dplyr::case_when(
        site %in% grp1 ~ 1,
        site %in% grp2 ~ 2,
        site %in% grp3 ~ 3,
        TRUE ~ NA_real_
      )
  )

# Join site_group with dat4_pooled
site_group_clean <- site_group %>%
  distinct(block, site_group)

dat4_pooled_UD <- dat4_pooled_UD %>%
  left_join(site_group_clean, by = "block")

# Save
write.csv(dat4_pooled_UD, "Output/Tabular Data/max_count_all_years_UDPC_SSM.csv")


### Check var:mean ratios

var(dat4_UD$ALFL) / mean(dat4_UD$ALFL)
var(dat4_UD$AMRE) / mean(dat4_UD$AMRE)
var(dat4_UD$AMRO) / mean(dat4_UD$AMRO)
var(dat4_UD$CHSP) / mean(dat4_UD$CHSP)
var(dat4_UD$DEJU) / mean(dat4_UD$DEJU)
var(dat4_UD$DUFL) / mean(dat4_UD$DUFL)
var(dat4_UD$HETH) / mean(dat4_UD$HETH)
var(dat4_UD$LISP) / mean(dat4_UD$LISP)
var(dat4_UD$MGWA) / mean(dat4_UD$MGWA)
var(dat4_UD$OCWA) / mean(dat4_UD$OCWA)
var(dat4_UD$OSFL) / mean(dat4_UD$OSFL)
var(dat4_UD$PISI) / mean(dat4_UD$PISI)
var(dat4_UD$RCKI) / mean(dat4_UD$RCKI)
var(dat4_UD$SWTH) / mean(dat4_UD$SWTH)
var(dat4_UD$WAVI) / mean(dat4_UD$WAVI)
var(dat4_UD$WETA) / mean(dat4_UD$WETA)
var(dat4_UD$WTSP) / mean(dat4_UD$WTSP)
var(dat4_UD$YRWA) / mean(dat4_UD$YRWA)




### INDIVIDUAL SPECIES MODELS


# Initialize an empty data frame to collect all BACI contrasts
all_contrasts_df_UD <- data.frame()

# --- ROBUST HELPER FUNCTION FOR LOG CONTRAST CALCULATION ---
calculate_baci_log_contrast <- function(model, species_name) {
  # EMMs for the interaction
  eg <- emmeans(model, ~ time_period * treatment)
  
  # BACI contrast: polynomial (Before/After) × pairwise treatments
  baci_contrast <- contrast(
    eg,
    method = "pairwise",
    interaction = c("poly", "pairwise"),
    adjust = "none"
  )
  
  # --- 1) SUMMARY on LINK (LOG) SCALE ---
  sum_tbl <- summary(baci_contrast) %>%  # default type = "link"
    as_tibble()
  
  # Expect 'estimate' and 'p.value'; keep only columns we need
  if (!"estimate" %in% names(sum_tbl)) {
    stop("Expected 'estimate' column (link scale) not found in emmeans summary().")
  }
  sum_tbl <- sum_tbl %>%
    select(any_of(c("time_period_poly", "treatment_pairwise", "estimate", "p.value")))
  
  # --- 2) CONFINT on LINK (LOG) SCALE ---
  ci_tbl <- confint(baci_contrast) %>%  # default type = "link"
    as_tibble()
  
  # Harmonize CI names across emmeans versions
  if (all(c(".lower", ".upper") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = .lower, UCL = .upper)
  } else if (all(c("lower.CL", "upper.CL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = lower.CL, UCL = upper.CL)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = asymp.LCL, UCL = asymp.UCL)
  } else if (!all(c("LCL","UCL") %in% names(ci_tbl))) {
    stop("Could not find CI columns in confint() output.")
  }
  
  ci_tbl <- ci_tbl %>%
    select(any_of(c("time_period_poly", "treatment_pairwise", "LCL", "UCL")))
  
  # --- 3) JOIN + FILTER to the BACI piece you want ---
  out <- sum_tbl %>%
    inner_join(ci_tbl, by = c("time_period_poly", "treatment_pairwise")) %>%
    # 'linear' is the BACI term; some versions label it 'poly1'
    filter(time_period_poly %in% c("linear", "poly1")) %>%
    # keep comparisons vs NT regardless of ordering/format
    filter(
      grepl("NT / (LR|HR|FR)$", treatment_pairwise) |
        grepl("^(LR|HR|FR) / NT$", treatment_pairwise) |
        grepl("^NT - (LR|HR|FR)$", treatment_pairwise) |
        grepl("^(LR|HR|FR) - NT$", treatment_pairwise)
    ) %>%
    mutate(
      Species = species_name,
      # Extract the non-NT treatment wherever it appears
      Treatment = case_when(
        grepl("NT / (LR|HR|FR)$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "NT / (LR|HR|FR)$")[,2],
        grepl("^(LR|HR|FR) / NT$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^(LR|HR|FR) / NT$")[,2],
        grepl("^NT - (LR|HR|FR)$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^NT - (LR|HR|FR)$")[,2],
        grepl("^(LR|HR|FR) - NT$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^(LR|HR|FR) - NT$")[,2],
        TRUE ~ NA_character_
      ),
      LogContrast = estimate # already on log scale (link)
    ) %>%
    select(Species, Treatment, LogContrast, LCL, UCL, PValue = p.value) %>%
    arrange(Species, Treatment)
  
  return(out)
}









###### MODEL FORMULA: Response ~ treatment * time_period + (1|site) + (1|block)

# ==============================================================================
# 1. AMRE (American Redstart) | VMR: 0.853 | Family: Com-Poisson
# ==============================================================================
AMRE_model_UD <- glmmTMB(
  AMRE ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link="log"),
  data = dat4_pooled_UD
)

AMRE_contrasts_UD <- calculate_baci_log_contrast(AMRE_model_UD, "AMRE")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, AMRE_contrasts_UD)
summary(AMRE_model_UD)

### Assess model residuals

set.seed(123)  # for reproducibility of the simulations
res_AMRE_UD <- simulateResiduals(
  fittedModel = AMRE_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_AMRE_UD)

# ==============================================================================
# 1. ALFL  | VMR: 0.853 | Family: Poisson
# ==============================================================================
ALFL_model_UD <- glmmTMB(
  ALFL ~ treatment * time_period + (1|site) + (1|block),
  family = compois(), 
  data = dat4_pooled_UD
)
ALFL_contrasts_UD <- calculate_baci_log_contrast(ALFL_model_UD, "ALFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, ALFL_contrasts_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_ALFL_UD <- simulateResiduals(
  fittedModel = ALFL_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_ALFL_UD)

# ==============================================================================
# 2. AMRO (American Robin) | VMR: 0.439 | Family: Com-Poisson
# ==============================================================================
AMRO_model_UD <- glmmTMB(
  AMRO ~ treatment * time_period + (1|site) + (1|block),
  family = compois(), 
  data = dat4_pooled_UD
)
AMRO_contrasts_UD <- calculate_baci_log_contrast(AMRO_model_UD, "AMRO")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, AMRO_contrasts_UD)
summary(AMRO_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_AMRO_UD <- simulateResiduals(
  fittedModel = AMRO_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_AMRO_UD)


# ==============================================================================
# 4. DEJU (Dark-eyed Junco) | VMR: 0.555 | Family: Com-Poisson
# ==============================================================================
DEJU_model_UD <- glmmTMB(
  DEJU ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
DEJU_contrasts_UD <- calculate_baci_log_contrast(DEJU_model_UD, "DEJU")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, DEJU_contrasts_UD)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DEJU_UD <- simulateResiduals(
  fittedModel = DEJU_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DEJU_UD)


# ==============================================================================
# 5. DUFL (Dusky Flycatcher) | VMR: 0.198 | Family: Zero-inflated compois
# ==============================================================================
DUFL_model_UD <- glmmTMB(
  DUFL ~ treatment * time_period + 
    (1 | site) + (1 | block),
  ziformula = ~ 1,
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

DUFL_contrasts_UD <- calculate_baci_log_contrast(DUFL_model_UD, "DUFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, DUFL_contrasts_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DUFL_UD <- simulateResiduals(
  fittedModel = DUFL_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DUFL_UD)

# ==============================================================================
# 7. LISP (Lincoln's Sparrow) | VMR: 0.568 | Family: Com-Poisson
# ==============================================================================
LISP_model_UD <- glmmTMB(
  LISP ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
LISP_contrasts_UD <- calculate_baci_log_contrast(LISP_model_UD, "LISP")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, LISP_contrasts_UD)
summary(LISP_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_LISP_UD <- simulateResiduals(
  fittedModel = LISP_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_LISP_UD)


# ==============================================================================
# 9. OCWA (Orange-crowned Warbler) | VMR: 0.337 | Family: Com-Poisson
# ==============================================================================
OCWA_model_UD <- glmmTMB(
  OCWA ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
OCWA_contrasts_UD <- calculate_baci_log_contrast(OCWA_model_UD, "OCWA")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, OCWA_contrasts_UD)
summary(OCWA_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OCWA_UD <- simulateResiduals(
  fittedModel = OCWA_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OCWA_UD)



# ==============================================================================
# 10. OSFL (OSFL) | VMR:  | Family: Com-Poisson
# ==============================================================================
OSFL_model_UD <- glmmTMB(
  OSFL ~ treatment * time_period + (1|site) + (1|block),
  family = poisson, 
  data = dat4_pooled_UD
)
OSFL_contrasts_UD <- calculate_baci_log_contrast(OSFL_model_UD, "OSFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, OSFL_contrasts_UD)
summary(OSFL_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OSFL_UD <- simulateResiduals(
  fittedModel = OSFL_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OSFL_UD)


# ==============================================================================
# 11. RCKI (Ruby-crowned Kinglet) | VMR: 0.638 | Family: Poisson
# ==============================================================================
RCKI_model_UD <- glmmTMB(
  RCKI ~ treatment * time_period + (1|site) + (1|block),
  family = compois(), # VMR >= 1.05
  data = dat4_pooled_UD
)
RCKI_contrasts_UD <- calculate_baci_log_contrast(RCKI_model_UD, "RCKI")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, RCKI_contrasts_UD)
summary(RCKI_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_RCKI_UD <- simulateResiduals(
  fittedModel = RCKI_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_RCKI_UD)

# ==============================================================================
# 12. SWTH (Swainson's Thrush) | VMR: 0.254 | Family: Com-Poisson
# ==============================================================================
SWTH_model_UD <- glmmTMB(
  SWTH ~ treatment * time_period + 
    (1 | site) + (1 | block),
  ziformula = ~ 1,
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

SWTH_contrasts_UD <- calculate_baci_log_contrast(SWTH_model_UD, "SWTH")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, SWTH_contrasts_UD)
summary(SWTH_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_SWTH_UD <- simulateResiduals(
  fittedModel = SWTH_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_SWTH_UD)

# ==============================================================================
# 13. WAVI | VMR: 0.725 | Family: Com-Poisson
# ==============================================================================
WAVI_model_UD <- glmmTMB(
  WETA ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
WAVI_contrasts_UD <- calculate_baci_log_contrast(WAVI_model_UD, "WAVI")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, WAVI_contrasts_UD)
summary(WAVI_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WAVI_UD <- simulateResiduals(
  fittedModel = WAVI_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WAVI_UD)

# ==============================================================================
# 14. WETA (Western Tanager) | VMR: 0.725 | Family: Com-Poisson
# ==============================================================================
WETA_model_UD <- glmmTMB(
  WETA ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
WETA_contrasts_UD <- calculate_baci_log_contrast(WETA_model_UD, "WETA")
all_contrasts_df <- bind_rows(all_contrasts_df_UD, WETA_contrasts_UD)
summary(WETA_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WETA_UD <- simulateResiduals(
  fittedModel = WETA_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WETA_UD)


# ==============================================================================
# 15. WTSP (White-throated Sparrow) | VMR: 0.324 | Family: Com-Poisson
# ==============================================================================
WTSP_model_UD <- glmmTMB(
  WTSP ~ treatment * time_period + 
    (1 | site) + (1 | block),
  ziformula = ~ 1,
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

WTSP_contrasts_UD <- calculate_baci_log_contrast(WTSP_model_UD, "WTSP")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, WTSP_contrasts_UD)
summary(WTSP_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WTSP_UD <- simulateResiduals(
  fittedModel = WTSP_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WTSP_UD)


# ==============================================================================
# 16. YRWA (Yellow-rumped Warbler) | VMR: 0.574 | Family: Com-Poisson
# ==============================================================================
YRWA_model_UD <- glmmTMB(
  YRWA ~ treatment * time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
YRWA_contrasts_UD <- calculate_baci_log_contrast(YRWA_model_UD, "YRWA")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, YRWA_contrasts_UD)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_YRWA_UD <- simulateResiduals(
  fittedModel = YRWA_model_UD,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_YRWA_UD)


print("--- BACI Contrast Ratios (Response Ratio) Summary (Exponentiated) ---")


# Arrange results for a cleaner printout - NOW USING LOGCONTRAST
print(all_contrasts_df_UD %>% 
        mutate(across(c(LogContrast, LCL, UCL), ~round(., 3))) %>%
        arrange(Species, Treatment))




########## PLOT RESPONSES BY TREATMENT
# Remove OSFL
all_contrasts_df_UD <- all_contrasts_df_UD %>%
  filter(!Species %in% c("OSFL"))

contrast_plot_log_treat_UD <- all_contrasts_df_UD %>%
  filter(!is.na(LogContrast)) %>%
  ggplot(aes(x = LogContrast, y = Species, color = Treatment)) +
  
  # Reference line at zero (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  
  # Whiskers and points
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  
  # Facet by treatment (3 panels)
  facet_wrap(~ Treatment, ncol = 3, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +
  
  labs(
    title = "BACI Contrasts by Treatment (Unlimited Distance)",
    x = "Contrasts (Log Scale)",
    y = "Species"
  ) +
  
  # Consistent treatment colors
  scale_color_viridis_d(option = "D", end = 0.8, direction = 1) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # color corresponds directly to facet
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

print("--- Plotting BACI Contrasts by Treatment (Log Scale) ---")
print(contrast_plot_log_treat_UD)

# Save
save_plot("Figures/Species-level models/BACI_SSM_UDPC.png", contrast_plot_log_treat_UD, base_width = 7, base_height = 5)




########## PLOT RESPONSES BY TREATMENT

# Remove spp w/ problematic CIs
all_contrasts_df_UD_filtered <- all_contrasts_df_UD %>%
  filter(!Species %in% c("OSFL")) %>%
  filter(!is.na(LogContrast)) %>% # Filter NA contrasts early
  # --- FIX: Explicitly set the factor levels for the desired order ---
  mutate(
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR"))
  )

# Contrasts
contrast_plot_log_treat_UD <- all_contrasts_df_UD_filtered %>%
  ggplot(aes(x = LogContrast, y = Species, color = Treatment)) +
  
  # Reference line at zero (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  
  # Whiskers and points
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  
  # Facet by treatment (3 panels). The order is now controlled by the factor levels set above.
  facet_wrap(~ Treatment, ncol = 3, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +
  
  labs(
    title = "BACI Contrasts by Treatment (Unlimited Distance)",
    x = "Contrasts (Log Scale)",
    y = "Species"
  ) +
  
  # Consistent treatment colors
  scale_color_viridis_d(option = "D", end = 0.8, direction = 1) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none", # color corresponds directly to facet
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

print("--- Plotting BACI Contrasts by Treatment (Log Scale) ---")
print(contrast_plot_log_treat_UD)

# Save
save_plot("Figures/Species-level models/BACI_SSM_UDPC.png", contrast_plot_log_treat_UD, base_width = 7, base_height = 5)




























########### 5B SINGLE SPECIES MODELS - LA - MAX COUNT

# ---
# title: "Single-species analysis limited amplitude"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script create extrapolation and rarefaction curves for Hill numbers of the order q = 0, 1, and 2
# ---


# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(tibble)
  library(viridis)
  library(DHARMa)
})

# Load data
dat1_LA <- read.csv("Output/Tabular Data/max_count_all_years_LA.csv")

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
spp <- c("AMRE", "ALFL", "AMRO", "CHSP", "DUFL","DEJU", "HETH", "LISP", "MGWA","OCWA", "OSFL", "PISI", "RCKI", "SWTH", "WAVI", "WETA", "WTSP", "YRWA")
spp <- names(dat3_LA)[grepl("^[A-Z]{4}$", names(dat3_LA))]
dat4_LA <- dat3_LA %>%
  select(block, site, treatment, treatment_applied, year, all_of(spp))

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


### Create RANDOM EFFECT for geographic grouping

# Create empty df
site_group <- data.frame(
  site  = factor(dat4_pooled_LA$site),
  block = factor(dat4_pooled_LA$block)
)

# Adding grouping by sites clusters
grp1 <- c(
  "MUS124-FR","MUS124-HR","MUS124-NT",
  "MUS113-LR","MUS113-NT",
  "MUS111-HR","MUS111-NT","MUS111-LR","MUS111-FR",
  "MUS061A-FR","MUS061A-HR","MUS061A-LR"
)

grp2 <- c(
  "MUS155-NT","MUS155-LR","MUS155-FR","MUS155-HR"
)

grp3 <- c(
  "DOC266-NT","DOC266-LR","DOC266-FR","DOC266-HR",
  "DOC213-HR","DOC213-LR","DOC213-NT","DOC213-FR"
)

site_group <- site_group %>%
  dplyr::mutate(
    site_group =
      dplyr::case_when(
        site %in% grp1 ~ 1,
        site %in% grp2 ~ 2,
        site %in% grp3 ~ 3,
        TRUE ~ NA_real_
      )
  )

# Join site_group with dat4_pooled
site_group_clean <- site_group %>%
  distinct(block, site_group)

dat4_pooled_LA <- dat4_pooled_LA %>%
  left_join(site_group_clean, by = "block")

# Save
write.csv(dat4_pooled_LA, "Output/Tabular Data/max_count_all_years_LAPC_SSM.csv")


### Check var:mean ratios

var(dat4_LA$ALFL) / mean(dat4_LA$ALFL)
var(dat4_LA$AMRE) / mean(dat4_LA$AMRE)
var(dat4_LA$AMRO) / mean(dat4_LA$AMRO)
var(dat4_LA$CHSP) / mean(dat4_LA$CHSP)
var(dat4_LA$DEJU) / mean(dat4_LA$DEJU)
var(dat4_LA$DUFL) / mean(dat4_LA$DUFL)
var(dat4_LA$HETH) / mean(dat4_LA$HETH)
var(dat4_LA$LISP) / mean(dat4_LA$LISP)
var(dat4_LA$OCWA) / mean(dat4_LA$OCWA)
var(dat4_LA$OSFL) / mean(dat4_LA$OSFL)
var(dat4_LA$RCKI) / mean(dat4_LA$RCKI)
var(dat4_LA$SWTH) / mean(dat4_LA$SWTH)
var(dat4_LA$WAVI) / mean(dat4_LA$WAVI)
var(dat4_LA$WETA) / mean(dat4_LA$WETA)
var(dat4_LA$WTSP) / mean(dat4_LA$WTSP)
var(dat4_LA$YRWA) / mean(dat4_LA$YRWA)



### INDIVIDUAL SPECIES MODELS


# Initialize an empty data frame to collect all BACI contrasts
all_contrasts_df_LA <- data.frame()

# --- ROBUST HELPER FUNCTION FOR LOG CONTRAST CALCULATION ---
calculate_baci_log_contrast <- function(model, species_name) {
  # EMMs for the interaction
  eg <- emmeans(model, ~ time_period * treatment)
  
  # BACI contrast: polynomial (Before/After) × pairwise treatments
  baci_contrast <- contrast(
    eg,
    method = "pairwise",
    interaction = c("poly", "pairwise"),
    adjust = "none"
  )
  
  # --- 1) SUMMARY on LINK (LOG) SCALE ---
  sum_tbl <- summary(baci_contrast) %>%  # default type = "link"
    as_tibble()
  
  # Expect 'estimate' and 'p.value'; keep only columns we need
  if (!"estimate" %in% names(sum_tbl)) {
    stop("Expected 'estimate' column (link scale) not found in emmeans summary().")
  }
  sum_tbl <- sum_tbl %>%
    select(any_of(c("time_period_poly", "treatment_pairwise", "estimate", "p.value")))
  
  # --- 2) CONFINT on LINK (LOG) SCALE ---
  ci_tbl <- confint(baci_contrast) %>%  # default type = "link"
    as_tibble()
  
  # Harmonize CI names across emmeans versions
  if (all(c(".lower", ".upper") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = .lower, UCL = .upper)
  } else if (all(c("lower.CL", "upper.CL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = lower.CL, UCL = upper.CL)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = asymp.LCL, UCL = asymp.UCL)
  } else if (!all(c("LCL","UCL") %in% names(ci_tbl))) {
    stop("Could not find CI columns in confint() output.")
  }
  
  ci_tbl <- ci_tbl %>%
    select(any_of(c("time_period_poly", "treatment_pairwise", "LCL", "UCL")))
  
  # --- 3) JOIN + FILTER to the BACI piece you want ---
  out <- sum_tbl %>%
    inner_join(ci_tbl, by = c("time_period_poly", "treatment_pairwise")) %>%
    # 'linear' is the BACI term; some versions label it 'poly1'
    filter(time_period_poly %in% c("linear", "poly1")) %>%
    # keep comparisons vs NT regardless of ordering/format
    filter(
      grepl("NT / (LR|HR|FR)$", treatment_pairwise) |
        grepl("^(LR|HR|FR) / NT$", treatment_pairwise) |
        grepl("^NT - (LR|HR|FR)$", treatment_pairwise) |
        grepl("^(LR|HR|FR) - NT$", treatment_pairwise)
    ) %>%
    mutate(
      Species = species_name,
      # Extract the non-NT treatment wherever it appears
      Treatment = case_when(
        grepl("NT / (LR|HR|FR)$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "NT / (LR|HR|FR)$")[,2],
        grepl("^(LR|HR|FR) / NT$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^(LR|HR|FR) / NT$")[,2],
        grepl("^NT - (LR|HR|FR)$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^NT - (LR|HR|FR)$")[,2],
        grepl("^(LR|HR|FR) - NT$", treatment_pairwise) ~
          stringr::str_match(treatment_pairwise, "^(LR|HR|FR) - NT$")[,2],
        TRUE ~ NA_character_
      ),
      LogContrast = estimate # already on log scale (link)
    ) %>%
    select(Species, Treatment, LogContrast, LCL, UCL, PValue = p.value) %>%
    arrange(Species, Treatment)
  
  return(out)
}









###### MODEL FORMULA: Response ~ treatment * time_period + (1|site) + (1|block)

# ==============================================================================
# 1. AMRE (American Redstart) | VMR: 0.879 | Family: Com-Poisson
# ==============================================================================
AMRE_model_LA <- glmmTMB(
  AMRE ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = compois(link="log"),
  data = dat4_pooled_LA
)

AMRE_contrasts_LA <- calculate_baci_log_contrast(AMRE_model_LA, "AMRE")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, AMRE_contrasts_LA)
summary(AMRE_model_LA)

### Assess model residuals

set.seed(123)  # for reproducibility of the simulations
res_AMRE_LA <- simulateResiduals(
  fittedModel = AMRE_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_AMRE_LA)

# ==============================================================================
# 1. ALFL  | VMR: 0.889 | Family: Poisson
# ==============================================================================
ALFL_model_LA <- glmmTMB(
  ALFL ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(),
  ziformula = ~ 1,
  data = dat4_pooled_LA
)
ALFL_contrasts_LA <- calculate_baci_log_contrast(ALFL_model_LA, "ALFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, ALFL_contrasts_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_ALFL_LA <- simulateResiduals(
  fittedModel = ALFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_ALFL_LA)

# ==============================================================================
# 2. AMRO (American Robin) | VMR: 0.831 | Family: Com-Poisson
# ==============================================================================
AMRO_model_LA <- glmmTMB(
  AMRO ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(),
  data = dat4_pooled_LA
)
AMRO_contrasts_LA <- calculate_baci_log_contrast(AMRO_model_LA, "AMRO")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, AMRO_contrasts_LA)
summary(AMRO_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_AMRO_LA <- simulateResiduals(
  fittedModel = AMRO_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_AMRO_LA)


# ==============================================================================
# 4. DEJU (Dark-eyed Junco) | VMR: 0.661 | Family: Com-Poisson
# ==============================================================================
DEJU_model_LA <- glmmTMB(
  DEJU ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_LA
)
DEJU_contrasts_LA <- calculate_baci_log_contrast(DEJU_model_LA, "DEJU")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, DEJU_contrasts_LA)
summary(DEJU_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DEJU_LA <- simulateResiduals(
  fittedModel = DEJU_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DEJU_LA)


# ==============================================================================
# 5. DUFL (Dusky Flycatcher) | VMR: 0.374 | Family: Zero-inflated compois
# ==============================================================================
DUFL_model_LA <- glmmTMB(
  DUFL ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(),
  ziformula = ~ 1,
  data   = dat4_pooled_LA
)

DUFL_contrasts_LA <- calculate_baci_log_contrast(DUFL_model_LA, "DUFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, DUFL_contrasts_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DUFL_LA <- simulateResiduals(
  fittedModel = DUFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DUFL_LA)


# ==============================================================================
# 7. LISP (Lincoln's Sparrow) | VMR: 0.732 | Family: Com-Poisson
# ==============================================================================
LISP_model_LA <- glmmTMB(
  LISP ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(), 
  data = dat4_pooled_LA
)

LISP_contrasts_LA <- calculate_baci_log_contrast(LISP_model_LA, "LISP")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, LISP_contrasts_LA)
summary(LISP_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_LISP_LA <- simulateResiduals(
  fittedModel = LISP_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_LISP_LA)


# ==============================================================================
# 9. OCWA (Orange-crowned Warbler) | VMR: 0.574 | Family: Com-Poisson
# ==============================================================================
OCWA_model_LA <- glmmTMB(
  OCWA ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_LA
)
OCWA_contrasts_LA <- calculate_baci_log_contrast(OCWA_model_LA, "OCWA")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, OCWA_contrasts_LA)
summary(OCWA_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OCWA_LA <- simulateResiduals(
  fittedModel = OCWA_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OCWA_LA)



# ==============================================================================
# 10. OSFL (OSFL) | VMR: 0.859 | Family: Com-Poisson
# ==============================================================================
OSFL_model_LA <- glmmTMB(
  OSFL ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(), 
  data = dat4_pooled_LA
)
OSFL_contrasts_LA <- calculate_baci_log_contrast(OSFL_model_LA, "OSFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, OSFL_contrasts_LA)
summary(OSFL_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OSFL_LA <- simulateResiduals(
  fittedModel = OSFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OSFL_LA)


# ==============================================================================
# 11. RCKI (Ruby-crowned Kinglet) | VMR: 0.930 | Family: Poisson
# ==============================================================================
RCKI_model_LA <- glmmTMB(
  RCKI ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(), 
  data = dat4_pooled_LA
)
RCKI_contrasts_LA <- calculate_baci_log_contrast(RCKI_model_LA, "RCKI")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, RCKI_contrasts_LA)
summary(RCKI_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_RCKI_LA <- simulateResiduals(
  fittedModel = RCKI_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_RCKI_LA)

# ==============================================================================
# 12. SWTH (Swainson's Thrush) | VMR: 0.878 | Family: Com-Poisson
# ==============================================================================
SWTH_model_LA <- glmmTMB(
  SWTH ~ treatment  + time_period + treatment*time_period + 
    (1 | site) + (1 | block),
  family = poisson(),
  data = dat4_pooled_LA
)

SWTH_contrasts_LA <- calculate_baci_log_contrast(SWTH_model_LA, "SWTH")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, SWTH_contrasts_LA)
summary(SWTH_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_SWTH_LA <- simulateResiduals(
  fittedModel = SWTH_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_SWTH_LA)

# ==============================================================================
# 13. WAVI | VMR: 0.939 | Family: Poisson
# ==============================================================================
WAVI_model_LA <- glmmTMB(
  WETA ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(), 
  data = dat4_pooled_LA
)
WAVI_contrasts_LA <- calculate_baci_log_contrast(WAVI_model_LA, "WAVI")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, WAVI_contrasts_LA)
summary(WAVI_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WAVI_LA <- simulateResiduals(
  fittedModel = WAVI_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WAVI_LA)

# ==============================================================================
# 14. WETA (Western Tanager) | VMR: 1.113 | Family: Com-Poisson
# ==============================================================================
WETA_model_LA <- glmmTMB(
  WETA ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = poisson(), 
  data = dat4_pooled_LA
)
WETA_contrasts_LA <- calculate_baci_log_contrast(WETA_model_LA, "WETA")
all_contrasts_df <- bind_rows(all_contrasts_df_LA, WETA_contrasts_LA)
summary(WETA_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WETA_LA <- simulateResiduals(
  fittedModel = WETA_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WETA_LA)


# ==============================================================================
# 15. WTSP (White-throated Sparrow) | VMR: 0.866 | Family: Com-Poisson
# ==============================================================================
WTSP_model_LA <- glmmTMB(
  WTSP ~ treatment  + time_period + treatment*time_period + 
    (1 | site) + (1 | block),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)

WTSP_contrasts_LA <- calculate_baci_log_contrast(WTSP_model_LA, "WTSP")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, WTSP_contrasts_LA)
summary(WTSP_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WTSP_LA <- simulateResiduals(
  fittedModel = WTSP_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WTSP_LA)


# ==============================================================================
# 16. YRWA (Yellow-rumped Warbler) | VMR: 0.631 | Family: Com-Poisson
# ==============================================================================
YRWA_model_LA <- glmmTMB(
  YRWA ~ treatment  + time_period + treatment*time_period + (1|site) + (1|block),
  family = compois(link = "log"), 
  data = dat4_pooled_LA
)
YRWA_contrasts_LA <- calculate_baci_log_contrast(YRWA_model_LA, "YRWA")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, YRWA_contrasts_LA)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_YRWA_LA <- simulateResiduals(
  fittedModel = YRWA_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_YRWA_LA)


print("--- BACI Contrast Ratios (Response Ratio) Summary (Exponentiated) ---")


# Arrange results for a cleaner printout - NOW USING LOGCONTRAST
print(all_contrasts_df_LA %>% 
        mutate(across(c(LogContrast, LCL, UCL), ~round(., 3))) %>%
        arrange(Species, Treatment))









########## PLOT RESPONSES BY TREATMENT

# Remove spp w/ problematic CIs
all_contrasts_df_LA_filtered <- all_contrasts_df_LA %>%
  filter(!Species %in% c("WAVI", "OSFL", "LISP", "AMRO", "RCKI")) %>%
  filter(!is.na(LogContrast)) %>% # Filter NA contrasts early
  # --- FIX: Explicitly set the factor levels for the desired order ---
  mutate(
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR"))
  )

# Contrasts
contrast_plot_log_treat_LA <- all_contrasts_df_LA_filtered %>%
  ggplot(aes(x = LogContrast, y = Species, color = Treatment)) +
  
  # Reference line at zero (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  
  # Whiskers and points
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  
  # Facet by treatment (3 panels). The order is now controlled by the factor levels set above.
  facet_wrap(~ Treatment, ncol = 3, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +
  
  labs(
    title = "BACI Contrasts by Treatment (Limited Amplitude)",
    x = "Contrasts (Log Scale)",
    y = "Species"
  ) +
  
  # Consistent treatment colors
  scale_color_viridis_d(option = "D", end = 0.8, direction = 1) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none", # color corresponds directly to facet
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

print("--- Plotting BACI Contrasts by Treatment (Log Scale) ---")
print(contrast_plot_log_treat_LA)

# Save
save_plot("Figures/Species-level models/BACI_SSM_LAPC.png", contrast_plot_log_treat_LA, base_width = 7, base_height = 5)
















############ 6a GAMMA DIV UD







###### GOOD START ##########


################################### MAX COUNT POOLING2024/2025 w/ BOOTSTRAPPING

rm(list = ls())  # Removes all objects from the environment

# Install and load packages
#install.packages("iNEXT")
library(tidyverse)
library(iNEXT)
library(cowplot)
library(forcats)
library(patchwork)

# Load data
abundance_all_years_UDPC <- read.csv("Output/Tabular Data/max_count_all_years_UDPC.csv")

# Check structure
str(abundance_all_years_UDPC)

# ---- Identify species columns ----
species_cols <- names(abundance_all_years_UDPC)[grepl("^[A-Z]{4}$", names(abundance_all_years_UDPC))]

# Coerce species columns to numeric and replace NAs with 0
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))

# Create treatment column
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location))

# Rearrange columns
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  select(location, treatment, year, all_of(species_cols))

# Split into a dataframe each for 2023, 2024, and 2025
abund_2023_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2023")
abund_2024_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2024")
abund_2025_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2025")

## 1. Sum species per treatment per year -------------------------------

sum_2023_UD <- abund_2023_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

sum_2024_UD <- abund_2024_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

sum_2025_UD <- abund_2025_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

## 2. Pooled AFTER = 2024 + 2025 per treatment ------------------------

sum_after_UD <- bind_rows(sum_2024_UD, sum_2025_UD) %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

## 3. Convert to iNEXT abundance lists --------------------------------

# BEFORE: 2023
iNext_before_UD <- sum_2023_UD %>%
  tidyr::nest(.by = treatment) %>%
  mutate(
    vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )
  ) %>%
  transmute(treatment, vec) %>%
  deframe()

# AFTER: pooled 2024+2025
iNext_after_UD <- sum_after_UD %>%
  tidyr::nest(.by = treatment) %>%
  mutate(
    vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )
  ) %>%
  transmute(treatment, vec) %>%
  deframe()

# Drop zero-abundance species 
iNext_before_UD <- lapply(iNext_before_UD, \(x) x[x > 0])
iNext_after_UD  <- lapply(iNext_after_UD,  \(x) x[x > 0])

## Optional iNEXT curves
inext_before_UD <- iNEXT(iNext_before_UD, datatype = "abundance", q = c(0, 1, 2))
inext_after_UD  <- iNEXT(iNext_after_UD,  datatype = "abundance", q = c(0, 1, 2))

# Plot (rarefaction/extrapolation)
RE_abund_before_UD <- ggiNEXT(inext_before_UD, type = 1, facet.var = "Order.q")
RE_abund_after_UD  <- ggiNEXT(inext_after_UD,  type = 1, facet.var = "Order.q")

# Save
save_plot("Figures/RE plots/RE_plot_before_UD.png",
          RE_abund_before_UD, base_width = 7, base_height = 5)
save_plot("Figures/RE plots/RE_plot_after_UD.png",
          RE_abund_after_UD,  base_width = 7, base_height = 5)

## 4. Find common sample size (min n across periods) ------------------

info_before <- DataInfo(iNext_before_UD, datatype = "abundance")
info_after  <- DataInfo(iNext_after_UD,  datatype = "abundance")

# n = total individuals (sum of abundances) in each assemblage
min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
min_n  # inspect

## 5. Standardize Hill numbers at this common sample size --------------

equal_before_UD <- estimateD(
  iNext_before_UD,
  datatype = "abundance",
  base     = "size",
  level    = min_n,
  conf     = 0.95
) %>%
  mutate(Period = "Before")

equal_after_UD <- estimateD(
  iNext_after_UD,
  datatype = "abundance",
  base     = "size",
  level    = min_n,
  conf     = 0.95
) %>%
  mutate(Period = "After")

equal_size_UD <- bind_rows(equal_before_UD, equal_after_UD)

## 6. BACI contrasts: (After-Before)_Treat - (After-Before)_NT ---------

df0_UD <- equal_size_UD %>%
  select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, Period) %>%
  mutate(
    Assemblage = factor(Assemblage, levels = c("NT", "LR", "HR", "FR")),
    Period     = factor(Period, levels = c("Before", "After")),
    Order.q    = as.integer(Order.q),
    SE         = (qD.UCL - qD.LCL) / (2 * 1.96),
    Var        = SE^2
  )

treatments <- c("LR", "HR", "FR")
control    <- "NT"

wide_UD <- df0_UD %>%
  select(Assemblage, Order.q, Period, qD, Var) %>%
  distinct() %>%  # guard against accidental duplicates
  tidyr::pivot_wider(
    names_from  = Period,
    values_from = c(qD, Var),
    names_sep   = "_"
  ) %>%
  mutate(
    d_Treat = qD_After - qD_Before,
    v_Treat = Var_After + Var_Before
  )

ctrl_UD <- wide_UD %>%
  filter(Assemblage == control) %>%
  select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)

baci_results_UD <- wide_UD %>%
  filter(Assemblage %in% treatments) %>%
  left_join(ctrl_UD, by = "Order.q", relationship = "many-to-one") %>%
  transmute(
    Order.q  = factor(Order.q, levels = c(0, 1, 2),
                      labels = c("q = 0", "q = 1", "q = 2")),
    Treatment = factor(Assemblage, levels = c("FR", "LR", "HR")),
    Estimate  = d_Treat - d_Ctrl,
    SE        = sqrt(v_Treat + v_Ctrl),
    LCL       = Estimate - 1.96 * SE,
    UCL       = Estimate + 1.96 * SE
  )


## 7. Plot -------------------------------------------------------------

p_baci_UD <- ggplot(baci_results_UD,
                    aes(x = Estimate, y = Treatment, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 height = 0.15, color = "black") +
  geom_point(size = 3) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Limited Amplitude (2024/2025 pooled)",
    x = "Change in Hill Numbers",
    y = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "right",
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold")
  )

print(p_baci_UD)

save_plot("Figures/RE plots/BACI_contrasts_UD_pooled_size.png",
          p_baci_UD, base_width = 7, base_height = 5)






############ GEMINI SIMPLIFICATION


## =====================================================================
## BOOTSTRAP BACI CONTRASTS (UNLIMITED DISTANCE, UDPC)
## =====================================================================

library(tidyverse)
library(iNEXT)
library(ggplot2) # Include plotting package

# ---------------------------------------------------------------------
# 1. Helper: compute BACI contrasts from a *given* raw data frame
#    (The definition of this function remains UNCHANGED from your working code)
# ---------------------------------------------------------------------

compute_baci_udpc <- function(dat) {
  # Identify species columns
  species_cols <- names(dat)[grepl("^[A-Z]{4}$", names(dat))]
  
  # Make sure species are numeric and NAs are zero
  dat <- dat %>%
    mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))
  
  # Split by year
  abund_2023 <- dat %>% filter(year == "2023")
  abund_2024 <- dat %>% filter(year == "2024")
  abund_2025 <- dat %>% filter(year == "2025")
  
  # If any period is empty, bail out with empty tibble
  if (nrow(abund_2023) == 0 || (nrow(abund_2024) == 0 & nrow(abund_2025) == 0)) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Sum species per treatment per year ----
  sum_2023 <- abund_2023 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  sum_2024 <- abund_2024 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  sum_2025 <- abund_2025 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  ## ---- Pooled AFTER = 2024 + 2025 per treatment ----
  sum_after <- bind_rows(sum_2024, sum_2025) %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  # If some treatment disappears completely, bail
  if (nrow(sum_2023) == 0 || nrow(sum_after) == 0) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Build iNEXT abundance lists (Before/After) ----
  iNext_before <- sum_2023 %>%
    tidyr::nest(.by = treatment) %>%
    mutate(vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )) %>%
    transmute(treatment, vec) %>%
    deframe()
  
  iNext_after <- sum_after %>%
    tidyr::nest(.by = treatment) %>%
    mutate(vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )) %>%
    transmute(treatment, vec) %>%
    deframe()
  
  # Drop zero-abundance species
  iNext_before <- lapply(iNext_before, \(x) x[x > 0])
  iNext_after  <- lapply(iNext_after,  \(x) x[x > 0])
  
  ## ---- Common sample size (min n across periods) ----
  info_before <- DataInfo(iNext_before, datatype = "abundance")
  info_after  <- DataInfo(iNext_after,  datatype = "abundance")
  
  min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
  
  # If min_n is invalid, bail
  if (!is.finite(min_n) || min_n <= 0) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Standardize Hill numbers at common sample size ----
  equal_before <- estimateD(
    iNext_before,
    datatype = "abundance",
    base     = "size",
    level    = min_n,
    conf     = 0.95
  ) %>%
    mutate(Period = "Before")
  
  equal_after <- estimateD(
    iNext_after,
    datatype = "abundance",
    base     = "size",
    level    = min_n,
    conf     = 0.95
  ) %>%
    mutate(Period = "After")
  
  equal_size <- bind_rows(equal_before, equal_after)
  
  ## ---- BACI contrasts: (After-Before)_Treat - (After-Before)_NT ----
  
  df0 <- equal_size %>%
    select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, Period) %>%
    mutate(
      Assemblage = factor(Assemblage, levels = c("NT", "LR", "HR", "FR")),
      Period     = factor(Period, levels = c("Before", "After")),
      Order.q    = as.integer(Order.q),
      SE         = (qD.UCL - qD.LCL) / (2 * 1.96),
      Var        = SE^2
    )
  
  treatments <- c("LR", "HR", "FR")
  control    <- "NT"
  
  wide <- df0 %>%
    select(Assemblage, Order.q, Period, qD, Var) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from  = Period,
      values_from = c(qD, Var),
      names_sep   = "_"
    ) %>%
    mutate(
      d_Treat = qD_After - qD_Before,
      v_Treat = Var_After + Var_Before
    )
  
  ctrl <- wide %>%
    filter(Assemblage == control) %>%
    select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)
  
  baci_results <- wide %>%
    filter(Assemblage %in% treatments) %>%
    left_join(ctrl, by = "Order.q", relationship = "many-to-one") %>%
    transmute(
      Order.q   = factor(Order.q, levels = c(0, 1, 2),
                         labels = c("q = 0", "q = 1", "q = 2")),
      Treatment = factor(Assemblage, levels = c("FR", "LR", "HR")),
      Estimate  = d_Treat - d_Ctrl,
      SE        = sqrt(v_Treat + v_Ctrl),
      LCL       = Estimate - 1.96 * SE,
      UCL       = Estimate + 1.96 * SE
    )
  
  return(baci_results)
}

## =====================================================================
## BOOTSTRAP BACI CONTRASTS (UNLIMITED DISTANCE, UDPC) – BOOTSTRAP PLOT
## =====================================================================

## 1. Bootstrap: resample locations within treatment, keep all years ----

set.seed(123)    # for reproducibility
B <- 1000        # number of bootstrap replicates

# Sampling units: locations within each treatment
loc_tbl <- abundance_all_years_UDPC %>%
  distinct(location, treatment)

boot_list <- vector("list", B)

# Define the structure for a failed run (9 rows: q=0,1,2 x T=FR,LR,HR)
failure_tibble_structure <- tibble(
  Order.q = factor(rep(c("q = 0", "q = 1", "q = 2"), each = 3), levels = c("q = 0", "q = 1", "q = 2")),
  Treatment = factor(rep(c("FR", "LR", "HR"), times = 3), levels = c("FR", "LR", "HR")),
  Estimate = NA_real_,
  SE = NA_real_,
  LCL = NA_real_,
  UCL = NA_real_
)


for (b in seq_len(B)) {
  # Resample locations *within* each treatment, with replacement
  sampled_locs <- loc_tbl %>%
    group_by(treatment) %>%
    summarise(
      location = sample(location, size = n(), replace = TRUE),
      .groups  = "drop"
    )
  
  # Join back to full data: keeps all years for each selected location
  boot_dat <- sampled_locs %>%
    left_join(abundance_all_years_UDPC,
              by = c("location", "treatment"))
  
  # Compute BACI for this bootstrap sample, with error handling
  boot_baci <- tryCatch(
    {
      # SUCCESS: Calculate BACI and add boot ID
      compute_baci_udpc(boot_dat) %>%
        mutate(boot = b)
    },
    error = function(e) {
      # FAILURE: Return the robust NA structure with boot ID
      # This ENSURES the factor levels are maintained for bind_rows
      failure_tibble_structure %>%
        mutate(boot = b)
    }
  )
  
  # IMPORTANT CHECK: If a successful run still didn't produce all 9 combinations (e.g., NT vanished)
  # Fill in the missing rows with NAs to prevent loss of factor levels during summarise.
  if (nrow(boot_baci) < 9) {
    expected_combinations <- failure_tibble_structure %>% select(Order.q, Treatment)
    missing_rows <- anti_join(expected_combinations, boot_baci, by = c("Order.q", "Treatment"))
    
    if (nrow(missing_rows) > 0) {
      missing_rows_full <- missing_rows %>% 
        mutate(Estimate = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_, boot = b)
      boot_baci <- bind_rows(boot_baci, missing_rows_full)
    }
  }
  
  boot_list[[b]] <- boot_baci
}

boot_baci <- bind_rows(boot_list) %>%
  filter(!is.na(Estimate))    # Now this only drops the rows where the calculation failed (Estimate = NA)

## 2. Summarise bootstrap distribution per Treatment × q ----------------

boot_summary_UD <- boot_baci %>%
  group_by(Order.q, Treatment) %>%
  summarise(
    Estimate_boot_mean = mean(Estimate, na.rm = TRUE),
    LCL_boot           = quantile(Estimate, probs = 0.025, na.rm = TRUE),
    UCL_boot           = quantile(Estimate, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR","LR","HR")), # Order for plot Y-axis
    Order.q   = factor(Order.q, levels = c("q = 0","q = 1","q = 2"))
  )

# Forest plot of bootstrap BACI CIs 
p_baci_UD_boot <- ggplot(boot_summary_UD,
                         aes(x = Estimate_boot_mean,
                             y = Treatment,
                             color = Treatment)) + # Re-added color aesthetic for legend/scales
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL_boot, xmax = UCL_boot), height = 0.2, color = "black", linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_color_brewer(palette = "Dark2") + # Re-added scale_color_brewer
  labs(
    title = "Gamma Diversity (Unlimited Distance)",
    x = "Change in Effective # of Species",
    y = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right" # Changed to 'right' to show color mapping
  )

print(p_baci_UD_boot)

# Save
save_plot("Figures/RE plots/BACI_bootstrap_UD.png",
          p_baci_UD_boot, base_width = 7, base_height = 5)

## --- BACI CONTRAST TABLE CODE ---

baci_table_UD <- boot_summary_UD %>%
  mutate(
    CI_width = UCL_boot - LCL_boot,
    Sig = case_when(
      LCL_boot > 0  ~ "Increase (CI > 0)",
      UCL_boot < 0  ~ "Decrease (CI < 0)",
      TRUE          ~ "NS (includes 0)"
    )
  ) %>%
  # Order by q, then HR, LR, FR for the table presentation
  arrange(Order.q, factor(Treatment, levels = c("HR","LR","FR"))) %>% 
  transmute(
    `Diversity order` = Order.q,
    Treatment,
    Estimate = round(Estimate_boot_mean, 2),
    `LCL (95% CI)` = round(LCL_boot, 2),
    `UCL (95% CI)` = round(UCL_boot, 2),
    `CI width`     = round(CI_width, 2),
    `Interpretation` = Sig
  )

baci_table_UD


############# GOOD END ############

































































# ---
# title: "Extrapolation and rarefaction of Hill numbers (q=0, 1, 2)"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script create extrapolation and rarefaction curves for Hill numbers of the order q = 0, 1, and 2
# ---

# Clear environment
rm(list = ls())  # Removes all objects from the environment

# Install and load packages
#install.packages("iNEXT")
library(tidyverse)
library(iNEXT)
library(cowplot)
library(forcats)
library(patchwork)


####### UNLIMITED DISTANCE

# Load data
abundance_all_years_UDPC <- read.csv("Output/Tabular Data/max_count_all_years_UDPC.csv")

# Check structure
str(abundance_all_years_UDPC)

# ---- Identify species columns ----
species_cols <- names(abundance_all_years_UDPC)[grepl("^[A-Z]{4}$", names(abundance_all_years_UDPC))]

# Coerce species columns to numeric and replace NAs with 0
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)) )) %>%
  mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))

# Create treatment columns
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location))

# Rearrange columns
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  select(location, treatment, year, all_of(species_cols))

# Split into a dataframe each for 2023, 2024, and 2025
abund_2023_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2023")
abund_2024_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2024")
abund_2025_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2025")




### Create list of abundances for each year

# Sum species by treatment 
iNext_2023_UDPC <- abund_2023_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
iNext_2024_UDPC <- abund_2024_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
iNext_2025_UDPC <- abund_2025_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# Convert each treatment row to an unnamed numeric vector and deframe to a named list
iNext_2023_UDPC <- iNext_2023_UDPC %>%
  tidyr::nest(.by = treatment) %>%
  mutate(vec = map(data, ~ .x %>%
                     select(all_of(species_cols)) %>%
                     unlist(use.names = FALSE) %>%
                     as.numeric())) %>%
  transmute(treatment, vec) %>%
  deframe()

iNext_2024_UDPC <- iNext_2024_UDPC %>%
  tidyr::nest(.by = treatment) %>%
  mutate(vec = map(data, ~ .x %>%
                     select(all_of(species_cols)) %>%
                     unlist(use.names = FALSE) %>%
                     as.numeric())) %>%
  transmute(treatment, vec) %>%
  deframe()

iNext_2025_UDPC <- iNext_2025_UDPC %>%
  tidyr::nest(.by = treatment) %>%
  mutate(vec = map(data, ~ .x %>%
                     select(all_of(species_cols)) %>%
                     unlist(use.names = FALSE) %>%
                     as.numeric())) %>%
  transmute(treatment, vec) %>%
  deframe()

# Remove 0 values
iNext_2023_UDPC <- lapply(iNext_2023_UDPC, \(x) x[x > 0])
iNext_2024_UDPC <- lapply(iNext_2024_UDPC, \(x) x[x > 0])
iNext_2025_UDPC <- lapply(iNext_2025_UDPC, \(x) x[x > 0])






### R/E model

# Run R/E ~ abundance curves for all years
inext_2023_UDPC <- iNEXT(iNext_2023_UDPC, datatype = "abundance", q = c(0, 1, 2))
inext_2024_UDPC <- iNEXT(iNext_2024_UDPC, datatype = "abundance", q = c(0, 1, 2))
inext_2025_UDPC <- iNEXT(iNext_2025_UDPC, datatype = "abundance", q = c(0, 1, 2))

# Plot
RE_abund_2023_UDPC <- ggiNEXT(inext_2023_UDPC, type=1, facet.var="Order.q") + facet_wrap(~Order.q, scales="free")
RE_abund_2024_UDPC <- ggiNEXT(inext_2024_UDPC, type=1, facet.var="Order.q") + facet_wrap(~Order.q, scales="free")
RE_abund_2025_UDPC <- ggiNEXT(inext_2025_UDPC, type=1, facet.var="Order.q") + facet_wrap(~Order.q, scales="free")

# Save plots
save_plot("Figures/RE plots/RE_plot_2023_UDPC.png", RE_abund_2023_UDPC, base_width = 7, base_height = 5)
save_plot("Figures/RE plots/RE_plot_2024_UDPC.png", RE_abund_2024_UDPC, base_width = 7, base_height = 5)
save_plot("Figures/RE plots/RE_plot_2025_UDPC.png", RE_abund_2025_UDPC, base_width = 7, base_height = 5)



### R/E Hill Numbers ~ COVERAGE

ggiNEXT(inext_2023_UDPC, type=3, facet.var="Order.q")
ggiNEXT(inext_2024_UDPC, type=3, facet.var="Order.q")
ggiNEXT(inext_2025_UDPC, type=3, facet.var="Order.q")


################# Point-based coverage estimation

# Check for minimum coverage across years
estimateD(iNext_2023_UDPC, datatype="abundance", base="coverage", conf=0.95)
estimateD(iNext_2024_UDPC, datatype="abundance", base="coverage", conf=0.95)
estimateD(iNext_2025_UDPC, datatype="abundance", base="coverage", conf=0.95)
# Lower coverage = 0.982

# Now run
equal_cov_2023_UDPC <- estimateD(iNext_2023_UDPC, datatype="abundance", base="coverage", level = 0.982, conf=0.95)
equal_cov_2024_UDPC <- estimateD(iNext_2024_UDPC, datatype="abundance", base="coverage", level = 0.982, conf=0.95)
equal_cov_2025_UDPC <- estimateD(iNext_2025_UDPC, datatype="abundance", base="coverage", level = 0.982, conf=0.95)

######### Plot PB coverage estimate w/ 95% CI

# Add year to each df so they can be bound
equal_cov_2023_UDPC <- equal_cov_2023_UDPC %>%
  mutate(year = "2023")
equal_cov_2024_UDPC <- equal_cov_2024_UDPC %>%
  mutate(year = "2024")
equal_cov_2025_UDPC <- equal_cov_2025_UDPC %>%
  mutate(year = "2025")

# Combine into single df
equal_cov_all_years_UDPC <- bind_rows(equal_cov_2023_UDPC, equal_cov_2024_UDPC, equal_cov_2025_UDPC)





############ BACI contrasts UNPOOLED


# --- BACI contrasts per year (no pooling) -------------------------------------

df0_UDPC <- equal_cov_all_years_UDPC %>%
  select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, year) %>%
  mutate(
    Assemblage = factor(Assemblage, levels = c("NT","LR","HR","FR")),
    year       = factor(as.character(year), levels = c("2023","2024","2025")),
    Order.q    = as.integer(Order.q),
    SE         = (qD.UCL - qD.LCL) / (2*1.96),
    Var        = SE^2
  )

before_year <- "2023"
after_years <- c("2024","2025")
treatments  <- c("LR","HR","FR")
control     <- "NT"

# --- Compute per-year BACI contrasts ------------------------------------------
baci_results_UDPC <- map_dfr(after_years, function(ayr) {
  wide_UDPC <- df0_UDPC %>%
    filter(year %in% c(before_year, ayr)) %>%
    select(Assemblage, Order.q, year, qD, Var) %>%
    pivot_wider(names_from = year, values_from = c(qD, Var), names_sep = "_") %>%
    mutate(
      d_Treat = .data[[paste0("qD_", ayr)]] - .data[[paste0("qD_", before_year)]],
      v_Treat = .data[[paste0("Var_", ayr)]] + .data[[paste0("Var_", before_year)]]
    )
  
  ctrl <- wide_UDPC %>% filter(Assemblage == control) %>%
    select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)
  
  wide_UDPC %>%
    filter(Assemblage %in% treatments) %>%
    left_join(ctrl, by = "Order.q") %>%
    transmute(
      Order.q,
      AfterYear = ayr,
      Treatment = Assemblage,
      Estimate  = d_Treat - d_Ctrl,
      SE        = sqrt(v_Treat + v_Ctrl),
      LCL       = Estimate - 1.96*SE,
      UCL       = Estimate + 1.96*SE
    )
}) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("LR","HR","FR")),
    Order.q   = factor(Order.q, levels = c(0,1,2), labels = c("q = 0","q = 1","q = 2")),
    AfterYear = factor(AfterYear, levels = c("2024","2025"))
  )

# --- Plot: BACI contrasts by year ---------------------------------------------
p_baci_UDPC <- ggplot(baci_results_UDPC, aes(x = Estimate, y = Treatment, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.15, color = "black", linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(Order.q ~ AfterYear, scales = "free_x") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Unlimited Distance",
    x = "Change in Hill Numbers relative to Pre-Treatment",
    y = "Treatment",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_baci_UDPC)


# Save
save_plot("Figures/RE plots/BACI_contrasts_UDPC.png", p_baci_UDPC, base_width = 7, base_height = 5)




















###### GOOD START ##########


################################### POOLING2024/2025 w/ BOOTSTRAPPING

rm(list = ls())  # Removes all objects from the environment

# Install and load packages
#install.packages("iNEXT")
library(tidyverse)
library(iNEXT)
library(cowplot)
library(forcats)
library(patchwork)

# Load data
abundance_all_years_UDPC <- read.csv("Output/Tabular Data/max_count_all_years_UDPC.csv")

# Check structure
str(abundance_all_years_UDPC)

# ---- Identify species columns ----
species_cols <- names(abundance_all_years_UDPC)[grepl("^[A-Z]{4}$", names(abundance_all_years_UDPC))]

# Coerce species columns to numeric and replace NAs with 0
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))

# Create treatment column
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  mutate(treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location))

# Rearrange columns
abundance_all_years_UDPC <- abundance_all_years_UDPC %>%
  select(location, treatment, year, all_of(species_cols))

# Split into a dataframe each for 2023, 2024, and 2025
abund_2023_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2023")
abund_2024_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2024")
abund_2025_UDPC <- abundance_all_years_UDPC %>%
  filter(year == "2025")

## 1. Sum species per treatment per year -------------------------------

sum_2023_UD <- abund_2023_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

sum_2024_UD <- abund_2024_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

sum_2025_UD <- abund_2025_UDPC %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

## 2. Pooled AFTER = 2024 + 2025 per treatment ------------------------

sum_after_UD <- bind_rows(sum_2024_UD, sum_2025_UD) %>%
  group_by(treatment) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

## 3. Convert to iNEXT abundance lists --------------------------------

# BEFORE: 2023
iNext_before_UD <- sum_2023_UD %>%
  tidyr::nest(.by = treatment) %>%
  mutate(
    vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )
  ) %>%
  transmute(treatment, vec) %>%
  deframe()

# AFTER: pooled 2024+2025
iNext_after_UD <- sum_after_UD %>%
  tidyr::nest(.by = treatment) %>%
  mutate(
    vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )
  ) %>%
  transmute(treatment, vec) %>%
  deframe()

# Drop zero-abundance species 
iNext_before_UD <- lapply(iNext_before_UD, \(x) x[x > 0])
iNext_after_UD  <- lapply(iNext_after_UD,  \(x) x[x > 0])

## Optional iNEXT curves
inext_before_UD <- iNEXT(iNext_before_UD, datatype = "abundance", q = c(0, 1, 2))
inext_after_UD  <- iNEXT(iNext_after_UD,  datatype = "abundance", q = c(0, 1, 2))

# Plot (rarefaction/extrapolation)
RE_abund_before_UD <- ggiNEXT(inext_before_UD, type = 1, facet.var = "Order.q")
RE_abund_after_UD  <- ggiNEXT(inext_after_UD,  type = 1, facet.var = "Order.q")

# Save
save_plot("Figures/RE plots/RE_plot_before_UD.png",
          RE_abund_before_UD, base_width = 7, base_height = 5)
save_plot("Figures/RE plots/RE_plot_after_UD.png",
          RE_abund_after_UD,  base_width = 7, base_height = 5)

## 4. Find common sample size (min n across periods) ------------------

info_before <- DataInfo(iNext_before_UD, datatype = "abundance")
info_after  <- DataInfo(iNext_after_UD,  datatype = "abundance")

# n = total individuals (sum of abundances) in each assemblage
min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
min_n  # inspect

## 5. Standardize Hill numbers at this common sample size --------------

equal_before_UD <- estimateD(
  iNext_before_UD,
  datatype = "abundance",
  base     = "size",
  level    = min_n,
  conf     = 0.95
) %>%
  mutate(Period = "Before")

equal_after_UD <- estimateD(
  iNext_after_UD,
  datatype = "abundance",
  base     = "size",
  level    = min_n,
  conf     = 0.95
) %>%
  mutate(Period = "After")

equal_size_UD <- bind_rows(equal_before_UD, equal_after_UD)

## 6. BACI contrasts: (After-Before)_Treat - (After-Before)_NT ---------

df0_UD <- equal_size_UD %>%
  select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, Period) %>%
  mutate(
    Assemblage = factor(Assemblage, levels = c("NT", "LR", "HR", "FR")),
    Period     = factor(Period, levels = c("Before", "After")),
    Order.q    = as.integer(Order.q),
    SE         = (qD.UCL - qD.LCL) / (2 * 1.96),
    Var        = SE^2
  )

treatments <- c("LR", "HR", "FR")
control    <- "NT"

wide_UD <- df0_UD %>%
  select(Assemblage, Order.q, Period, qD, Var) %>%
  distinct() %>%  # guard against accidental duplicates
  tidyr::pivot_wider(
    names_from  = Period,
    values_from = c(qD, Var),
    names_sep   = "_"
  ) %>%
  mutate(
    d_Treat = qD_After - qD_Before,
    v_Treat = Var_After + Var_Before
  )

ctrl_UD <- wide_UD %>%
  filter(Assemblage == control) %>%
  select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)

baci_results_UD <- wide_UD %>%
  filter(Assemblage %in% treatments) %>%
  left_join(ctrl_UD, by = "Order.q", relationship = "many-to-one") %>%
  transmute(
    Order.q  = factor(Order.q, levels = c(0, 1, 2),
                      labels = c("q = 0", "q = 1", "q = 2")),
    Treatment = factor(Assemblage, levels = c("FR", "LR", "HR")),
    Estimate  = d_Treat - d_Ctrl,
    SE        = sqrt(v_Treat + v_Ctrl),
    LCL       = Estimate - 1.96 * SE,
    UCL       = Estimate + 1.96 * SE
  )


## 7. Plot -------------------------------------------------------------

p_baci_UD <- ggplot(baci_results_UD,
                    aes(x = Estimate, y = Treatment, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 height = 0.15, color = "black") +
  geom_point(size = 3) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Limited Amplitude (2024/2025 pooled)",
    x = "Change in Hill Numbers",
    y = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "right",
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold")
  )

print(p_baci_UD)

save_plot("Figures/RE plots/BACI_contrasts_UD_pooled_size.png",
          p_baci_UD, base_width = 7, base_height = 5)

























############ GEMINI SIMPLIFICATION


## =====================================================================
## BOOTSTRAP BACI CONTRASTS (UNLIMITED DISTANCE, UDPC)
## =====================================================================

library(tidyverse)
library(iNEXT)
library(ggplot2) # Include plotting package

# ---------------------------------------------------------------------
# 1. Helper: compute BACI contrasts from a *given* raw data frame
#    (The definition of this function remains UNCHANGED from your working code)
# ---------------------------------------------------------------------

compute_baci_udpc <- function(dat) {
  # Identify species columns
  species_cols <- names(dat)[grepl("^[A-Z]{4}$", names(dat))]
  
  # Make sure species are numeric and NAs are zero
  dat <- dat %>%
    mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))
  
  # Split by year
  abund_2023 <- dat %>% filter(year == "2023")
  abund_2024 <- dat %>% filter(year == "2024")
  abund_2025 <- dat %>% filter(year == "2025")
  
  # If any period is empty, bail out with empty tibble
  if (nrow(abund_2023) == 0 || (nrow(abund_2024) == 0 & nrow(abund_2025) == 0)) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Sum species per treatment per year ----
  sum_2023 <- abund_2023 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  sum_2024 <- abund_2024 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  sum_2025 <- abund_2025 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  ## ---- Pooled AFTER = 2024 + 2025 per treatment ----
  sum_after <- bind_rows(sum_2024, sum_2025) %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
  
  # If some treatment disappears completely, bail
  if (nrow(sum_2023) == 0 || nrow(sum_after) == 0) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Build iNEXT abundance lists (Before/After) ----
  iNext_before <- sum_2023 %>%
    tidyr::nest(.by = treatment) %>%
    mutate(vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )) %>%
    transmute(treatment, vec) %>%
    deframe()
  
  iNext_after <- sum_after %>%
    tidyr::nest(.by = treatment) %>%
    mutate(vec = purrr::map(
      data,
      ~ .x %>%
        select(all_of(species_cols)) %>%
        unlist(use.names = FALSE) %>%
        as.numeric()
    )) %>%
    transmute(treatment, vec) %>%
    deframe()
  
  # Drop zero-abundance species
  iNext_before <- lapply(iNext_before, \(x) x[x > 0])
  iNext_after  <- lapply(iNext_after,  \(x) x[x > 0])
  
  ## ---- Common sample size (min n across periods) ----
  info_before <- DataInfo(iNext_before, datatype = "abundance")
  info_after  <- DataInfo(iNext_after,  datatype = "abundance")
  
  min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
  
  # If min_n is invalid, bail
  if (!is.finite(min_n) || min_n <= 0) {
    return(tibble(
      Order.q   = factor(),
      Treatment = factor(),
      Estimate  = numeric(),
      SE        = numeric(),
      LCL       = numeric(),
      UCL       = numeric()
    ))
  }
  
  ## ---- Standardize Hill numbers at common sample size ----
  equal_before <- estimateD(
    iNext_before,
    datatype = "abundance",
    base     = "size",
    level    = min_n,
    conf     = 0.95
  ) %>%
    mutate(Period = "Before")
  
  equal_after <- estimateD(
    iNext_after,
    datatype = "abundance",
    base     = "size",
    level    = min_n,
    conf     = 0.95
  ) %>%
    mutate(Period = "After")
  
  equal_size <- bind_rows(equal_before, equal_after)
  
  ## ---- BACI contrasts: (After-Before)_Treat - (After-Before)_NT ----
  
  df0 <- equal_size %>%
    select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, Period) %>%
    mutate(
      Assemblage = factor(Assemblage, levels = c("NT", "LR", "HR", "FR")),
      Period     = factor(Period, levels = c("Before", "After")),
      Order.q    = as.integer(Order.q),
      SE         = (qD.UCL - qD.LCL) / (2 * 1.96),
      Var        = SE^2
    )
  
  treatments <- c("LR", "HR", "FR")
  control    <- "NT"
  
  wide <- df0 %>%
    select(Assemblage, Order.q, Period, qD, Var) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from  = Period,
      values_from = c(qD, Var),
      names_sep   = "_"
    ) %>%
    mutate(
      d_Treat = qD_After - qD_Before,
      v_Treat = Var_After + Var_Before
    )
  
  ctrl <- wide %>%
    filter(Assemblage == control) %>%
    select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)
  
  baci_results <- wide %>%
    filter(Assemblage %in% treatments) %>%
    left_join(ctrl, by = "Order.q", relationship = "many-to-one") %>%
    transmute(
      Order.q   = factor(Order.q, levels = c(0, 1, 2),
                         labels = c("q = 0", "q = 1", "q = 2")),
      Treatment = factor(Assemblage, levels = c("FR", "LR", "HR")),
      Estimate  = d_Treat - d_Ctrl,
      SE        = sqrt(v_Treat + v_Ctrl),
      LCL       = Estimate - 1.96 * SE,
      UCL       = Estimate + 1.96 * SE
    )
  
  return(baci_results)
}

## =====================================================================
## BOOTSTRAP BACI CONTRASTS (UNLIMITED DISTANCE, UDPC) – BOOTSTRAP PLOT
## =====================================================================

## 1. Bootstrap: resample locations within treatment, keep all years ----

set.seed(123)    # for reproducibility
B <- 1000        # number of bootstrap replicates

# Sampling units: locations within each treatment
loc_tbl <- abundance_all_years_UDPC %>%
  distinct(location, treatment)

boot_list <- vector("list", B)

# Define the structure for a failed run (9 rows: q=0,1,2 x T=FR,LR,HR)
failure_tibble_structure <- tibble(
  Order.q = factor(rep(c("q = 0", "q = 1", "q = 2"), each = 3), levels = c("q = 0", "q = 1", "q = 2")),
  Treatment = factor(rep(c("FR", "LR", "HR"), times = 3), levels = c("FR", "LR", "HR")),
  Estimate = NA_real_,
  SE = NA_real_,
  LCL = NA_real_,
  UCL = NA_real_
)


for (b in seq_len(B)) {
  # Resample locations *within* each treatment, with replacement
  sampled_locs <- loc_tbl %>%
    group_by(treatment) %>%
    summarise(
      location = sample(location, size = n(), replace = TRUE),
      .groups  = "drop"
    )
  
  # Join back to full data: keeps all years for each selected location
  boot_dat <- sampled_locs %>%
    left_join(abundance_all_years_UDPC,
              by = c("location", "treatment"))
  
  # Compute BACI for this bootstrap sample, with error handling
  boot_baci <- tryCatch(
    {
      # SUCCESS: Calculate BACI and add boot ID
      compute_baci_udpc(boot_dat) %>%
        mutate(boot = b)
    },
    error = function(e) {
      # FAILURE: Return the robust NA structure with boot ID
      # This ENSURES the factor levels are maintained for bind_rows
      failure_tibble_structure %>%
        mutate(boot = b)
    }
  )
  
  # IMPORTANT CHECK: If a successful run still didn't produce all 9 combinations (e.g., NT vanished)
  # Fill in the missing rows with NAs to prevent loss of factor levels during summarise.
  if (nrow(boot_baci) < 9) {
    expected_combinations <- failure_tibble_structure %>% select(Order.q, Treatment)
    missing_rows <- anti_join(expected_combinations, boot_baci, by = c("Order.q", "Treatment"))
    
    if (nrow(missing_rows) > 0) {
      missing_rows_full <- missing_rows %>% 
        mutate(Estimate = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_, boot = b)
      boot_baci <- bind_rows(boot_baci, missing_rows_full)
    }
  }
  
  boot_list[[b]] <- boot_baci
}

boot_baci <- bind_rows(boot_list) %>%
  filter(!is.na(Estimate))    # Now this only drops the rows where the calculation failed (Estimate = NA)

## 2. Summarise bootstrap distribution per Treatment × q ----------------

boot_summary_UD <- boot_baci %>%
  group_by(Order.q, Treatment) %>%
  summarise(
    Estimate_boot_mean = mean(Estimate, na.rm = TRUE),
    LCL_boot           = quantile(Estimate, probs = 0.025, na.rm = TRUE),
    UCL_boot           = quantile(Estimate, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR","LR","HR")), # Order for plot Y-axis
    Order.q   = factor(Order.q, levels = c("q = 0","q = 1","q = 2"))
  )

# Forest plot of bootstrap BACI CIs 
p_baci_UD_boot <- ggplot(boot_summary_UD,
                         aes(x = Estimate_boot_mean,
                             y = Treatment,
                             color = Treatment)) + # Re-added color aesthetic for legend/scales
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL_boot, xmax = UCL_boot), height = 0.2, color = "black", linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_color_brewer(palette = "Dark2") + # Re-added scale_color_brewer
  labs(
    title = "Gamma Diversity (Unlimited Distance)",
    x = "Change in Effective # of Species",
    y = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right" # Changed to 'right' to show color mapping
  )

print(p_baci_UD_boot)

# Save
save_plot("Figures/RE plots/BACI_bootstrap_UD.png",
          p_baci_UD_boot, base_width = 7, base_height = 5)

## --- BACI CONTRAST TABLE CODE ---

baci_table_UD <- boot_summary_UD %>%
  mutate(
    CI_width = UCL_boot - LCL_boot,
    Sig = case_when(
      LCL_boot > 0  ~ "Increase (CI > 0)",
      UCL_boot < 0  ~ "Decrease (CI < 0)",
      TRUE          ~ "NS (includes 0)"
    )
  ) %>%
  # Order by q, then HR, LR, FR for the table presentation
  arrange(Order.q, factor(Treatment, levels = c("HR","LR","FR"))) %>% 
  transmute(
    `Diversity order` = Order.q,
    Treatment,
    Estimate = round(Estimate_boot_mean, 2),
    `LCL (95% CI)` = round(LCL_boot, 2),
    `UCL (95% CI)` = round(UCL_boot, 2),
    `CI width`     = round(CI_width, 2),
    `Interpretation` = Sig
  )

baci_table_UD


############# GOOD END ############















############### 6b GAMMA DIV MAX COUNT



# Clear environment
rm(list = ls())  # Removes all objects from the environment

library(tidyverse)
library(iNEXT)
library(cowplot)
library(forcats)
library(patchwork)


# Load data 
abundance_all_years_LA <- read.csv("Output/Tabular Data/max_count_all_years_LA.csv")

# ---------------------------------------------------------------------
# Data Pre-processing (One-time setup)
# ---------------------------------------------------------------------

# ---- Identify species columns (4 uppercase letters) ----
species_cols <- names(abundance_all_years_LA)[grepl("^[A-Z]{4}$", names(abundance_all_years_LA))]

# Coerce species columns to numeric and replace NAs with 0
abundance_all_years_LA <- abundance_all_years_LA %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(species_cols), ~ tidyr::replace_na(.x, 0)))

# Create treatment column (extracting LR, HR, FR, NT from location name)
abundance_all_years_LA <- abundance_all_years_LA %>%
  mutate(treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location))

# Rearrange columns
abundance_all_years_LA <- abundance_all_years_LA %>%
  select(location, treatment, year, all_of(species_cols))

# ---------------------------------------------------------------------
# 1. Helper Function: Compute BACI contrasts for a single dataset (LA)
# ---------------------------------------------------------------------
# Note: This function uses the same logic (pooling 2024+2025) as the UDPC version
# but is renamed to fit the LA context.

compute_baci_LA <- function(dat) {
  # Identify species columns
  species_cols_func <- names(dat)[grepl("^[A-Z]{4}$", names(dat))]
  
  # Ensure data is clean (though done globally, safe to repeat in helper)
  dat <- dat %>%
    mutate(across(all_of(species_cols_func), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(species_cols_func), ~ tidyr::replace_na(.x, 0)))
  
  # Split by year
  abund_2023 <- dat %>% filter(year == "2023")
  abund_after <- dat %>% filter(year %in% c("2024", "2025"))
  
  # If any period is empty, bail out with empty tibble
  if (nrow(abund_2023) == 0 || nrow(abund_after) == 0) {
    return(tibble(Order.q = factor(), Treatment = factor(), Estimate = numeric(), 
                  SE = numeric(), LCL = numeric(), UCL = numeric()))
  }
  
  ## ---- Sum species per treatment per period (Pooling 2024+2025 in 'After') ----
  sum_2023 <- abund_2023 %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols_func), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  sum_after <- abund_after %>%
    group_by(treatment) %>%
    summarise(across(all_of(species_cols_func), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  # If some treatment disappears completely, bail
  if (nrow(sum_2023) == 0 || nrow(sum_after) == 0) {
    return(tibble(Order.q = factor(), Treatment = factor(), Estimate = numeric(), 
                  SE = numeric(), LCL = numeric(), UCL = numeric()))
  }
  
  ## ---- Build iNEXT abundance lists (Before/After) ----
  
  # Helper to convert summary table to iNEXT list format
  to_iNext_list <- function(df) {
    df %>%
      tidyr::nest(.by = treatment) %>%
      mutate(vec = purrr::map(
        data,
        ~ .x %>% select(all_of(species_cols_func)) %>% unlist(use.names = FALSE) %>% as.numeric()
      )) %>%
      transmute(treatment, vec) %>%
      deframe() %>%
      lapply(\(x) x[x > 0]) # Drop zero-abundance species
  }
  
  iNext_before <- to_iNext_list(sum_2023)
  iNext_after  <- to_iNext_list(sum_after)
  
  ## ---- Common sample size (min n across periods) ----
  info_before <- DataInfo(iNext_before, datatype = "abundance")
  info_after  <- DataInfo(iNext_after,  datatype = "abundance")
  
  min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
  
  if (!is.finite(min_n) || min_n <= 0) {
    return(tibble(Order.q = factor(), Treatment = factor(), Estimate = numeric(), 
                  SE = numeric(), LCL = numeric(), UCL = numeric()))
  }
  
  ## ---- Standardize Hill numbers at common sample size ----
  estimate_hill <- function(data_list, period) {
    estimateD(data_list, datatype = "abundance", base = "size", level = min_n, conf = 0.95) %>%
      mutate(Period = period)
  }
  
  equal_before <- estimate_hill(iNext_before, "Before")
  equal_after  <- estimate_hill(iNext_after, "After")
  
  equal_size <- bind_rows(equal_before, equal_after)
  
  ## ---- BACI contrasts: (After-Before)_Treat - (After-Before)_NT ----
  
  df0 <- equal_size %>%
    select(Assemblage, Order.q, qD, qD.LCL, qD.UCL, Period) %>%
    mutate(
      Order.q    = as.integer(Order.q),
      SE         = (qD.UCL - qD.LCL) / (2 * 1.96),
      Var        = SE^2
    )
  
  treatments <- c("LR", "HR", "FR")
  control    <- "NT"
  
  wide <- df0 %>%
    select(Assemblage, Order.q, Period, qD, Var) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from  = Period,
      values_from = c(qD, Var),
      names_sep   = "_"
    ) %>%
    mutate(
      d_Treat = qD_After - qD_Before,
      v_Treat = Var_After + Var_Before
    )
  
  ctrl <- wide %>%
    filter(Assemblage == control) %>%
    select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)
  
  baci_results <- wide %>%
    filter(Assemblage %in% treatments) %>%
    left_join(ctrl, by = "Order.q", relationship = "many-to-one") %>%
    transmute(
      Order.q   = factor(Order.q, levels = c(0, 1, 2), labels = c("q = 0", "q = 1", "q = 2")),
      Treatment = factor(Assemblage, levels = c("FR", "LR", "HR")),
      Estimate  = d_Treat - d_Ctrl,
      SE        = sqrt(v_Treat + v_Ctrl),
      LCL       = Estimate - 1.96 * SE,
      UCL       = Estimate + 1.96 * SE
    )
  
  return(baci_results)
}


## =====================================================================
## 2. BOOTSTRAP EXECUTION (LA)
## =====================================================================

set.seed(123)    # for reproducibility
B <- 1000        # number of bootstrap replicates

# Sampling units: locations within each treatment
loc_tbl <- abundance_all_years_LA %>%
  distinct(location, treatment)

boot_list <- vector("list", B)

# Define the structure for a failed run (9 rows: q=0,1,2 x T=FR,LR,HR)
failure_tibble_structure <- tibble(
  Order.q = factor(rep(c("q = 0", "q = 1", "q = 2"), each = 3), levels = c("q = 0", "q = 1", "q = 2")),
  Treatment = factor(rep(c("FR", "LR", "HR"), times = 3), levels = c("FR", "LR", "HR")),
  Estimate = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_
)

for (b in seq_len(B)) {
  # Resample locations *within* each treatment, with replacement
  sampled_locs <- loc_tbl %>%
    group_by(treatment) %>%
    summarise(location = sample(location, size = n(), replace = TRUE), .groups  = "drop")
  
  # Join back to full data: keeps all years for each selected location
  boot_dat <- sampled_locs %>%
    left_join(abundance_all_years_LA, by = c("location", "treatment"))
  
  # Compute BACI for this bootstrap sample, with robust error handling
  boot_baci <- tryCatch(
    {
      compute_baci_LA(boot_dat) %>% mutate(boot = b)
    },
    error = function(e) {
      # FAILURE: Return the robust NA structure with boot ID
      failure_tibble_structure %>% mutate(boot = b)
    }
  )
  
  # Check for missing combinations in successful run (if not 9 rows) and fill with NA
  if (nrow(boot_baci) < 9) {
    expected_combinations <- failure_tibble_structure %>% select(Order.q, Treatment)
    missing_rows <- anti_join(expected_combinations, boot_baci, by = c("Order.q", "Treatment"))
    
    if (nrow(missing_rows) > 0) {
      missing_rows_full <- missing_rows %>% 
        mutate(Estimate = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_, boot = b)
      boot_baci <- bind_rows(boot_baci, missing_rows_full)
    }
  }
  
  boot_list[[b]] <- boot_baci
}

# Combine and filter out replicates that failed (Estimate is NA)
boot_baci <- bind_rows(boot_list) %>%
  filter(!is.na(Estimate))

## =====================================================================
## 3. SUMMARISE AND PLOT
## =====================================================================

# Summarise bootstrap distribution per Treatment × q 
boot_summary_LA <- boot_baci %>%
  group_by(Order.q, Treatment) %>%
  summarise(
    Estimate_boot_mean = mean(Estimate, na.rm = TRUE),
    LCL_boot           = quantile(Estimate, probs = 0.025, na.rm = TRUE),
    UCL_boot           = quantile(Estimate, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR","LR","HR")), # Order for plot Y-axis
    Order.q   = factor(Order.q, levels = c("q = 0","q = 1","q = 2"))
  )

# Forest plot of bootstrap BACI CIs 
p_baci_LA_boot <- ggplot(boot_summary_LA,
                         aes(x = Estimate_boot_mean,
                             y = Treatment,
                             color = Treatment)) + # Re-added color aesthetic for legend/scales
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL_boot, xmax = UCL_boot), height = 0.2, color = "black", linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_color_brewer(palette = "Dark2") + # Re-added scale_color_brewer
  labs(
    title = "Gamma Diversity (Limited Amplitude)",
    x = "Change in Effective # of Species",
    y = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right" # Changed to 'right' to show color mapping
  )

print(p_baci_LA_boot)

# Save
save_plot("Figures/RE plots/BACI_bootstrap_LA.png",
          p_baci_LA_boot, base_width = 7, base_height = 5)

## --- BACI CONTRAST TABLE CODE ---

baci_table_LA <- boot_summary_LA %>%
  mutate(
    CI_width = UCL_boot - LCL_boot,
    Sig = case_when(
      LCL_boot > 0  ~ "Increase (CI > 0)",
      UCL_boot < 0  ~ "Decrease (CI < 0)",
      TRUE          ~ "NS (includes 0)"
    )
  ) %>%
  # Order by q, then HR, LR, FR for the table presentation
  arrange(Order.q, factor(Treatment, levels = c("HR","LR","FR"))) %>% 
  transmute(
    `Diversity order` = Order.q,
    Treatment,
    Estimate = round(Estimate_boot_mean, 2),
    `LCL (95% CI)` = round(LCL_boot, 2),
    `UCL (95% CI)` = round(UCL_boot, 2),
    `CI width`     = round(CI_width, 2),
    `Interpretation` = Sig
  )

baci_table_LA


########## GOOD END ############
















############# 8a alpha diversity - MAX COUNT - UD


# ---
# title: "Alpha diversity using UD"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: 
# ---

## ==========================================
## 0. Packages
## ==========================================
library(dplyr)
library(tidyr)
library(vegan)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(cowplot)   # for save_plot if you use it elsewhere
library(DHARMa)   # needed for simulateResiduals()

# Data 
dat1_alp_UD <- read.csv("Output/Tabular Data/max_count_all_years_UDPC.csv")

# Identify species columns in the filtered dataset
species_cols <- setdiff(names(dat1_alp_UD), c("location", "year"))

# Build species matrix
sp_mat <- dat1_alp_UD %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# q = 0: richness
alpha_q0 <- specnumber(sp_mat)

# q = 1: exp(Shannon)
alpha_shannon <- diversity(sp_mat, index = "shannon")
alpha_q1 <- exp(alpha_shannon)

# q = 2: inverse Simpson (Hill q = 2)
alpha_q2 <- diversity(sp_mat, index = "invsimpson")

# Attach alpha metrics to the data
alpha_dat <- dat1_alp_UD %>%
  mutate(
    alpha_q0 = alpha_q0,
    alpha_q1 = alpha_q1,
    alpha_q2 = alpha_q2
  )


## Create site, block, treatment, time_period (Before vs After only)
alpha_dat <- alpha_dat %>%
  mutate(
    # site = location
    site = location,
    
    # block = prefix before first "-"
    block = sub("-.*$", "", location),
    
    # treatment = middle bit between first and second "-"
    # e.g., DOC213-FR-C -> "FR"
    treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location),
    
    # time period from year: 2023 = Before; 2024 & 2025 = After
    time_period = case_when(
      year == 2023 ~ "Before",
      year %in% c(2024, 2025) ~ "After",
      TRUE ~ NA_character_
    ),
    
    time_period = factor(time_period, levels = c("Before", "After")),
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR"))
  )


######### Fit GLMMs for alpha Hill numbers

### q = 0 (alpha richness)
alpha_q0_model_UD <- glmmTMB(
  alpha_q0 ~ treatment  + time_period + treatment * time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q0_model_UD)

# Assess model residuals
set.seed(123)  # for reproducibility of the simulations
alpha_q0_res <- simulateResiduals(
  fittedModel = alpha_q0_model_UD,
  n = 1000
)
plot(alpha_q0_res)



### q = 1
alpha_q1_model_UD <- glmmTMB(
  alpha_q1 ~ treatment * time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q1_model_UD)

# Assess model residuals
set.seed(123)
alpha_q1_res <- simulateResiduals(
  fittedModel = alpha_q1_model_UD,
  n = 1000
)
plot(alpha_q1_res)



### q = 2
alpha_q2_model_UD <- glmmTMB(
  alpha_q2 ~ treatment * time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q2_model_UD)

# Assess model residuals
set.seed(123)
alpha_q2_res <- simulateResiduals(
  fittedModel = alpha_q2_model_UD,
  n = 1000
)
plot(alpha_q2_res)



## ==========================================
## 4. Get EMMs (expected alpha Hill numbers)
##    for treatment * time_period
## ==========================================

# Cell means on the linear scale (gaussian link is identity)
emm_q0_UD <- emmeans(alpha_q0_model_UD, ~ treatment * time_period)
emm_q1_UD <- emmeans(alpha_q1_model_UD, ~ treatment * time_period)
emm_q2_UD <- emmeans(alpha_q2_model_UD, ~ treatment * time_period)


## ==========================================
## 5. BACI-style contrasts via emmeans::contrast
##    (Using original, mathematically verified coefficient matrix)
## ==========================================

get_baci_UD <- function(emm_obj, order_q) {
  
  # Coefficient vectors defined for the confirmed row order:
  # 1 NT Before, 2 LR Before, 3 HR Before, 4 FR Before, 
  # 5 NT After, 6 LR After, 7 HR After, 8 FR After
  L <- list(
    LR = c( 1, -1,  0,  0, -1,  1,  0,  0),
    HR = c( 1,  0, -1,  0, -1,  0,  1,  0),
    FR = c( 1,  0,  0, -1, -1,  0,  0,  1)
  )
  
  # Linear contrasts
  out <- contrast(emm_obj, method = L)
  
  # Get CIs and convert to a data frame. 
  # The diagnostic confirmed this output has 'lower.CL' and 'upper.CL'.
  out_df <- as.data.frame(confint(out))
  
  # Format the output using the definitive column names.
  out_df %>%
    dplyr::transmute(
      Order.q   = order_q,
      Treatment = contrast,
      Estimate  = estimate,
      SE        = SE,
      # These are the correct, confirmed column names
      LCL       = lower.CL, 
      UCL       = upper.CL
    )
}

# Apply to each q
baci_q0_UD <- get_baci_UD(emm_q0_UD, 0)
baci_q1_UD <- get_baci_UD(emm_q1_UD, 1)
baci_q2_UD <- get_baci_UD(emm_q2_UD, 2)

# Combine and set factor orders for plotting
baci_results_alpha_UD <- bind_rows(baci_q0_UD, baci_q1_UD, baci_q2_UD) %>%
  mutate(
    Order.q = factor(
      Order.q,
      levels = c(0, 1, 2),
      labels = c("q = 0", "q = 1", "q = 2")
    ),
    # HR at top, LR middle, FR bottom
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR"))
  )



## ==========================================
## 6. Plot BACI contrasts (alpha Hill numbers)
## ==========================================

p_baci_alpha_UD <- ggplot(baci_results_alpha_UD,
                          aes(x = Estimate, y = Treatment, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 height = 0.15, color = "black", linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  # Legend in same order: HR, LR, FR
  scale_color_brewer(palette = "Dark2",
                     breaks = c("HR", "LR", "FR")) +
  labs(
    title = "Alpha Diversity (Unlimited Distance)",
    x = "Change in Effective # of Species",
    y = "Treatment",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "right",
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold")
  )

print(p_baci_alpha_UD)












############### 8b alpha div LA

# ---
# title: "Alpha diversity using LA"
# author: "Leonard Patterson"
# created: "2025-07-02"
# description: 
# ---

## ==========================================
## 0. Packages
## ==========================================
library(dplyr)
library(tidyr)
library(vegan)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(cowplot)   # for save_plot if you use it elsewhere
library(DHARMa)   # needed for simulateResiduals()

# Data 
dat1_alp_LA <- read.csv("Output/Tabular Data/max_count_all_years_LA.csv")

# Identify species columns in the filtered dataset
species_cols <- setdiff(names(dat1_alp_LA), c("location", "year"))

# Build species matrix
sp_mat <- dat1_alp_LA %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# q = 0: richness
alpha_q0 <- specnumber(sp_mat)

# q = 1: exp(Shannon)
alpha_shannon <- diversity(sp_mat, index = "shannon")
alpha_q1 <- exp(alpha_shannon)

# q = 2: inverse Simpson (Hill q = 2)
alpha_q2 <- diversity(sp_mat, index = "invsimpson")

# Attach alpha metrics to the data
alpha_dat <- dat1_alp_LA %>%
  mutate(
    alpha_q0 = alpha_q0,
    alpha_q1 = alpha_q1,
    alpha_q2 = alpha_q2
  )


## Create site, block, treatment, time_period (Before vs After only)

alpha_dat <- alpha_dat %>%
  mutate(
    # site = location (per your instruction)
    site = location,
    
    # block = prefix before first "-"
    block = sub("-.*$", "", location),
    
    # treatment = middle bit between first and second "-"
    # e.g., DOC213-FR-C -> "FR"
    treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location),
    
    # time period from year: 2023 = Before; 2024 & 2025 = After
    time_period = case_when(
      year == 2023 ~ "Before",
      year %in% c(2024, 2025) ~ "After",
      TRUE ~ NA_character_
    ),
    
    time_period = factor(time_period, levels = c("Before", "After")),
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR"))
  )

# Check data
head(alpha_dat)






######### Fit GLMMs for alpha Hill numbers

### q = 0 (alpha richness)
alpha_q0_model_LA <- glmmTMB(
  alpha_q0 ~ treatment + time_period + treatment*time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q0_model_LA)

# Assess model residuals
set.seed(123)  # for reproducibility of the simulations
alpha_q0_res <- simulateResiduals(
  fittedModel = alpha_q0_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(alpha_q0_res)

### q = 1
alpha_q1_model_LA <- glmmTMB(
  alpha_q1 ~ treatment + time_period + treatment*time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q1_model_LA)

# Assess model residuals
set.seed(123)
alpha_q1_res <- simulateResiduals(
  fittedModel = alpha_q1_model_LA,
  n = 1000
)
plot(alpha_q1_res)

### q = 2
alpha_q2_model_LA <- glmmTMB(
  alpha_q2 ~ treatment + time_period + treatment*time_period + (1 | site) + (1 | block),
  family = gaussian(),
  data   = alpha_dat
)
summary(alpha_q2_model_LA)

# Assess model residuals
set.seed(123)
alpha_q2_res <- simulateResiduals(
  fittedModel = alpha_q2_model_LA,
  n = 1000
)
plot(alpha_q2_res)

hist(alpha_dat$alpha_q0)
hist(alpha_dat$alpha_q1)
hist(alpha_dat$alpha_q2)


## ==========================================
## 4. Get EMMs (expected alpha Hill numbers)
##    for treatment * time_period
## ==========================================

# Cell means on the linear scale (gaussian link is identity, so this is fine)
emm_q0_LA <- emmeans(alpha_q0_model_LA, ~ treatment * time_period)
emm_q1_LA <- emmeans(alpha_q1_model_LA, ~ treatment * time_period)
emm_q2_LA <- emmeans(alpha_q2_model_LA, ~ treatment * time_period)




## ==========================================
## 5. BACI-style contrasts via emmeans::contrast
##    (Using original, mathematically verified coefficient matrix)
## ==========================================

get_baci_LA <- function(emm_obj, order_q) {
  
  # Coefficient vectors defined for the confirmed row order:
  # 1 NT Before, 2 LR Before, 3 HR Before, 4 FR Before, 
  # 5 NT After, 6 LR After, 7 HR After, 8 FR After
  L <- list(
    LR = c( 1, -1,  0,  0, -1,  1,  0,  0),
    HR = c( 1,  0, -1,  0, -1,  0,  1,  0),
    FR = c( 1,  0,  0, -1, -1,  0,  0,  1)
  )
  
  # Linear contrasts
  out <- contrast(emm_obj, method = L)
  
  # Get CIs and convert to a data frame. 
  # The diagnostic confirmed this output has 'lower.CL' and 'upper.CL'.
  out_df <- as.data.frame(confint(out))
  
  # Format the output using the definitive column names.
  out_df %>%
    dplyr::transmute(
      Order.q   = order_q,
      Treatment = contrast,
      Estimate  = estimate,
      SE        = SE,
      # These are the correct, confirmed column names
      LCL       = lower.CL, 
      UCL       = upper.CL
    )
}

# Apply to each q
baci_q0_LA <- get_baci_LA(emm_q0_LA, 0)
baci_q1_LA <- get_baci_LA(emm_q1_LA, 1)
baci_q2_LA <- get_baci_LA(emm_q2_LA, 2)

# Combine and set factor orders for plotting
baci_results_alpha_LA <- bind_rows(baci_q0_LA, baci_q1_LA) %>%
  mutate(
    Order.q = factor(
      Order.q,
      levels = c(0, 1),
      labels = c("Species richness", "Shannon Diversity")
    ),
    # FR at bottom, LR middle, HR top in facet
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR"))
  )


## ==========================================
## 6. Plot BACI contrasts (alpha Hill numbers)
## ==========================================

# Compute symmetric x-axis limits around 0
x_lim <- max(abs(c(
  baci_results_alpha_LA$LCL,
  baci_results_alpha_LA$UCL
)), na.rm = TRUE)

p_baci_alpha_LA <- ggplot(
  baci_results_alpha_LA,
  aes(x = Estimate, y = Treatment, color = Treatment)
) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 1) +
  geom_errorbarh(
    aes(xmin = LCL, xmax = UCL),
    height = 0.15, color = "black", linewidth = 0.9
  ) +
  geom_point(size = 2.8) +
  facet_wrap(~ Order.q, ncol = 1, scales = "fixed") +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +   # 👈 key line
  scale_color_brewer(
    palette = "Dark2",
    breaks  = c("HR", "LR", "FR")
  ) +
  labs(
    title = "Alpha Diversity (Limited Amplitude)",
    x     = "Change in Effective # of Species",
    y     = "Treatment",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "right",
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank(),
    strip.text        = element_text(face = "bold")
  )

print(p_baci_alpha_LA)







########################### 6a MEAN COUNT



# ---
# title: "Extrapolation and rarefaction of Gamma Diversity Hill numbers (q=0, 1, 2) - UD"
# author: "Leonard Patterson"
# created: "2025-12-12"
# description: "This script create extrapolation and rarefaction curves for Hill numbers of the order q = 0, 1, and 2
# ---

# Clear environment
rm(list = ls())  # Removes all objects from the environment

# Install and load packages
#install.packages("iNEXT")
library(tidyverse)
library(iNEXT)
library(cowplot)
library(forcats)
library(patchwork)

# Load data 
abundance_all_years_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

# ---------------------
# Data Pre-processing 
# ---------------------

# ---- Identify species columns (4 uppercase letters) ----
species_cols <- names(abundance_all_years_UD)[grepl("^[A-Z]{4}$", names(abundance_all_years_UD))]

# Coerce species columns to numeric 
abundance_all_years_UD <- abundance_all_years_UD %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x))))

# OPTIONAL but recommended: check that species columns are integers
if (any(unlist(abundance_all_years_UD[species_cols]) %% 1 != 0, na.rm = TRUE)) {
  stop("Species columns contain non-integer values. iNEXT abundance data must be integer counts.")
}

# Make sure year is a character
abundance_all_years_UD <- abundance_all_years_UD %>%
  mutate(year = as.character(year))

# Create treatment column (extracting LR, HR, FR, NT from location name)
abundance_all_years_UD <- abundance_all_years_UD %>%
  mutate(treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location))

# Rearrange columns
abundance_all_years_UD <- abundance_all_years_UD %>%
  select(location, treatment, year, all_of(species_cols))

# Condense to one row per site-year by summing counts across surveys
abundance_all_years_UD <- abundance_all_years_UD %>%
  group_by(location, treatment, year) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")







# ---------------------------------------------------------------------
# BACI contrasts function
# ---------------------------------------------------------------------

compute_baci_UD <- function(dat) {
  
  sp_cols <- names(dat)[grepl("^[A-Z]{4}$", names(dat))]
  if (length(sp_cols) == 0) stop("No species columns found (expected 4-letter codes).")
  
  dat <- dat %>%
    mutate(across(all_of(sp_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(sp_cols), ~ replace_na(.x, 0))) %>%
    mutate(year = as.character(year))
  
  abund_2023  <- dat %>% filter(year == "2023")
  abund_after <- dat %>% filter(year %in% c("2024","2025"))
  
  # Always return the right object type
  empty_out <- list(cell = tibble(), baci = tibble())
  if (nrow(abund_2023) == 0 || nrow(abund_after) == 0) return(empty_out)
  
  sum_2023 <- abund_2023 %>%
    group_by(treatment) %>%
    summarise(across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  sum_after <- abund_after %>%
    group_by(treatment) %>%
    summarise(across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  if (nrow(sum_2023) == 0 || nrow(sum_after) == 0) return(empty_out)
  
  to_iNext_list <- function(df) {
    df %>%
      nest(.by = treatment) %>%
      mutate(vec = map(
        data,
        ~ .x %>% select(all_of(sp_cols)) %>% unlist(use.names = FALSE) %>% as.numeric()
      )) %>%
      transmute(treatment, vec) %>%
      deframe() %>%
      lapply(\(x) x[x > 0])
  }
  
  iNext_before <- to_iNext_list(sum_2023)
  iNext_after  <- to_iNext_list(sum_after)
  
  info_before <- DataInfo(iNext_before, datatype = "abundance")
  info_after  <- DataInfo(iNext_after,  datatype = "abundance")
  min_n <- min(c(info_before$n, info_after$n), na.rm = TRUE)
  
  if (!is.finite(min_n) || min_n <= 0) return(empty_out)
  
  est <- function(x, period) {
    estimateD(x, datatype = "abundance", base = "size", level = min_n, conf = 0.95) %>%
      mutate(Period = period)
  }
  
  equal_size <- bind_rows(
    est(iNext_before, "Before"),
    est(iNext_after,  "After")
  )
  
  # ---- Standardized gamma cell means (ALL treatments incl NT) ----
  cell <- equal_size %>%
    transmute(
      Treatment = Assemblage,
      Order.q   = as.integer(Order.q),
      Period,
      qD,
      LCL = qD.LCL,
      UCL = qD.UCL
    )
  
  # If NT missing in either period, we can't compute BACI contrasts
  if (!("NT" %in% cell$Treatment)) {
    return(list(cell = cell, baci = tibble()))
  }
  
  # ---- BACI contrasts: (After-Before)_Treat - (After-Before)_NT ----
  df0 <- cell %>%
    mutate(SE = (UCL - LCL) / (2 * 1.96), Var = SE^2)
  
  wide <- df0 %>%
    select(Treatment, Order.q, Period, qD, Var) %>%
    distinct() %>%
    pivot_wider(names_from = Period, values_from = c(qD, Var), names_sep = "_") %>%
    mutate(
      d_Treat = qD_After - qD_Before,
      v_Treat = Var_After + Var_Before
    )
  
  ctrl <- wide %>%
    filter(Treatment == "NT") %>%
    select(Order.q, d_Ctrl = d_Treat, v_Ctrl = v_Treat)
  
  baci <- wide %>%
    filter(Treatment %in% c("FR","LR","HR")) %>%
    left_join(ctrl, by = "Order.q") %>%
    transmute(
      Order.q   = factor(Order.q, levels = c(0,1,2), labels = c("q = 0","q = 1","q = 2")),
      Treatment = factor(Treatment, levels = c("FR","LR","HR")),
      Estimate  = d_Treat - d_Ctrl,
      SE        = sqrt(v_Treat + v_Ctrl),
      LCL       = Estimate - 1.96 * SE,
      UCL       = Estimate + 1.96 * SE
    )
  
  list(cell = cell, baci = baci)
}






## ============================================================
## 2) Bootstrap (UD): keep RAW replicates for BOTH cell + baci
## ============================================================

set.seed(123)
B <- 1000

loc_tbl   <- abundance_all_years_UD %>% distinct(location, treatment)
loc_sizes <- loc_tbl %>% count(treatment, name = "n_locs")

boot_list_baci <- vector("list", B)
boot_list_cell <- vector("list", B)

failure_baci <- tibble(
  Order.q = factor(rep(c("q = 0", "q = 1", "q = 2"), each = 3),
                   levels = c("q = 0", "q = 1", "q = 2")),
  Treatment = factor(rep(c("FR", "LR", "HR"), times = 3),
                     levels = c("FR", "LR", "HR")),
  Estimate = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_
)

empty_cell <- tibble(
  Treatment = character(),
  Order.q   = integer(),
  Period    = character(),
  qD        = numeric(),
  LCL       = numeric(),
  UCL       = numeric()
)

for (b in seq_len(B)) {
  
  sampled_locs <- loc_tbl %>%
    left_join(loc_sizes, by = "treatment") %>%
    group_by(treatment) %>%
    reframe(location = sample(location, size = n_locs[1], replace = TRUE)) %>%
    ungroup()
  
  boot_dat <- sampled_locs %>%
    left_join(abundance_all_years_UD,
              by = c("location", "treatment"),
              relationship = "many-to-many")
  
  res <- tryCatch(compute_baci_UD(boot_dat), error = function(e) NULL)
  
  baci_b <- if (is.list(res) && all(c("cell", "baci") %in% names(res)) && nrow(res$baci) > 0) res$baci else failure_baci
  cell_b <- if (is.list(res) && all(c("cell", "baci") %in% names(res)) && nrow(res$cell) > 0) res$cell else empty_cell
  
  boot_list_baci[[b]] <- baci_b %>% mutate(boot = b)
  boot_list_cell[[b]] <- cell_b %>% mutate(boot = b)
}

# RAW replicate objects (DO NOT overwrite these with summaries)
boot_baci_raw <- bind_rows(boot_list_baci) %>% filter(!is.na(Estimate))
boot_cell_raw <- bind_rows(boot_list_cell) %>% filter(!is.na(qD))

stopifnot(nrow(boot_baci_raw) > 0, nrow(boot_cell_raw) > 0)










## ============================================================
## 3) Gamma "EMM-like" estimates table (q=0/1): Before/After/Total
## ============================================================

gamma_emm_table_UD <- boot_cell_raw %>%
  filter(Order.q %in% c(0, 1),
         Period %in% c("Before", "After"),
         Treatment %in% c("NT", "FR", "LR", "HR")) %>%
  mutate(
    Metric = factor(Order.q, levels = c(0, 1), labels = c("Species richness", "Shannon diversity")),
    Treatment = factor(Treatment, levels = c("NT", "HR", "LR", "FR")),
    Period    = factor(Period, levels = c("Before", "After"))
  ) %>%
  group_by(Metric, Treatment, Period) %>%
  summarise(
    Mean = mean(qD, na.rm = TRUE),
    LCL  = quantile(qD, 0.025, na.rm = TRUE),
    UCL  = quantile(qD, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CI_half = (UCL - LCL) / 2,
    Cell = paste0(round(Mean, 2), " \u00B1 ", round(CI_half, 2))
  )

wide_cells_gamma_UD <- gamma_emm_table_UD %>%
  select(Metric, Treatment, Period, Cell) %>%
  pivot_wider(names_from = Period, values_from = Cell) %>%
  left_join(
    gamma_emm_table_UD %>%
      group_by(Metric, Treatment) %>%
      summarise(Mean = mean(Mean, na.rm = TRUE),
                CI_half = mean(CI_half, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(`Grand Total` = paste0(round(Mean, 2), " \u00B1 ", round(CI_half, 2))) %>%
      select(Metric, Treatment, `Grand Total`),
    by = c("Metric", "Treatment")
  ) %>%
  arrange(Metric, Treatment)

stopifnot(all(c("Before", "After", "Grand Total") %in% names(wide_cells_gamma_UD)))

## ============================================================
## 4) BACI contrasts table (q=0/1): mean + bootstrap CI + p-values
##    p-value = 2 * min( P(Estimate<=0), P(Estimate>=0) )
## ============================================================

boot_pvals_baci_UD <- boot_baci_raw %>%
  filter(Order.q %in% c("q = 0", "q = 1"),
         Treatment %in% c("FR", "LR", "HR"),
         !is.na(Estimate)) %>%
  group_by(Order.q, Treatment) %>%
  summarise(
    p.value = 2 * min(mean(Estimate <= 0), mean(Estimate >= 0)),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(recode(Order.q, "q = 0" = "Species richness", "q = 1" = "Shannon diversity"),
                    levels = c("Species richness", "Shannon diversity")),
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR"))
  ) %>%
  select(Metric, Treatment, p.value)

## ============================================================
## FIXED: Gamma BACI contrasts WITH NON-ZERO CIs
## Uses RAW bootstrap distribution correctly
## ============================================================

gamma_baci_table_UD <- boot_baci_raw %>%
  filter(
    Order.q %in% c("q = 0","q = 1"),
    Treatment %in% c("FR","LR","HR"),
    !is.na(Estimate)
  ) %>%
  mutate(
    Metric = factor(
      recode(Order.q,
             "q = 0" = "Species richness",
             "q = 1" = "Shannon diversity"),
      levels = c("Species richness","Shannon diversity")
    ),
    Treatment = factor(Treatment, levels = c("FR","LR","HR"))
  ) %>%
  group_by(Metric, Treatment) %>%
  summarise(
    LCL = quantile(Estimate, 0.025, na.rm = TRUE),
    UCL = quantile(Estimate, 0.975, na.rm = TRUE),
    Estimate = mean(Estimate, na.rm = TRUE),
    p.value = 2 * min(
      mean(Estimate <= 0, na.rm = TRUE),
      mean(Estimate >= 0, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    CI_half = (UCL - LCL)/2,
    `Estimate ± half CI` = paste0(round(Estimate, 2), " \u00B1 ", round(CI_half, 2)),
    `p-value`  = signif(p.value, 3),
    `p (Holm)` = signif(p.adjust(p.value, method = "holm"), 3)
  ) %>%
  transmute(
    `Diversity order` = Metric,
    Treatment,
    `Estimate ± half CI`,
    Estimate = round(Estimate, 2),
    `LCL (95% CI)` = round(LCL, 2),
    `UCL (95% CI)` = round(UCL, 2),
    `p-value`,
    `p (Holm)`,
    `CI half-width` = round(CI_half, 2)
  ) %>%
  arrange(`Diversity order`, Treatment)



write.csv(gamma_baci_table_UD, "Output/Tables/GD_UD_BACI_long.csv", row.names = FALSE)

## ============================================================
## 5) FINAL report table (header + estimates + contrasts)
##    (Same structure as your alpha summary table)
## ============================================================

gamma_report_table_UD <- wide_cells_gamma_UD %>%
  mutate(RowType = "estimate", Row = as.character(Treatment)) %>%
  select(Metric, RowType, Row, Before, After, `Grand Total`) %>%
  bind_rows(
    wide_cells_gamma_UD %>% distinct(Metric) %>%
      transmute(Metric, RowType = "header", Row = as.character(Metric),
                Before = "", After = "", `Grand Total` = ""),
    wide_cells_gamma_UD %>% distinct(Metric) %>%
      transmute(Metric, RowType = "subheader", Row = "BACI contrasts (vs NT)",
                Before = "", After = "", `Grand Total` = ""),
    gamma_baci_table_UD %>%
      transmute(
        Metric = `Diversity order`,
        RowType = "contrast",
        Row = paste0(as.character(Treatment), " vs NT"),
        Before = "",
        After = paste0(`Estimate ± half CI`,
                       " (95% CI [", `LCL (95% CI)`, ", ", `UCL (95% CI)`, "]",
                       "; p=", `p-value`, ", Holm=", `p (Holm)`, ")"),
        `Grand Total` = ""
      )
  ) %>%
  arrange(
    Metric,
    factor(RowType, levels = c("header", "estimate", "subheader", "contrast")),
    factor(Row, levels = c("Species richness", "Shannon diversity",
                           "FR", "LR", "HR", "NT",
                           "BACI contrasts (vs NT)",
                           "FR vs NT", "LR vs NT", "HR vs NT"))
  ) %>%
  select(Row, Before, After, `Grand Total`) %>%
  bind_rows(tibble(Row = "Grand Total", Before = "", After = "", `Grand Total` = ""))

# Print
gamma_report_table_UD

# Save
write.csv(gamma_report_table_UD, "Output/Tables/GD_UD_EMMplusBACI_report.csv", row.names = FALSE)








################### 9b beta div LA

# ---
# title: "Beta diversity LIMITED AMPLITUDE"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script create extrapolation and rarefaction curves for Hill numbers of the order q = 0, 1, and 2
# ---

# Load package
library(betapart)
library(tidyverse)

# Load data
dat1_div_LA <- read.csv("Output/Tabular Data//max_count_all_years_LA.csv")

# (Corrected: removed duplicate UDPC load)
# dat1_div_LA <- read.csv("Output/Tabular Data//max_count_all_years_UDPC.csv")

# Create treatment column without remove original location column
dat2_div_LA <- dat1_div_LA %>%
  mutate(site = location)

dat3_div_LA <- dat2_div_LA %>%
  separate_wider_delim(site, delim = "-", names = c("block", "treatment", "plot"))

# Identify species columns
species_cols <- setdiff(names(dat3_div_LA), c("location", "year", "block", "treatment", "plot"))

# Remove unneeded columns
dat4_div_LA <- dat3_div_LA %>%
  select(location, treatment, year, all_of(species_cols))

# Convert from counts to presence/absence
dat5_div_LA <- dat4_div_LA %>%
  mutate(across(
    all_of(species_cols),
    ~ ifelse(is.na(.) | . == 0, 0, 1)
  ))

### Subset df containing all years and treatment into individual df,
### one per treatment-year combination

#FR
FR_2023 <- dat5_div_LA %>% filter(treatment == "FR" & year == "2023")
FR_2024 <- dat5_div_LA %>% filter(treatment == "FR" & year == "2024")
FR_2025 <- dat5_div_LA %>% filter(treatment == "FR" & year == "2025")

#LR
LR_2023 <- dat5_div_LA %>% filter(treatment == "LR" & year == "2023")
LR_2024 <- dat5_div_LA %>% filter(treatment == "LR" & year == "2024")
LR_2025 <- dat5_div_LA %>% filter(treatment == "LR" & year == "2025")

#HR
HR_2023 <- dat5_div_LA %>% filter(treatment == "HR" & year == "2023")
HR_2024 <- dat5_div_LA %>% filter(treatment == "HR" & year == "2024")
HR_2025 <- dat5_div_LA %>% filter(treatment == "HR" & year == "2025")

#NT
NT_2023 <- dat5_div_LA %>% filter(treatment == "NT" & year == "2023")
NT_2024 <- dat5_div_LA %>% filter(treatment == "NT" & year == "2024")
NT_2025 <- dat5_div_LA %>% filter(treatment == "NT" & year == "2025")

### Assign location to row names
rownames(FR_2023) <- FR_2023$location
rownames(FR_2024) <- FR_2024$location
rownames(FR_2025) <- FR_2025$location

rownames(LR_2023) <- LR_2023$location
rownames(LR_2024) <- LR_2024$location
rownames(LR_2025) <- LR_2025$location

rownames(HR_2023) <- HR_2023$location
rownames(HR_2024) <- HR_2024$location
rownames(HR_2025) <- HR_2025$location

rownames(NT_2023) <- NT_2023$location
rownames(NT_2024) <- NT_2024$location
rownames(NT_2025) <- NT_2025$location

## Remove unneeded columns and convert to matrix
#FR
FR_2023_matrix <- FR_2023 %>% select(all_of(species_cols)) %>% as.matrix()
FR_2024_matrix <- FR_2024 %>% select(all_of(species_cols)) %>% as.matrix()
FR_2025_matrix <- FR_2025 %>% select(all_of(species_cols)) %>% as.matrix()

#LR
LR_2023_matrix <- LR_2023 %>% select(all_of(species_cols)) %>% as.matrix()
LR_2024_matrix <- LR_2024 %>% select(all_of(species_cols)) %>% as.matrix()
LR_2025_matrix <- LR_2025 %>% select(all_of(species_cols)) %>% as.matrix()

#HR
HR_2023_matrix <- HR_2023 %>% select(all_of(species_cols)) %>% as.matrix()
HR_2024_matrix <- HR_2024 %>% select(all_of(species_cols)) %>% as.matrix()
HR_2025_matrix <- HR_2025 %>% select(all_of(species_cols)) %>% as.matrix()

#NT
NT_2023_matrix <- NT_2023 %>% select(all_of(species_cols)) %>% as.matrix()
NT_2024_matrix <- NT_2024 %>% select(all_of(species_cols)) %>% as.matrix()
NT_2025_matrix <- NT_2025 %>% select(all_of(species_cols)) %>% as.matrix()

## Convert from matrices to betapart objects
FR_2023_BPC <- betapart.core(FR_2023_matrix)
FR_2024_BPC <- betapart.core(FR_2024_matrix)
FR_2025_BPC <- betapart.core(FR_2025_matrix)

LR_2023_BPC <- betapart.core(LR_2023_matrix)
LR_2024_BPC <- betapart.core(LR_2024_matrix)
LR_2025_BPC <- betapart.core(LR_2025_matrix)

HR_2023_BPC <- betapart.core(HR_2023_matrix)
HR_2024_BPC <- betapart.core(HR_2024_matrix)
HR_2025_BPC <- betapart.core(HR_2025_matrix)

NT_2023_BPC <- betapart.core(NT_2023_matrix)
NT_2024_BPC <- betapart.core(NT_2024_matrix)
NT_2025_BPC <- betapart.core(NT_2025_matrix)

### Calculate beta diversity
#2023 vs 2024
NT_2023_2024 <- beta.temp(NT_2023_matrix, NT_2024_matrix, index.family="sorensen") %>% mutate(treatment = "NT_2023_2024")
FR_2023_2024 <- beta.temp(FR_2023_matrix, FR_2024_matrix, index.family="sorensen") %>% mutate(treatment = "FR_2023_2024")
LR_2023_2024 <- beta.temp(LR_2023_matrix, LR_2024_matrix, index.family="sorensen") %>% mutate(treatment = "LR_2023_2024")
HR_2023_2024 <- beta.temp(HR_2023_matrix, HR_2024_matrix, index.family="sorensen") %>% mutate(treatment = "HR_2023_2024")

# Bind
all_treatments_2023_2024_LA <- bind_rows(NT_2023_2024, FR_2023_2024, HR_2023_2024, LR_2023_2024)

########## Stacked bar plot 2024
beta_summary_2024_LA <- all_treatments_2023_2024_LA %>%
  group_by(treatment) %>%
  summarise(
    mean_turnover   = mean(beta.sim, na.rm = TRUE),
    mean_nestedness = mean(beta.sne, na.rm = TRUE),
    mean_total      = mean(beta.sor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = factor(treatment, levels = sort(unique(treatment))))

beta_long_2024_LA <- beta_summary_2024_LA %>%
  select(treatment, mean_turnover, mean_nestedness) %>%
  pivot_longer(cols = c(mean_turnover, mean_nestedness),
               names_to = "Component",
               values_to = "Mean_beta") %>%
  mutate(Component = recode(
    Component,
    mean_turnover   = "Turnover (β_sim)",
    mean_nestedness = "Nestedness (β_sne)"
  ))

p_beta_stack_2024_LA <- ggplot(beta_long_2024_LA,
                               aes(x = treatment, y = Mean_beta, fill = Component)) +
  geom_col(color = "black") +
  labs(
    x = "Treatment–year comparison",
    y = "Mean β-diversity (Sørensen family)",
    fill = "Component",
    title = "β by turnover and nestedness - Limited Amplitude"
  ) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

print(p_beta_stack_2024_LA)

### Calculate beta diversity
#2023 vs 2025
NT_2023_2025 <- beta.temp(NT_2023_matrix, NT_2025_matrix, index.family="sorensen") %>% mutate(treatment = "NT_2023_2025")
FR_2023_2025 <- beta.temp(FR_2023_matrix, FR_2025_matrix, index.family="sorensen") %>% mutate(treatment = "FR_2023_2025")
LR_2023_2025 <- beta.temp(LR_2023_matrix, LR_2025_matrix, index.family="sorensen") %>% mutate(treatment = "LR_2023_2025")
HR_2023_2025 <- beta.temp(HR_2023_matrix, HR_2025_matrix, index.family="sorensen") %>% mutate(treatment = "HR_2023_2025")

# Bind
all_treatments_2023_2025_LA <- bind_rows(NT_2023_2025, FR_2023_2025, HR_2023_2025, LR_2023_2025)

#### Stacked bar plot 2025
beta_summary_2025_LA <- all_treatments_2023_2025_LA %>%
  group_by(treatment) %>%
  summarise(
    mean_turnover   = mean(beta.sim, na.rm = TRUE),
    mean_nestedness = mean(beta.sne, na.rm = TRUE),
    mean_total      = mean(beta.sor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = factor(treatment, levels = sort(unique(treatment))))

beta_long_2025_LA <- beta_summary_2025_LA %>%
  select(treatment, mean_turnover, mean_nestedness) %>%
  pivot_longer(cols = c(mean_turnover, mean_nestedness),
               names_to = "Component",
               values_to = "Mean_beta") %>%
  mutate(Component = recode(
    Component,
    mean_turnover   = "Turnover (β_sim)",
    mean_nestedness = "Nestedness (β_sne)"
  ))

p_beta_stack_2025_LA <- ggplot(beta_long_2025_LA,
                               aes(x = treatment, y = Mean_beta, fill = Component)) +
  geom_col(color = "black") +
  labs(
    x = "Treatment–year comparison",
    y = "Mean β-diversity (Sørensen family)",
    fill = "Component",
    title = "β by turnover and nestedness - Limited Amplitude"
  ) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

print(p_beta_stack_2025_LA)







## ================================================================
## Community trajectories in NMDS space
## ================================================================

library(vegan)
library(tidyverse)
library(grid)

comm_mat_LA <- dat5_div_LA %>%
  select(all_of(species_cols)) %>%
  as.matrix()

meta_LA <- dat5_div_LA %>%
  select(location, treatment, year) %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    year      = factor(year, levels = c("2023", "2024", "2025"))
  )

set.seed(123)
nmds_LA <- metaMDS(comm_mat_LA, distance = "bray", k = 2, trymax = 100)

site_scores_LA <- as.data.frame(scores(nmds_LA, display = "sites")) %>%
  mutate(
    location  = meta_LA$location,
    treatment = meta_LA$treatment,
    year      = meta_LA$year
  )

centroids_LA <- site_scores_LA %>%
  group_by(treatment, year) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  ) %>%
  arrange(treatment, year)

traj_length_LA <- centroids_LA %>%
  arrange(treatment, year) %>%
  group_by(treatment) %>%
  summarise(
    traj_length = sum(
      sqrt(diff(NMDS1)^2 + diff(NMDS2)^2)
    ),
    .groups = "drop"
  )
print(traj_length_LA)

gg_traj_LA <- ggplot(centroids_LA,
                     aes(x = NMDS1, y = NMDS2,
                         colour = treatment, group = treatment)) +
  geom_path(arrow = arrow(type = "closed", length = unit(0.2, "cm")),
            linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = year),
            nudge_y = 0.03, size = 3) +
  scale_colour_brewer(palette = "Dark2",
                      name = "Treatment") +
  coord_equal() +
  labs(
    title = "Community trajectories in NMDS space (Limited Amplitude)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(gg_traj_LA)




































############## GOOD ##################


################## POOLING 2024/2025 for betapart


# ---
# title: "Beta diversity unlimited distance (Averaged 2024 & 2025 Contrast)"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script calculates beta diversity turnover and nestedness components,
#              averaging the 2023-2024 and 2023-2025 contrasts for a robust Before vs After estimate."
# ---

rm(list = ls())  # Removes all objects from the environment

# Load package
library(betapart)
library(tidyverse)
library(vegan) # For NMDS later

# Load data
dat1_div_LA<- read.csv("Output/Tabular Data//max_count_all_years_LA.csv")

# Create treatment column without remove original location column
dat2_div_LA <- dat1_div_LA %>%
  mutate(site = location)

dat3_div_LA <- dat2_div_LA %>%
  separate_wider_delim(site, delim = "-", names = c("block", "treatment", "plot"))

# Identify species columns
species_cols <- setdiff(names(dat3_div_LA), c("location", "year", "block", "treatment", "plot"))

# Remove unneeded columns
dat4_div_LA <- dat3_div_LA %>%
  select(location, treatment, year, all_of(species_cols))

# Convert from counts to presence/absence
dat5_div_LA <- dat4_div_LA %>%
  mutate(across(
    all_of(species_cols),
    ~ ifelse(is.na(.) | . == 0, 0, 1)
  ))


### Subset df containing all years and treatment into individual df,
### one per treatment-year combination 

#FR
FR_2023 <- dat5_div_LA %>%
  filter(treatment == "FR" & year == "2023")

# ... (ALL other subsets for FR_2024, FR_2025, LR_2023 to NT_2025 remain here) ...
FR_2024 <- dat5_div_LA %>%
  filter(treatment == "FR" & year == "2024")

FR_2025 <- dat5_div_LA %>%
  filter(treatment == "FR" & year == "2025")

#LR
LR_2023 <- dat5_div_LA %>%
  filter(treatment == "LR" & year == "2023")

LR_2024 <- dat5_div_LA %>%
  filter(treatment == "LR" & year == "2024")

LR_2025 <- dat5_div_LA %>%
  filter(treatment == "LR" & year == "2025")


#HR
HR_2023 <- dat5_div_LA %>%
  filter(treatment == "HR" & year == "2023")

HR_2024 <- dat5_div_LA %>%
  filter(treatment == "HR" & year == "2024")

HR_2025 <- dat5_div_LA %>%
  filter(treatment == "HR" & year == "2025")


#NT
NT_2023 <- dat5_div_LA %>%
  filter(treatment == "NT" & year == "2023")

NT_2024 <- dat5_div_LA %>%
  filter(treatment == "NT" & year == "2024")

NT_2025 <- dat5_div_LA %>%
  filter(treatment == "NT" & year == "2025")

### Assign location to row names 
rownames(FR_2023) <- FR_2023$location
rownames(FR_2024) <- FR_2024$location
rownames(FR_2025) <- FR_2025$location

rownames(LR_2023) <- LR_2023$location
rownames(LR_2024) <- LR_2024$location
rownames(LR_2025) <- LR_2025$location

rownames(HR_2023) <- HR_2023$location
rownames(HR_2024) <- HR_2024$location
rownames(HR_2025) <- HR_2025$location

rownames(NT_2023) <- NT_2023$location
rownames(NT_2024) <- NT_2024$location
rownames(NT_2025) <- NT_2025$location


## Remove unneeded columns and convert to matrix 
#FR
FR_2023_matrix <- FR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

FR_2024_matrix <- FR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

FR_2025_matrix <- FR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#LR
LR_2023_matrix <- LR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

LR_2024_matrix <- LR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

LR_2025_matrix <- LR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#HR
HR_2023_matrix <- HR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

HR_2024_matrix <- HR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

HR_2025_matrix <- HR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#NT
NT_2023_matrix <- NT_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

NT_2024_matrix <- NT_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

NT_2025_matrix <- NT_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()


## Calculation of Beta Diversity Contrasts

### Contrast 1: 2023 vs 2024
NT_2023_2024 <- beta.temp(NT_2023_matrix, NT_2024_matrix, index.family="jaccard")
FR_2023_2024 <- beta.temp(FR_2023_matrix, FR_2024_matrix, index.family="jaccard")
LR_2023_2024 <- beta.temp(LR_2023_matrix, LR_2024_matrix, index.family="jaccard")
HR_2023_2024 <- beta.temp(HR_2023_matrix, HR_2024_matrix, index.family="jaccard")

# Add treatment and contrast group identifier
all_treatments_2023_2024_LA <- bind_rows(
  NT_2023_2024 %>% mutate(treatment = "NT"),
  FR_2023_2024 %>% mutate(treatment = "FR"),
  LR_2023_2024 %>% mutate(treatment = "LR"),
  HR_2023_2024 %>% mutate(treatment = "HR")
) %>% mutate(Contrast_Group = "2023_vs_2024")


### Contrast 2: 2023 vs 2025
NT_2023_2025 <- beta.temp(NT_2023_matrix, NT_2025_matrix, index.family="jaccard")
FR_2023_2025 <- beta.temp(FR_2023_matrix, FR_2025_matrix, index.family="jaccard")
LR_2023_2025 <- beta.temp(LR_2023_matrix, LR_2025_matrix, index.family="jaccard")
HR_2023_2025 <- beta.temp(HR_2023_matrix, HR_2025_matrix, index.family="jaccard")

# Add treatment and contrast group identifier
all_treatments_2023_2025_LA <- bind_rows(
  NT_2023_2025 %>% mutate(treatment = "NT"),
  FR_2023_2025 %>% mutate(treatment = "FR"),
  LR_2023_2025 %>% mutate(treatment = "LR"),
  HR_2023_2025 %>% mutate(treatment = "HR")
) %>% mutate(Contrast_Group = "2023_vs_2025")



## ----------------------------------------------------------------
## Averaging the Contrasts for a Final Before vs. After Estimate
## ----------------------------------------------------------------

# 1. Combine the two contrast tables (No change)
final_beta_data_LA <- bind_rows(all_treatments_2023_2024_LA, all_treatments_2023_2025_LA)

# --- NEW: Calculate 95% CIs based on the Standard Error of the Mean (SEM) ---
# Since N = 12 (6 sites * 2 contrasts) per treatment, df = 11.
# The 95% t-score for df=11 is qt(0.975, 11) ≈ 2.201. We'll use this.
t_score <- qt(0.975, 11) # Degrees of Freedom = N - 1 = 12 - 1 = 11

# 2. Summarise mean components and calculate CIs
beta_summary_FINAL_LA <- final_beta_data_LA %>%
  group_by(treatment) %>%
  summarise(
    # Turnover Component (beta.jtu)
    mean_turnover = mean(beta.jtu, na.rm = TRUE),
    sd_turnover = sd(beta.jtu, na.rm = TRUE),
    # Nestedness Component (beta.jne)
    mean_nestedness = mean(beta.jne, na.rm = TRUE),
    sd_nestedness = sd(beta.jne, na.rm = TRUE),
    # Total Beta (beta.jac)
    mean_total = mean(beta.jac, na.rm = TRUE),
    sd_total = sd(beta.jac, na.rm = TRUE),
    n_contrasts = n(), # Should be 12 for all treatments
    .groups = "drop"
  ) %>%
  # Calculate 95% CI bounds for the TOTAL beta diversity (mean_total)
  mutate(
    SE_total = sd_total / sqrt(n_contrasts),
    LCL_total = mean_total - t_score * SE_total,
    UCL_total = mean_total + t_score * SE_total
  )


# 3. Prepare for plotting: Pivot to long format (No change to logic)
beta_long_FINAL_LA <- beta_summary_FINAL_LA %>%
  select(treatment, mean_turnover, mean_nestedness, mean_total, LCL_total, UCL_total) %>%
  pivot_longer(
    cols = c(mean_turnover, mean_nestedness),
    names_to = "Component",
    values_to = "Mean_beta"
  ) %>%
  mutate(
    Component = recode(
      Component,
      mean_turnover = "Species replacement",
      mean_nestedness = "Species loss"
    ),
    # Ensure treatment order is logical for the plot
    treatment = factor(treatment, levels = c("NT", "FR", "HR", "LR"))
  )

# 4. Final Stacked Bar Plot with CIs on Total Mean
p_beta_stack_FINAL_LA <- ggplot(beta_long_FINAL_LA,
                                aes(x = treatment,
                                    y = Mean_beta,
                                    fill = Component)) +
  geom_col(color = "black") +
  # --- FIX APPLIED HERE: Added x = treatment to the aes() mapping ---
  geom_errorbar(aes(x = treatment, # <-- ADDED X-MAPPING
                    ymin = LCL_total,
                    ymax = UCL_total, 
                    y = mean_total),
                width = 0.2, 
                linewidth = 0.8,
                color = "black",
                inherit.aes = FALSE, 
                data = beta_summary_FINAL_LA
  ) +
  labs(
    x = "Treatment",
    y = "Mean β-diversity (Jaccard family)",
    fill = "Component",
    title = "β-Diversity: Mean Change (2023 vs. 2024/2025)(Limited Amplitude"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_beta_stack_FINAL_LA)








## ================================================================
## Community trajectories in NMDS space 
## ================================================================

# ... (NMDS section remains exactly as it was, using all three years for trajectory plotting) ...

## 1. Build community matrix and metadata ---------------------------

# Site × species matrix (can be presence/absence or abundance)
comm_mat_LA <- dat5_div_LA %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# Metadata for each row in the matrix
meta_LA <- dat5_div_LA %>%
  select(location, treatment, year) %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    year      = factor(year, levels = c("2023", "2024", "2025"))
  )

## 2. NMDS ordination (Bray–Curtis) --------------------------------

set.seed(123)
# Note: Using 'bray' on P/A data is equivalent to Jaccard-like measures.
nmds_LA <- metaMDS(comm_mat_LA, distance = "bray", k = 2, trymax = 1000)

# Site scores
site_scores_LA <- as.data.frame(scores(nmds_LA, display = "sites")) %>%
  mutate(
    location  = meta_LA$location,
    treatment = meta_LA$treatment,
    year      = meta_LA$year
  )

## 3. Treatment × year centroids in ordination space ---------------

centroids_LA <- site_scores_LA %>%
  group_by(treatment, year) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  ) %>%
  arrange(treatment, year)

# Optional: trajectory length (amount of compositional change)
traj_length_LA <- centroids_LA %>%
  arrange(treatment, year) %>%
  group_by(treatment) %>%
  summarise(
    traj_length = sum(
      sqrt(diff(NMDS1)^2 + diff(NMDS2)^2)
    ),
    .groups = "drop"
  )
print(traj_length_LA)

## 4. Trajectory plot (arrows between years per treatment) ---------

gg_traj_LA <- ggplot(centroids_LA,
                     aes(x = NMDS1, y = NMDS2,
                         colour = treatment, group = treatment)) +
  # arrows showing direction of change across years within treatment
  geom_path(arrow = arrow(type = "closed", length = unit(0.2, "cm")),
            linewidth = 0.8) +
  # centroid points
  geom_point(size = 3) +
  # year labels next to points
  geom_text(aes(label = year),
            nudge_y = 0.03, size = 3) +
  scale_colour_brewer(palette = "Dark2",
                      name = "Treatment") +
  coord_equal() +
  labs(
    title = "Community trajectories in NMDS space (Limited Amplitude)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(gg_traj_LA)















########################### 9a beta div UD

# ---
# title: "Beta diversity unlimited distance"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script create extrapolation and rarefaction curves for Hill numbers of the order q = 0, 1, and 2
# ---


# Load package
library(betapart)
library(tidyverse)

# Load data
dat1_div_UD <- read.csv("Output/Tabular Data//max_count_all_years_UDPC.csv")

# Create treatment column without remove original location column
dat2_div_UD <- dat1_div_UD %>%
  mutate(site = location)

dat3_div_UD <- dat2_div_UD %>%
  separate_wider_delim(site, delim = "-", names = c("block", "treatment", "plot"))

# Identify species columns
species_cols <- setdiff(names(dat3_div_UD), c("location", "year", "block", "treatment", "plot"))

# Remove unneeded columns
dat4_div_UD <- dat3_div_UD %>%
  select(location, treatment, year, all_of(species_cols))

# Convert from counts to presence/absence
dat5_div_UD <- dat4_div_UD %>%
  mutate(across(
    all_of(species_cols),
    ~ ifelse(is.na(.) | . == 0, 0, 1)
  ))


### Subset df containing all years and treatment into individual df,
### one per treatment-year combination

#FR
FR_2023 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2023")

FR_2024 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2024")

FR_2025 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2025")

#LR
LR_2023 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2023")

LR_2024 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2024")

LR_2025 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2025")


#HR
HR_2023 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2023")

HR_2024 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2024")

HR_2025 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2025")


#NT
NT_2023 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2023")

NT_2024 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2024")

NT_2025 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2025")




### Assign location to row names
rownames(FR_2023) <- FR_2023$location
rownames(FR_2024) <- FR_2024$location
rownames(FR_2025) <- FR_2025$location

rownames(LR_2023) <- LR_2023$location
rownames(LR_2024) <- LR_2024$location
rownames(LR_2025) <- LR_2025$location

rownames(HR_2023) <- HR_2023$location
rownames(HR_2024) <- HR_2024$location
rownames(HR_2025) <- HR_2025$location

rownames(NT_2023) <- NT_2023$location
rownames(NT_2024) <- NT_2024$location
rownames(NT_2025) <- NT_2025$location



## Remove unneeded columns and convert to matrix
#FR
FR_2023_matrix <- FR_2023 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

FR_2024_matrix <- FR_2024 %>%
  select(all_of(species_cols)) %>% 
  as.matrix()

FR_2025_matrix <- FR_2025 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

#LR
LR_2023_matrix <- LR_2023 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

LR_2024_matrix <- LR_2024 %>%
  select(all_of(species_cols)) %>% 
  as.matrix()

LR_2025_matrix <- LR_2025 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

#HR
HR_2023_matrix <- HR_2023 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

HR_2024_matrix <- HR_2024 %>%
  select(all_of(species_cols)) %>% 
  as.matrix()

HR_2025_matrix <- HR_2025 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

#NT
NT_2023_matrix <- NT_2023 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()

NT_2024_matrix <- NT_2024 %>%
  select(all_of(species_cols)) %>% 
  as.matrix()

NT_2025_matrix <- NT_2025 %>%
  select(all_of(species_cols)) %>%  
  as.matrix()


## Convert from matrices to betapart objects
FR_2023_BPC <- betapart.core(FR_2023_matrix)
FR_2024_BPC <- betapart.core(FR_2024_matrix)
FR_2025_BPC <- betapart.core(FR_2025_matrix)

LR_2023_BPC <- betapart.core(LR_2023_matrix)
LR_2024_BPC <- betapart.core(LR_2024_matrix)
LR_2025_BPC <- betapart.core(LR_2025_matrix)

HR_2023_BPC <- betapart.core(HR_2023_matrix)
HR_2024_BPC <- betapart.core(HR_2024_matrix)
HR_2025_BPC <- betapart.core(HR_2025_matrix)

NT_2023_BPC <- betapart.core(NT_2023_matrix)
NT_2024_BPC <- betapart.core(NT_2024_matrix)
NT_2025_BPC <- betapart.core(NT_2025_matrix)










### Calculate beta diversity
#2023 vs 2024
NT_2023_2024 <- beta.temp(NT_2023_matrix, NT_2024_matrix, index.family="sorensen")
FR_2023_2024 <- beta.temp(FR_2023_matrix, FR_2024_matrix, index.family="sorensen")
LR_2023_2024 <- beta.temp(LR_2023_matrix, LR_2024_matrix, index.family="sorensen")
HR_2023_2024 <- beta.temp(HR_2023_matrix, HR_2024_matrix, index.family="sorensen")

# Add treatment back to beta.part output
NT_2023_2024 <- NT_2023_2024 %>%
  mutate(treatment = "NT_2023_2024")

FR_2023_2024 <- FR_2023_2024 %>%
  mutate(treatment = "FR_2023_2024")

LR_2023_2024 <- LR_2023_2024 %>%
  mutate(treatment = "LR_2023_2024")

HR_2023_2024 <- HR_2023_2024 %>%
  mutate(treatment = "HR_2023_2024")

# Bind
all_treatments_2023_2024_UD <- bind_rows(NT_2023_2024, FR_2023_2024, HR_2023_2024, LR_2023_2024)







########## Stacked bar plot 2024


# 1. Summarise mean components by treatment
beta_summary_2024_UD <- all_treatments_2023_2024_UD %>%
  group_by(treatment) %>%
  summarise(
    mean_turnover   = mean(beta.sim, na.rm = TRUE),  # β_sim
    mean_nestedness = mean(beta.sne, na.rm = TRUE),  # β_sne
    mean_total      = mean(beta.sor, na.rm = TRUE),  # β_sor (≈ turnover + nestedness)
    .groups = "drop"
  )

# (Optional) clean up treatment order if you like
# e.g., put NT first, then FR, HR, LR, etc.
beta_summary_2024_UD <- beta_summary_2024_UD %>%
  mutate(
    treatment = factor(
      treatment,
      levels = sort(unique(treatment))  # or specify manually
      # levels = c("NT_2023_2024","FR_2023_2024","HR_2023_2024","LR_2023_2024")
    )
  )

# 2. Pivot to long format for stacked bars (turnover + nestedness)
beta_long_2024_UD <- beta_summary_2024_UD %>%
  select(treatment, mean_turnover, mean_nestedness) %>%
  pivot_longer(
    cols = c(mean_turnover, mean_nestedness),
    names_to = "Component",
    values_to = "Mean_beta"
  ) %>%
  mutate(
    Component = recode(
      Component,
      mean_turnover   = "Turnover (β_sim)",
      mean_nestedness = "Nestedness (β_sne)"
    )
  )

# 3. Stacked bar plot
p_beta_stack_2024_UD <- ggplot(beta_long_2024_UD,
                               aes(x = treatment,
                                   y = Mean_beta,
                                   fill = Component)) +
  geom_col(color = "black") +
  labs(
    x = "Treatment–year comparison",
    y = "Mean β-diversity (Sørensen family)",
    fill = "Component",
    title = "β by turnover and nestedness - Unlimited Distance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_beta_stack_2024_UD)
















### Calculate beta diversity
#2023 vs 2025
NT_2023_2025 <- beta.temp(NT_2023_matrix, NT_2025_matrix, index.family="sorensen")
FR_2023_2025 <- beta.temp(FR_2023_matrix, FR_2025_matrix, index.family="sorensen")
LR_2023_2025 <- beta.temp(LR_2023_matrix, LR_2025_matrix, index.family="sorensen")
HR_2023_2025 <- beta.temp(HR_2023_matrix, HR_2025_matrix, index.family="sorensen")

# Add treatment back to beta.part output
NT_2023_2025 <- NT_2023_2025 %>%
  mutate(treatment = "NT_2023_2025")

FR_2023_2025 <- FR_2023_2025 %>%
  mutate(treatment = "FR_2023_2025")

LR_2023_2025 <- LR_2023_2025 %>%
  mutate(treatment = "LR_2023_2025")

HR_2023_2025 <- HR_2023_2025 %>%
  mutate(treatment = "HR_2023_2025")

# Bind
all_treatments_2023_2025_UD <- bind_rows(NT_2023_2025, FR_2023_2025, HR_2023_2025, LR_2023_2025)







#### Stacked bar plot

## 2025

# Summarise mean components by treatment
beta_summary_2025_UD <- all_treatments_2023_2025_UD %>%
  group_by(treatment) %>%
  summarise(
    mean_turnover   = mean(beta.sim, na.rm = TRUE),  # β_sim
    mean_nestedness = mean(beta.sne, na.rm = TRUE),  # β_sne
    mean_total      = mean(beta.sor, na.rm = TRUE),  # β_sor (≈ turnover + nestedness)
    .groups = "drop"
  )

# (Optional) clean up treatment order if you like
# e.g., put NT first, then FR, HR, LR, etc.
beta_summary_2025_UD <- beta_summary_2025_UD %>%
  mutate(
    treatment = factor(
      treatment,
      levels = sort(unique(treatment))  # or specify manually
      # levels = c("NT_2023_2025","FR_2023_2025","HR_2023_2025","LR_2023_2025")
    )
  )

# Pivot to long format for stacked bars (turnover + nestedness)
beta_long_2025_UD <- beta_summary_2025_UD %>%
  select(treatment, mean_turnover, mean_nestedness) %>%
  pivot_longer(
    cols = c(mean_turnover, mean_nestedness),
    names_to = "Component",
    values_to = "Mean_beta"
  ) %>%
  mutate(
    Component = recode(
      Component,
      mean_turnover   = "Turnover (β_sim)",
      mean_nestedness = "Nestedness (β_sne)"
    )
  )

# Stacked bar plot
p_beta_stack_2025_UD <- ggplot(beta_long,
                               aes(x = treatment,
                                   y = Mean_beta,
                                   fill = Component)) +
  geom_col(color = "black") +
  labs(
    x = "Treatment–year comparison",
    y = "Mean β-diversity (Sørensen family)",
    fill = "Component",
    title = "β by turnover and nestedness - Unlimited Distance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_beta_stack_2025_UD)




















## ================================================================
## Community trajectories in NMDS space (unlimited distance example)
## Requires:
##   dat5_div_UD : data frame with cols location, treatment, year, species_cols
##   species_cols: character vector of species column names
## ================================================================

library(vegan)
library(tidyverse)
library(grid)   # for arrow()

## 1. Build community matrix and metadata ---------------------------

# Site × species matrix (can be presence/absence or abundance)
comm_mat_UD <- dat5_div_UD %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# Metadata for each row in the matrix
meta_UD <- dat5_div_UD %>%
  select(location, treatment, year) %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    year      = factor(year, levels = c("2023", "2024", "2025"))
  )

## 2. NMDS ordination (Bray–Curtis) --------------------------------

set.seed(123)
nmds_UD <- metaMDS(comm_mat_UD, distance = "bray", k = 2, trymax = 1000)

# Site scores
site_scores_UD <- as.data.frame(scores(nmds_UD, display = "sites")) %>%
  mutate(
    location  = meta_UD$location,
    treatment = meta_UD$treatment,
    year      = meta_UD$year
  )

## 3. Treatment × year centroids in ordination space ---------------

centroids_UD <- site_scores_UD %>%
  group_by(treatment, year) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  ) %>%
  arrange(treatment, year)

# Optional: trajectory length (amount of compositional change)
traj_length_UD <- centroids_UD %>%
  arrange(treatment, year) %>%
  group_by(treatment) %>%
  summarise(
    traj_length = sum(
      sqrt(diff(NMDS1)^2 + diff(NMDS2)^2)
    ),
    .groups = "drop"
  )
print(traj_length_UD)

## 4. Trajectory plot (arrows between years per treatment) ---------

gg_traj_UD <- ggplot(centroids_UD,
                     aes(x = NMDS1, y = NMDS2,
                         colour = treatment, group = treatment)) +
  # arrows showing direction of change across years within treatment
  geom_path(arrow = arrow(type = "closed", length = unit(0.2, "cm")),
            linewidth = 0.8) +
  # centroid points
  geom_point(size = 3) +
  # year labels next to points
  geom_text(aes(label = year),
            nudge_y = 0.03, size = 3) +
  scale_colour_brewer(palette = "Dark2",
                      name = "Treatment") +
  coord_equal() +
  labs(
    title = "Community trajectories in NMDS space (unlimited distance)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(gg_traj_UD)
















############## GOOD ##################


################## POOLING 2024/2025 for betapart


# ---
# title: "Beta diversity unlimited distance (Averaged 2024 & 2025 Contrast)"
# author: "Leonard Patterson"
# created: "2025-07-04"
# description: "This script calculates beta diversity turnover and nestedness components,
#              averaging the 2023-2024 and 2023-2025 contrasts for a robust Before vs After estimate."
# ---

rm(list = ls())  # Removes all objects from the environment

# Load package
library(betapart)
library(tidyverse)
library(vegan) # For NMDS later

# Load data
dat1_div_UD <- read.csv("Output/Tabular Data//max_count_all_years_UDPC.csv")

# Create treatment column without remove original location column
dat2_div_UD <- dat1_div_UD %>%
  mutate(site = location)

dat3_div_UD <- dat2_div_UD %>%
  separate_wider_delim(site, delim = "-", names = c("block", "treatment", "plot"))

# Identify species columns
species_cols <- setdiff(names(dat3_div_UD), c("location", "year", "block", "treatment", "plot"))

# Remove unneeded columns
dat4_div_UD <- dat3_div_UD %>%
  select(location, treatment, year, all_of(species_cols))

# Convert from counts to presence/absence
dat5_div_UD <- dat4_div_UD %>%
  mutate(across(
    all_of(species_cols),
    ~ ifelse(is.na(.) | . == 0, 0, 1)
  ))


### Subset df containing all years and treatment into individual df,
### one per treatment-year combination 

#FR
FR_2023 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2023")

# ... (ALL other subsets for FR_2024, FR_2025, LR_2023 to NT_2025 remain here) ...
FR_2024 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2024")

FR_2025 <- dat5_div_UD %>%
  filter(treatment == "FR" & year == "2025")

#LR
LR_2023 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2023")

LR_2024 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2024")

LR_2025 <- dat5_div_UD %>%
  filter(treatment == "LR" & year == "2025")


#HR
HR_2023 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2023")

HR_2024 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2024")

HR_2025 <- dat5_div_UD %>%
  filter(treatment == "HR" & year == "2025")


#NT
NT_2023 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2023")

NT_2024 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2024")

NT_2025 <- dat5_div_UD %>%
  filter(treatment == "NT" & year == "2025")

### Assign location to row names 
rownames(FR_2023) <- FR_2023$location
rownames(FR_2024) <- FR_2024$location
rownames(FR_2025) <- FR_2025$location

rownames(LR_2023) <- LR_2023$location
rownames(LR_2024) <- LR_2024$location
rownames(LR_2025) <- LR_2025$location

rownames(HR_2023) <- HR_2023$location
rownames(HR_2024) <- HR_2024$location
rownames(HR_2025) <- HR_2025$location

rownames(NT_2023) <- NT_2023$location
rownames(NT_2024) <- NT_2024$location
rownames(NT_2025) <- NT_2025$location


## Remove unneeded columns and convert to matrix 
#FR
FR_2023_matrix <- FR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

FR_2024_matrix <- FR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

FR_2025_matrix <- FR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#LR
LR_2023_matrix <- LR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

LR_2024_matrix <- LR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

LR_2025_matrix <- LR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#HR
HR_2023_matrix <- HR_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

HR_2024_matrix <- HR_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

HR_2025_matrix <- HR_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

#NT
NT_2023_matrix <- NT_2023 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

NT_2024_matrix <- NT_2024 %>%
  select(all_of(species_cols)) %>%
  as.matrix()

NT_2025_matrix <- NT_2025 %>%
  select(all_of(species_cols)) %>%
  as.matrix()


## Calculation of Beta Diversity Contrasts

### Contrast 1: 2023 vs 2024
NT_2023_2024 <- beta.temp(NT_2023_matrix, NT_2024_matrix, index.family="jaccard")
FR_2023_2024 <- beta.temp(FR_2023_matrix, FR_2024_matrix, index.family="jaccard")
LR_2023_2024 <- beta.temp(LR_2023_matrix, LR_2024_matrix, index.family="jaccard")
HR_2023_2024 <- beta.temp(HR_2023_matrix, HR_2024_matrix, index.family="jaccard")

# Add treatment and contrast group identifier
all_treatments_2023_2024_UD <- bind_rows(
  NT_2023_2024 %>% mutate(treatment = "NT"),
  FR_2023_2024 %>% mutate(treatment = "FR"),
  LR_2023_2024 %>% mutate(treatment = "LR"),
  HR_2023_2024 %>% mutate(treatment = "HR")
) %>% mutate(Contrast_Group = "2023_vs_2024")


### Contrast 2: 2023 vs 2025
NT_2023_2025 <- beta.temp(NT_2023_matrix, NT_2025_matrix, index.family="jaccard")
FR_2023_2025 <- beta.temp(FR_2023_matrix, FR_2025_matrix, index.family="jaccard")
LR_2023_2025 <- beta.temp(LR_2023_matrix, LR_2025_matrix, index.family="jaccard")
HR_2023_2025 <- beta.temp(HR_2023_matrix, HR_2025_matrix, index.family="jaccard")

# Add treatment and contrast group identifier
all_treatments_2023_2025_UD <- bind_rows(
  NT_2023_2025 %>% mutate(treatment = "NT"),
  FR_2023_2025 %>% mutate(treatment = "FR"),
  LR_2023_2025 %>% mutate(treatment = "LR"),
  HR_2023_2025 %>% mutate(treatment = "HR")
) %>% mutate(Contrast_Group = "2023_vs_2025")




## ----------------------------------------------------------------
## Averaging the Contrasts for a Final Before vs. After Estimate
## ----------------------------------------------------------------

# 1. Combine the two contrast tables (No change)
final_beta_data_UD <- bind_rows(all_treatments_2023_2024_UD, all_treatments_2023_2025_UD)

# --- NEW: Calculate 95% CIs based on the Standard Error of the Mean (SEM) ---
# Since N = 12 (6 sites * 2 contrasts) per treatment, df = 11.
# The 95% t-score for df=11 is qt(0.975, 11) ≈ 2.201. We'll use this.
t_score <- qt(0.975, 11) # Degrees of Freedom = N - 1 = 12 - 1 = 11

# 2. Summarise mean components and calculate CIs
beta_summary_FINAL_UD <- final_beta_data_UD %>%
  group_by(treatment) %>%
  summarise(
    # Turnover Component (beta.jtu)
    mean_turnover = mean(beta.jtu, na.rm = TRUE),
    sd_turnover = sd(beta.jtu, na.rm = TRUE),
    # Nestedness Component (beta.jne)
    mean_nestedness = mean(beta.jne, na.rm = TRUE),
    sd_nestedness = sd(beta.jne, na.rm = TRUE),
    # Total Beta (beta.jac)
    mean_total = mean(beta.jac, na.rm = TRUE),
    sd_total = sd(beta.jac, na.rm = TRUE),
    n_contrasts = n(), # Should be 12 for all treatments
    .groups = "drop"
  ) %>%
  # Calculate 95% CI bounds for the TOTAL beta diversity (mean_total)
  mutate(
    SE_total = sd_total / sqrt(n_contrasts),
    LCL_total = mean_total - t_score * SE_total,
    UCL_total = mean_total + t_score * SE_total
  )


# 3. Prepare for plotting: Pivot to long format (No change to logic)
beta_long_FINAL_UD <- beta_summary_FINAL_UD %>%
  select(treatment, mean_turnover, mean_nestedness, mean_total, LCL_total, UCL_total) %>%
  pivot_longer(
    cols = c(mean_turnover, mean_nestedness),
    names_to = "Component",
    values_to = "Mean_beta"
  ) %>%
  mutate(
    Component = recode(
      Component,
      mean_turnover = "Species replacement",
      mean_nestedness = "Species loss"
    ),
    # Ensure treatment order is logical for the plot
    treatment = factor(treatment, levels = c("NT", "FR", "HR", "LR"))
  )

# 4. Final Stacked Bar Plot with CIs on Total Mean
p_beta_stack_FINAL_UD <- ggplot(beta_long_FINAL_UD,
                                aes(x = treatment,
                                    y = Mean_beta,
                                    fill = Component)) +
  geom_col(color = "black") +
  # --- FIX APPLIED HERE: Added x = treatment to the aes() mapping ---
  geom_errorbar(aes(x = treatment, # <-- ADDED X-MAPPING
                    ymin = LCL_total,
                    ymax = UCL_total, 
                    y = mean_total),
                width = 0.2, 
                linewidth = 0.8,
                color = "black",
                inherit.aes = FALSE, 
                data = beta_summary_FINAL_UD
  ) +
  labs(
    x = "Treatment",
    y = "Mean β-diversity (Jaccard family)",
    fill = "Component",
    title = "β-Diversity: Mean Change (2023 vs. 2024/2025)(Unlimited Distance)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p_beta_stack_FINAL_UD)








## ================================================================
## Community trajectories in NMDS space (unlimited distance example)
## ================================================================

# ... (NMDS section remains exactly as it was, using all three years for trajectory plotting) ...

## 1. Build community matrix and metadata ---------------------------

# Site × species matrix (can be presence/absence or abundance)
comm_mat_UD <- dat5_div_UD %>%
  select(all_of(species_cols)) %>%
  as.matrix()

# Metadata for each row in the matrix
meta_UD <- dat5_div_UD %>%
  select(location, treatment, year) %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    year      = factor(year, levels = c("2023", "2024", "2025"))
  )

## 2. NMDS ordination (Bray–Curtis) --------------------------------

set.seed(123)
# Note: Using 'bray' on P/A data is equivalent to Jaccard-like measures.
nmds_UD <- metaMDS(comm_mat_UD, distance = "bray", k = 2, trymax = 1000)

# Site scores
site_scores_UD <- as.data.frame(scores(nmds_UD, display = "sites")) %>%
  mutate(
    location  = meta_UD$location,
    treatment = meta_UD$treatment,
    year      = meta_UD$year
  )

## 3. Treatment × year centroids in ordination space ---------------

centroids_UD <- site_scores_UD %>%
  group_by(treatment, year) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  ) %>%
  arrange(treatment, year)

# Optional: trajectory length (amount of compositional change)
traj_length_UD <- centroids_UD %>%
  arrange(treatment, year) %>%
  group_by(treatment) %>%
  summarise(
    traj_length = sum(
      sqrt(diff(NMDS1)^2 + diff(NMDS2)^2)
    ),
    .groups = "drop"
  )
print(traj_length_UD)

## 4. Trajectory plot (arrows between years per treatment) ---------

gg_traj_UD <- ggplot(centroids_UD,
                     aes(x = NMDS1, y = NMDS2,
                         colour = treatment, group = treatment)) +
  # arrows showing direction of change across years within treatment
  geom_path(arrow = arrow(type = "closed", length = unit(0.2, "cm")),
            linewidth = 0.8) +
  # centroid points
  geom_point(size = 3) +
  # year labels next to points
  geom_text(aes(label = year),
            nudge_y = 0.03, size = 3) +
  scale_colour_brewer(palette = "Dark2",
                      name = "Treatment") +
  coord_equal() +
  labs(
    title = "Community trajectories in NMDS space (Unlimited Distance)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(gg_traj_UD)


























