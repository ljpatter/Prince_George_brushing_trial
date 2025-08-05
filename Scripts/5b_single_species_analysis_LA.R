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

# Create Site_ID column for random effects structure (site nested within block)
dat3_LA <- dat2_LA %>%
  mutate(site = paste(block, treatment, sep = "-"))

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
    # Define site and block as factors for random effects
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













### CREATE FUNCTION FOR SINGLE-SPECIES MODELS

# Initialize an empty data frame to collect all BACI contrasts
all_contrasts_df_LA <- data.frame()

# --- ROBUST HELPER FUNCTION FOR LOG CONTRAST CALCULATION ---
calculate_baci_log_contrast <- function(model, species_name) {
  
  # Get treatment x time-period marginal means on the log scale
  emm <- emmeans(model, ~ treatment * time_period, type = "link")
  
  # BACI contrasts using the row order:
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
  
  # Summary table
  sum_tbl <- as.data.frame(summary(baci))
  
  # Confidence intervals
  ci_tbl <- as.data.frame(confint(baci))
  
  # Standardize CI column names
  if ("lower.CL" %in% names(ci_tbl)) {
    ci_tbl$LCL <- ci_tbl$lower.CL
    ci_tbl$UCL <- ci_tbl$upper.CL
  } else if (".lower" %in% names(ci_tbl)) {
    ci_tbl$LCL <- ci_tbl$.lower
    ci_tbl$UCL <- ci_tbl$.upper
  } else if ("asymp.LCL" %in% names(ci_tbl)) {
    ci_tbl$LCL <- ci_tbl$asymp.LCL
    ci_tbl$UCL <- ci_tbl$asymp.UCL
  } else {
    stop("Could not find CI columns in confint(baci).")
  }
  
  # Join summary + CI
  out <- dplyr::left_join(
    sum_tbl[, c("contrast", "estimate", "p.value")],
    ci_tbl[, c("contrast", "LCL", "UCL")],
    by = "contrast"
  ) %>%
    dplyr::transmute(
      Species = species_name,
      Treatment = contrast,
      LogContrast = estimate,
      LCL = LCL,
      UCL = UCL,
      PValue = p.value
    ) %>%
    dplyr::arrange(Species, Treatment)
  
  return(out)
}








###### MODELs
### FORMULA: Response ~ treatment * time_period + (1|site) + (1|block)

# =============================
# 1. AMRE (American Redstart) 
# =============================
AMRE_model_LA <- glmmTMB(
  AMRE ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)

AMRE_contrasts_LA <- calculate_baci_log_contrast(AMRE_model_LA, "AMRE")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, AMRE_contrasts_LA)
summary(AMRE_model_LA)
AICc(AMRE_model_LA)

### Assess model residuals

set.seed(123)  # for reproducibility of the simulations
res_AMRE_LA <- simulateResiduals(
  fittedModel = AMRE_model_LA,
  n = 1000    
)
plot(res_AMRE_LA)


# =========
# 2. ALFL  
# =========
ALFL_model_LA <- glmmTMB(
  ALFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)
ALFL_contrasts_LA <- calculate_baci_log_contrast(ALFL_model_LA, "ALFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, ALFL_contrasts_LA)
AICc(ALFL_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_ALFL_LA <- simulateResiduals(
  fittedModel = ALFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_ALFL_LA)


# ==============================================================================
# 3. DEJU (Dark-eyed Junco) | VMR: 1.03 | Family: Poisson
# ==============================================================================
DEJU_model_LA <- glmmTMB(
  DEJU ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_LA
)
DEJU_contrasts_LA <- calculate_baci_log_contrast(DEJU_model_LA, "DEJU")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, DEJU_contrasts_LA)
AICc(DEJU_model_LA)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DEJU_LA <- simulateResiduals(
  fittedModel = DEJU_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DEJU_LA)


# ==============================================================================
# 4. DUFL (Dusky Flycatcher) | VMR: 0.74 | Family: Compois
# ==============================================================================
DUFL_model_LA <- glmmTMB(
  DUFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)

DUFL_contrasts_LA <- calculate_baci_log_contrast(DUFL_model_LA, "DUFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, DUFL_contrasts_LA)
AICc(DUFL_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_DUFL_LA <- simulateResiduals(
  fittedModel = DUFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_DUFL_LA)


# ==============================================================================
# 5. LISP (Lincoln's Sparrow) | VMR: 0.91 | Family: Com-Poisson
# ==============================================================================
LISP_model_LA <- glmmTMB(
  LISP ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_LA
)
LISP_contrasts_LA <- calculate_baci_log_contrast(LISP_model_LA, "LISP")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, LISP_contrasts_LA)
summary(LISP_model_LA)
AICc(LISP_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_LISP_LA <- simulateResiduals(
  fittedModel = LISP_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_LISP_LA)


# ==============================================================================
# 6. OCWA (Orange-crowned Warbler) | VMR: 0.92 | Family: Com-Poisson
# ==============================================================================
OCWA_model_LA <- glmmTMB(
  OCWA ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)
OCWA_contrasts_LA <- calculate_baci_log_contrast(OCWA_model_LA, "OCWA")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, OCWA_contrasts_LA)
summary(OCWA_model_LA)
AICc(OCWA_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OCWA_LA <- simulateResiduals(
  fittedModel = OCWA_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OCWA_LA)



# ==============================================================================
# 7. OSFL (OSFL) | VMR: 0.94 | Family: Poisson
# ==============================================================================
OSFL_model_LA <- glmmTMB(
  OSFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_LA
)
OSFL_contrasts_LA <- calculate_baci_log_contrast(OSFL_model_LA, "OSFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, OSFL_contrasts_LA)
summary(OSFL_model_LA)
AICc(OSFL_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OSFL_LA <- simulateResiduals(
  fittedModel = OSFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_OSFL_LA)


# ==============================================================================
# 8. SWTH (Swainson's Thrush) | VMR: 1.03 | Family: Poisson
# ==============================================================================
SWTH_model_LA <- glmmTMB(
  SWTH ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"),
  data = dat4_pooled_LA
)

SWTH_contrasts_LA <- calculate_baci_log_contrast(SWTH_model_LA, "SWTH")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, SWTH_contrasts_LA)
AICc(SWTH_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_SWTH_LA <- simulateResiduals(
  fittedModel = SWTH_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_SWTH_LA)


# ==============================================================================
# 9. WAVI | VMR: 1.05 | Family: Com-Poisson
# ==============================================================================
WAVI_model_LA <- glmmTMB(
  WAVI ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)

WAVI_contrasts_LA <- calculate_baci_log_contrast(WAVI_model_LA, "WAVI")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, WAVI_contrasts_LA)
summary(WAVI_model_LA)
AICc(WAVI_model_LA)



# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WAVI_LA <- simulateResiduals(
  fittedModel = WAVI_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WAVI_LA)


# ==============================================================================
# 10. WTSP (White-throated Sparrow) | VMR: 1.30 | Family: Poisson
# ==============================================================================
WTSP_model_LA <- glmmTMB(
  WTSP ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"),
  data = dat4_pooled_LA
)

WTSP_contrasts_LA <- calculate_baci_log_contrast(WTSP_model_LA, "WTSP")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, WTSP_contrasts_LA)
summary(WTSP_model_LA)
AICc(WTSP_model_LA)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_WTSP_LA <- simulateResiduals(
  fittedModel = WTSP_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_WTSP_LA)


# ==============================================================================
# 11. YRWA (Yellow-rumped Warbler) | VMR: 0.81 | Family: Com-Poisson
# ==============================================================================
YRWA_model_LA <- glmmTMB(
  YRWA ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_LA
)
YRWA_contrasts_LA <- calculate_baci_log_contrast(YRWA_model_LA, "YRWA")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, YRWA_contrasts_LA)
AICc(YRWA_model_LA)


# Residuals
set.seed(123)  # for reproducibility of the simulations
res_YRWA_LA <- simulateResiduals(
  fittedModel = YRWA_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_YRWA_LA)


# Arrange results for a cleaner printout - NOW USING LOGCONTRAST
print(all_contrasts_df_LA %>% 
        mutate(across(c(LogContrast, LCL, UCL), ~round(., 3))) %>%
        arrange(Species, Treatment))
















############ Generate EMM for each treatment

get_log_emm_table <- function(model, species_name) {
  
  emm <- emmeans(model, ~ treatment * time_period, type = "link")
  
  # Summary (log scale): emmean + SE
  sum_tbl <- summary(emm) %>%
    as_tibble() %>%
    select(treatment, time_period, emmean, SE)
  
  # CI (log scale)
  ci_tbl <- confint(emm) %>%
    as_tibble()
  
  # Standardize CI column names robustly
  if (all(c("lower.CL", "upper.CL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = lower.CL, UCL = upper.CL)
  } else if (all(c(".lower", ".upper") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = .lower, UCL = .upper)
  } else if (all(c("asymp.LCL", "asymp.UCL") %in% names(ci_tbl))) {
    ci_tbl <- ci_tbl %>% rename(LCL = asymp.LCL, UCL = asymp.UCL)
  } else if (!all(c("LCL", "UCL") %in% names(ci_tbl))) {
    stop("Could not find CI columns in confint(emm) output. Names were: ",
         paste(names(ci_tbl), collapse = ", "))
  }
  
  ci_tbl <- ci_tbl %>% select(treatment, time_period, LCL, UCL)
  
  # Join + label
  out <- sum_tbl %>%
    left_join(ci_tbl, by = c("treatment", "time_period")) %>%
    transmute(
      Species     = species_name,
      treatment   = treatment,
      time_period = time_period,
      LogMean     = emmean,
      SE          = SE,
      LCL         = LCL,
      UCL         = UCL
    ) %>%
    arrange(Species, treatment, time_period)
  
  out
}

## ---------- Put your fitted models in a named list ----------
model_list_LA <- list(
  AMRE = AMRE_model_LA,
  ALFL = ALFL_model_LA,
  DEJU = DEJU_model_LA,
  DUFL = DUFL_model_LA,
  OCWA = OCWA_model_LA,
  SWTH = SWTH_model_LA,
  WAVI = WAVI_model_LA,
  WTSP = WTSP_model_LA,
  YRWA = YRWA_model_LA
)

## ---------- Extract EMMs for ALL species (and keep going if one fails) ----------
all_emmeans_log_LA <- purrr::imap_dfr(
  model_list_LA,
  ~ tryCatch(
    get_log_emm_table(.x, .y),
    error = function(e) {
      message("FAILED for ", .y, ": ", conditionMessage(e))
      tibble(
        Species     = .y,
        treatment   = factor(NA, levels = c("NT","LR","HR","FR")),
        time_period = factor(NA, levels = c("Before","After")),
        LogMean     = NA_real_, SE = NA_real_, LCL = NA_real_, UCL = NA_real_
      )
    }
  )
)

## ---------- Print ----------
print(
  all_emmeans_log_LA %>%
    mutate(across(c(LogMean, SE, LCL, UCL), ~ round(.x, 3))) %>%
    arrange(Species, treatment, time_period)
)

## ---------- Save ----------
write.csv(
  all_emmeans_log_LA,
  "Output/Tables/SSM_EMMeans_LogScale_LA.csv",
  row.names = FALSE
)











## ==========================================
## Create table of model outputs
## ==========================================

# Convert log contrasts to IRRs
all_contrasts_IRR_LA <- all_contrasts_df_LA %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>%
  mutate(
    IRR        = exp(LogContrast),
    IRR_LCL    = exp(LCL),
    IRR_UCL    = exp(UCL),
    Delta      = IRR - 1,          # absolute change (multiplicative - 1)
    Delta_pct  = 100 * (IRR - 1)   # percent change
  )

# Tidy summary table for reporting
species_IRR_table_LA <- all_contrasts_IRR_LA %>%
  mutate(
    # Half CI width on log scale
    LogCI_half = (UCL - LCL) / 2,
    
    # Rounded values
    LogContrast = round(LogContrast, 3),
    LogCI_half  = round(LogCI_half, 3),
    LCL         = round(LCL, 3),
    UCL         = round(UCL, 3),
    IRR         = round(IRR, 2),
    IRR_LCL     = round(IRR_LCL, 2),
    IRR_UCL     = round(IRR_UCL, 2),
    Delta_pct   = round(Delta_pct, 1),
    
    Delta_pct_text = dplyr::case_when(
      Delta_pct > 0  ~ paste0(Delta_pct, "% increase"),
      Delta_pct < 0  ~ paste0(abs(Delta_pct), "% decrease"),
      TRUE           ~ "0.0% (no change)"
    ),
    
    Interpretation = case_when(
      IRR_LCL > 1 ~ "Increase (CI > 1)",
      IRR_UCL < 1 ~ "Decrease (CI < 1)",
      TRUE        ~ "NS (includes 1)"
    )
  ) %>%
  arrange(Species, Treatment) %>%
  select(
    Species,
    Treatment,
    LogContrast, LogCI_half, LCL, UCL,
    IRR, IRR_LCL, IRR_UCL,
    Delta_pct,
    Delta_pct_text,
    PValue,
    Interpretation
  )

# Look at the table
species_IRR_table_LA

# Save to CSV
write.csv(
  species_IRR_table_LA,
  "Output/Tables/SSM_BACI_IRR_LA.csv",
  row.names = FALSE
)








########## PLOT RESPONSES BY TREATMENT (log coefficients)

all_contrasts_df_LA_filtered <- all_contrasts_df_LA %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>%  # same exclusions as before
  filter(!is.na(LogContrast)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR")),
    Species   = factor(Species),
    # order species A→Z, then flip so A is at top
    Species   = forcats::fct_rev(forcats::fct_relevel(Species, sort(levels(Species))))
  )

contrast_plot_log_treat_LA <- all_contrasts_df_LA_filtered %>%
  ggplot(aes(x = LogContrast, y = Species)) +
  
  # Reference line at 0 (no BACI effect on log scale)
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 1) +
  
  # Symmetric CIs on log scale
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 height = 0.2, linewidth = 1, color = "black") +
  geom_point(aes(color = Treatment), size = 3) +
  
  facet_wrap(~ Treatment, ncol = 3, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +   # adjust limits if needed
  
  labs(
    title = "Limited Amplitude",
    x = "BACI contrasts (log scale)",
    y = "Species"
  ) +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y      = element_text(face = "italic"),
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1)
  )


print("--- Plotting BACI Contrasts by Treatment (log scale) ---")
print(contrast_plot_log_treat_LA)

#cowplot::save_plot("Figures/Species-level models/BACI_SSM_LA_logcoef.png",
#                   contrast_plot_log_treat_LA, base_width = 7, base_height = 5)


out_fig <- "Figures/Species-level models/BACI_SSM_LA_logcoef.tiff"
ggsave(
  filename    = out_fig,
  plot        = contrast_plot_log_treat_LA,
  device      = "tiff",
  width       = 42,
  height      = 18,
  units       = "cm",
  dpi         = 600,
  compression = "lzw")













########## PLOT MARGINAL MEANS BY TREATMENT (FULL SPECIES NAMES)

# Species code -> full name lookup
species_lookup <- c(
  YRWA = "Yellow-rumped warbler",
  AMRE = "American redstart",
  ALFL = "Alder Flycatcher",
  DEJU = "Dark-eyed junco",
  DUFL = "Dusky flycatcher",
  OCWA = "Orange-crowned warbler",
  OSFL = "Olive-sided flycatcher",
  SWTH = "Swainson's thrush",
  WAVI = "Warbling vireo",
  WTSP = "White-throated sparrow"
)

# Convert log-scale EMMs to response scale for plotting
all_emmeans_resp_LA <- all_emmeans_log_LA %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>%
  filter(!is.na(LogMean)) %>%
  mutate(
    Species     = unname(species_lookup[as.character(Species)]),
    treatment   = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After")),
    Mean        = exp(LogMean),
    LCL_resp    = exp(LCL),
    UCL_resp    = exp(UCL)
  )

# Plot 
p_ssm_emm_LA <- ggplot(
  all_emmeans_resp_LA,
  aes(x = treatment, y = Mean, color = time_period, group = time_period)
) +
  geom_errorbar(
    aes(ymin = LCL_resp, ymax = UCL_resp, group = time_period),
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.5
  ) +
  facet_wrap(~ Species, ncol = 3, scales = "free_y") +
  labs(
    title = "Limited Amplitude",
    x     = "Treatment",
    y     = "Estimated marginal mean",
    color = "Time period"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    strip.text       = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "right",
    plot.background  = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x      = element_text(hjust = 1, color = "black"),
    axis.text.y      = element_text(color = "black"),
    axis.title.x     = element_text(color = "black"),
    axis.title.y     = element_text(color = "black")
  )

print("--- Plotting marginal means by treatment and time period ---")
print(p_ssm_emm_LA)

cowplot::save_plot(
  "Figures/Species-level models/SSM_LA_marginal_means.png",
  p_ssm_emm_LA,
  base_width = 12,
  base_height = 12
)


















############# EMM for individual treatments

# 1. Prepare data
all_emmeans_filtered <- all_emmeans_log_LA %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>%
  mutate(
    treatment = factor(treatment, levels = c("NT", "LR", "HR", "FR")),
    time_period = factor(time_period, levels = c("Before", "After"))
  )

# 2. Create the Forest Plot
emm_forest_plot_LA <- ggplot(all_emmeans_filtered, 
                             aes(x = treatment, y = LogMean, color = time_period)) +
  
  # Horizontal reference line at 0 (optional, helps ground the log scale)
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  
  # Use position_dodge to prevent dots and lines from overlapping
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2, linewidth = 0.8) +
  
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  
  # Facet by species
  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
  
  labs(
    title = "Limited Amplitude",
    x = "Treatment",
    y = "Estimated Log Mean Count (± 95% CI)",
    color = "Time Period"
  ) +
  
  # High-contrast colors for Before vs After
  scale_color_brewer(palette = "Dark2") +
  
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "italic", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "gray90", fill = NA),
    legend.position = "bottom"
  )

# 3. Print and Save
print(emm_forest_plot_LA)

cowplot::save_plot("Figures/Species-level models/EMM_ForestPlot_LA.png", 
                   emm_forest_plot_LA, base_width = 13, base_height = 11)







bird_spp <- mvabund(dat4_pooled_LA[, 5:25])

meanvar.plot(bird_spp)

plot(bird_spp ~ as.factor(dat4_pooled_LA$treatment), cex.axis = 0.8, cex = 0.8)

mod1 <- manyglm(bird_spp ~ treatment*time_period, family = "poisson", data = dat4_pooled_LA)
plot(mod1)

mod2 <- manyglm(bird_spp ~ treatment*time_period, family = "negative_binomial", data = dat4_pooled_LA)
plot(mod2)

# 
anova(mod1, p.uni = "adjusted")






### Test code

View(dat4_pooled_LA)

# =============================
# 1. AMRE (American Redstart) 
# =============================
AMRE_model_LA <- glmmTMB(
  AMRE ~ treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)

AMRE_contrasts_LA <- calculate_baci_log_contrast(AMRE_model_LA, "AMRE")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, AMRE_contrasts_LA)
summary(AMRE_model_LA)
AICc(AMRE_model_LA)

### Assess model residuals

set.seed(123)  # for reproducibility of the simulations
res_AMRE_LA <- simulateResiduals(
  fittedModel = AMRE_model_LA,
  n = 1000    
)
plot(res_AMRE_LA)



### AMRE EMMS
emm_AMRE_LA <- emmeans(
  AMRE_model_LA,
  ~ treatment * time_period,
  type = "link"
)

# View table
emm_AMRE_LA_table <- as.data.frame(emm_AMRE_LA)
print(emm_AMRE_LA_table)


### AMRE BACI contrasts
baci_AMRE_LA <- contrast(emm_AMRE_LA, method = list(
    LR_vs_NT = c( 1, -1,  0,  0, -1,  1,  0,  0),
    HR_vs_NT = c( 1,  0, -1,  0, -1,  0,  1,  0),
    FR_vs_NT = c( 1,  0,  0, -1, -1,  0,  0,  1)))

# Summary on link scale
baci_AMRE_LA_df <- as.data.frame(summary(baci_AMRE_LA))
print(baci_AMRE_LA_df)

# Calculate IRRs
baci_AMRE_LA_df_IRR <- baci_AMRE_LA_df %>%
  mutate(
    IRR = exp(estimate))
print(baci_AMRE_LA_df_IRR)


## Calculate manually
print(emm_AMRE_LA_table)

# EMM output
treatment time_period     emmean        SE  df asymp.LCL  asymp.UCL
NT        Before      -0.4961477 0.4964102 Inf -1.469094  0.4767985
LR        Before      -2.0111074 0.6535253 Inf -3.291993 -0.7302214
HR        Before      -1.3048450 0.5694162 Inf -2.420880 -0.1888098
FR        Before      -2.6146107 0.8245873 Inf -4.230772 -0.9984494
NT        After       -0.6442686 0.4829297 Inf -1.590793  0.3022562
LR        After       -1.4862392 0.5752937 Inf -2.613794 -0.3586842
HR        After       -1.4685787 0.5382572 Inf -2.523543 -0.4136140
FR        After       -1.5357192 0.5476752 Inf -2.609143 -0.4622956

# FR contrast (After-treatment - before-treatment) - (After-control - Before-control)
FR_BACI_AMRE <- (-1.5357192 - (-2.6146107)) - (-0.6442686 - (-0.4961477)) 
print(FR_BACI_AMRE)
exp(FR_BACI_AMRE)





# =========
# 2. ALFL  
# =========
ALFL_model_LA <- glmmTMB(
  ALFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_LA
)
ALFL_contrasts_LA <- calculate_baci_log_contrast(ALFL_model_LA, "ALFL")
all_contrasts_df_LA <- bind_rows(all_contrasts_df_LA, ALFL_contrasts_LA)
AICc(ALFL_model_LA)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_ALFL_LA <- simulateResiduals(
  fittedModel = ALFL_model_LA,
  n = 1000    # number of simulations; increase if you want more stable tests
)
plot(res_ALFL_LA)


