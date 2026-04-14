# ---
# title: "Singe species analysis - UD - MEAN COUNT"
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
dat1_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

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
spp <- names(dat3_UD)[grepl("^[A-Z]{4}$", names(dat3_UD))]
dat4_UD <- dat3_UD %>%
  select(block, site, treatment, year, all_of(spp))

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

# Some sites are spread over different two blocks so block-site structure doesn't currently reflect
# nested structure. Hard code correct block-site structure

dat4_pooled_UD <- dat4_pooled_UD %>%
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
write.csv(dat4_pooled_UD, "Output/Tabular Data/mean_count_all_years_UD_SSM.csv")






### INDIVIDUAL SPECIES MODELS

# Initialize an empty data frame to collect all BACI contrasts
all_contrasts_df_UD <- data.frame()

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






###### MODELS
### FORMULA: Response ~ treatment * time_period + (1|site) + (1|block)

# ===================================================
# 1. AMRE (American Redstart) | Family: Com-Poisson
# ===================================================
AMRE_model_UD <- glmmTMB(
  AMRE ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

AMRE_contrasts_UD <- calculate_baci_log_contrast(AMRE_model_UD, "AMRE")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, AMRE_contrasts_UD)
AICc(AMRE_model_UD)

### Assess model residuals

#set.seed(123)  # for reproducibility of the simulations
#res_AMRE_UD <- simulateResiduals(
#  fittedModel = AMRE_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_AMRE_UD)


# ===============================
# 2. ALFL | Family: Com-Poisson
# ===============================
ALFL_model_UD <- glmmTMB(
  ALFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
ALFL_contrasts_UD <- calculate_baci_log_contrast(ALFL_model_UD, "ALFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, ALFL_contrasts_UD)
AICc(ALFL_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_ALFL_UD <- simulateResiduals(
#  fittedModel = ALFL_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_ALFL_UD)

# ===============================================
# 3. AMRO (American Robin) | Family: Com-Poisson
# ===============================================
AMRO_model_UD <- glmmTMB(
  AMRO ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
AMRO_contrasts_UD <- calculate_baci_log_contrast(AMRO_model_UD, "AMRO")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, AMRO_contrasts_UD)
AICc(AMRO_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_AMRO_UD <- simulateResiduals(
#  fittedModel = AMRO_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_AMRO_UD)


# ===============================================
# 4. DEJU (Dark-eyed Junco) | Family: Poisson
# ===============================================
DEJU_model_UD <- glmmTMB(
  DEJU ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_UD
)

DEJU_contrasts_UD <- calculate_baci_log_contrast(DEJU_model_UD, "DEJU")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, DEJU_contrasts_UD)
AICc(DEJU_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_DEJU_UD <- simulateResiduals(
#  fittedModel = DEJU_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_DEJU_UD)


# ==================================================
# 5. DUFL (Dusky Flycatcher) | Family: Com-Poisson
# ==================================================
DUFL_model_UD <- glmmTMB(
  DUFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

DUFL_contrasts_UD <- calculate_baci_log_contrast(DUFL_model_UD, "DUFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, DUFL_contrasts_UD)
AICc(DUFL_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_DUFL_UD <- simulateResiduals(
#  fittedModel = DUFL_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_DUFL_UD)

# ================================================================
# 6. LISP (Lincoln's Sparrow) | VMR: 0.82 | Family: Com-Poisson
# ================================================================
LISP_model_UD <- glmmTMB(
  LISP ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
LISP_contrasts_UD <- calculate_baci_log_contrast(LISP_model_UD, "LISP")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, LISP_contrasts_UD)
AICc(LISP_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_LISP_UD <- simulateResiduals(
#  fittedModel = LISP_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_LISP_UD)


# ===========================================================
# 7. OCWA (Orange-crowned Warbler) | Family: ZI-Com-Poisson
# ===========================================================
OCWA_model_UD <- glmmTMB(
  OCWA ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  ziformula = ~ 1,
  data = dat4_pooled_UD
)
OCWA_contrasts_UD <- calculate_baci_log_contrast(OCWA_model_UD, "OCWA")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, OCWA_contrasts_UD)
AICc(OCWA_model_UD)

# Residuals
set.seed(123)  # for reproducibility of the simulations
res_OCWA_UD <- simulateResiduals(
  fittedModel = OCWA_model_UD,
  n = 1000)    # number of simulations; increase if you want more stable tests
plot(res_OCWA_UD)



# ==================================
# 8. OSFL (OSFL) | Family: Poisson
# ==================================
OSFL_model_UD <- glmmTMB(
  OSFL ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_UD
)
OSFL_contrasts_UD <- calculate_baci_log_contrast(OSFL_model_UD, "OSFL")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, OSFL_contrasts_UD)
AICc(OSFL_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_OSFL_UD <- simulateResiduals(
#  fittedModel = OSFL_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_OSFL_UD)


# =====================================================
# 9. RCKI (Ruby-crowned Kinglet) | Family: Com Poisson
# =====================================================
RCKI_model_UD <- glmmTMB(
  RCKI ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
RCKI_contrasts_UD <- calculate_baci_log_contrast(RCKI_model_UD, "RCKI")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, RCKI_contrasts_UD)
AICc(RCKI_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_RCKI_UD <- simulateResiduals(
#  fittedModel = RCKI_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_RCKI_UD)

# ===================================================
# 10. SWTH (Swainson's Thrush) | Family: Com-Poisson
# ===================================================
SWTH_model_UD <- glmmTMB(
  SWTH ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

SWTH_contrasts_UD <- calculate_baci_log_contrast(SWTH_model_UD, "SWTH")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, SWTH_contrasts_UD)
AICc(SWTH_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_SWTH_UD <- simulateResiduals(
#  fittedModel = SWTH_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_SWTH_UD)

# ================================
# 11. WAVI | Family: Com-Poisson
# ================================
WAVI_model_UD <- glmmTMB(
  WAVI ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
WAVI_contrasts_UD <- calculate_baci_log_contrast(WAVI_model_UD, "WAVI")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, WAVI_contrasts_UD)
AICc(WAVI_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_WAVI_UD <- simulateResiduals(
#  fittedModel = WAVI_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_WAVI_UD)

# =============================================
# 12. WETA (Western Tanager) | Family: Poisson
# =============================================
WETA_model_UD <- glmmTMB(
  WETA ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = poisson(link = "log"), 
  data = dat4_pooled_UD
)
WETA_contrasts_UD <- calculate_baci_log_contrast(WETA_model_UD, "WETA")
all_contrasts_df <- bind_rows(all_contrasts_df_UD, WETA_contrasts_UD)
AICc(WETA_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_WETA_UD <- simulateResiduals(
#  fittedModel = WETA_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_WETA_UD)


# =========================================================
# 13. WTSP (White-throated Sparrow) | Family: Com-Poisson
# =========================================================
WTSP_model_UD <- glmmTMB(
  WTSP ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"),
  data = dat4_pooled_UD
)

WTSP_contrasts_UD <- calculate_baci_log_contrast(WTSP_model_UD, "WTSP")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, WTSP_contrasts_UD)
AICc(WTSP_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_WTSP_UD <- simulateResiduals(
#  fittedModel = WTSP_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_WTSP_UD)


# =======================================================
# 14. YRWA (Yellow-rumped Warbler) | Family: Com-Poisson
# =======================================================
YRWA_model_UD <- glmmTMB(
  YRWA ~ treatment + time_period + treatment*time_period + (1|block/site),
  family = compois(link = "log"), 
  data = dat4_pooled_UD
)
YRWA_contrasts_UD <- calculate_baci_log_contrast(YRWA_model_UD, "YRWA")
all_contrasts_df_UD <- bind_rows(all_contrasts_df_UD, YRWA_contrasts_UD)
AICc(YRWA_model_UD)

# Residuals
#set.seed(123)  # for reproducibility of the simulations
#res_YRWA_UD <- simulateResiduals(
#  fittedModel = YRWA_model_UD,
#  n = 1000    # number of simulations; increase if you want more stable tests
#)
#plot(res_YRWA_UD)


# Arrange results for a cleaner printout - NOW USING LOGCONTRAST
print(all_contrasts_df_UD %>% 
        mutate(across(c(LogContrast, LCL, UCL), ~round(., 3))) %>%
        arrange(Species, Treatment))







############ Generate EMM for each treatment

# Create function
get_log_emm_table <- function(model, species_name) {
  
  emm <- emmeans(model, ~ treatment * time_period, type = "link")
  
  # Summary (log scale): emmean + SE
  sum_tbl <- summary(emm) %>%
    as_tibble() %>%
    select(treatment, time_period, emmean, SE)
  
  # CI (log scale): confint() output has version-dependent column names
  ci_tbl <- confint(emm) %>%
    as_tibble()
  
  # Standardize CI column names 
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

## ---------- Put fitted models in a named list ----------
model_list_UD <- list(
  AMRE = AMRE_model_UD,
  ALFL = ALFL_model_UD,
  AMRO = AMRO_model_UD,
  DEJU = DEJU_model_UD,
  DUFL = DUFL_model_UD,
  LISP = LISP_model_UD,
  OCWA = OCWA_model_UD,
  OSFL = OSFL_model_UD,
  RCKI = RCKI_model_UD,
  SWTH = SWTH_model_UD,
  WAVI = WAVI_model_UD,
  WETA = WETA_model_UD,
  WTSP = WTSP_model_UD,
  YRWA = YRWA_model_UD
)

## ---------- Extract EMMs for ALL species ----------
all_emmeans_log_UD <- purrr::imap_dfr(
  model_list_UD,
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
  all_emmeans_log_UD %>%
    mutate(across(c(LogMean, SE, LCL, UCL), ~ round(.x, 3))) %>%
    arrange(Species, treatment, time_period)
)

## ---------- Save ----------
write.csv(
  all_emmeans_log_UD,
  "Output/Tables/SSM_EMMeans_LogScale_UD.csv",
  row.names = FALSE
)





## ==========================================
## Create table of model outputs
## ==========================================

# Convert log contrasts to IRRs
all_contrasts_IRR_UD <- all_contrasts_df_UD %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>%
  mutate(
    IRR        = exp(LogContrast),
    IRR_LCL    = exp(LCL),
    IRR_UCL    = exp(UCL),
    Delta      = IRR - 1,          # absolute change (multiplicative - 1)
    Delta_pct  = 100 * (IRR - 1)   # percent change
  )


# Tidy summary table for reporting
species_IRR_table_UD <- all_contrasts_IRR_UD %>%
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
    
    # Text versions for reporting
    LogContrast_text = paste0(
      LogContrast, " \u00B1 ", LogCI_half
    ),
    
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
    LogContrast_text,      # <- main reported value
    LogContrast, LogCI_half, LCL, UCL,
    IRR, IRR_LCL, IRR_UCL,
    Delta_pct,
    Delta_pct_text,
    PValue,
    Interpretation
  )

# Look at the table
species_IRR_table_UD

# Save to CSV
write.csv(
  species_IRR_table_UD,
  "Output/Tables/SSM_BACI_IRR_UD.csv",
  row.names = FALSE
)






###### PLOT RESPONSES BY TREATMENT (log coefficients) 

# Remove models that threw errors / untrusted CIs
all_contrasts_df_UD_filtered <- all_contrasts_df_UD %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO")) %>% 
  filter(!is.na(LogContrast)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR")),
    Species   = factor(Species),
    # order species A→Z, then flip so A is at top
    Species   = forcats::fct_rev(forcats::fct_relevel(Species, sort(levels(Species))))
  )

contrast_plot_log_treat_UD <- all_contrasts_df_UD_filtered %>%
  ggplot(aes(x = LogContrast, y = Species, color = Treatment)) +
  
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
    title = "Unlimited Distance",
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
print(contrast_plot_log_treat_UD)

cowplot::save_plot("Figures/Species-level models/BACI_SSM_UD_logcoef.png",
                   contrast_plot_log_treat_UD, base_width = 7, base_height = 5)

out_fig <- "Figures/Species-level models/BACI_SSM_UD_logcoef.tiff"
ggsave(
  filename    = out_fig,
  plot        = contrast_plot_log_treat_UD,
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
  SWTH = "Swainson's thrush",
  WAVI = "Warbling vireo",
  WTSP = "White-throated sparrow"
)

# Convert log-scale EMMs to response scale for plotting
all_emmeans_resp_UD <- all_emmeans_log_UD %>%
  filter(!Species %in% c("OSFL", "RCKI", "LISP", "AMRO", "WETA")) %>%
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
p_ssm_emm_UD <- ggplot(
  all_emmeans_resp_UD,
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
    title = "Unlimited Distance",
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
print(p_ssm_emm_UD)

cowplot::save_plot(
  "Figures/Species-level models/SSM_UD_marginal_means.png",
  p_ssm_emm_UD,
  base_width = 12,
  base_height = 12
)

