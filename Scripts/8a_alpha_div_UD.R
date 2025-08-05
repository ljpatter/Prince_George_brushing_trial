# ---
# title: "Alpha diversity using UD"
# author: "Leonard Patterson"
# created: "2025-12-12"
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
library(cowplot)   
library(DHARMa) 
library(MuMIn)


rm(list = ls())

# Data 
dat1_alp_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

# Identify species columns in the filtered dataset
species_cols <- setdiff(names(dat1_alp_UD), c("location", "year", "recording_date_time"))

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
alpha_dat_1 <- alpha_dat %>%
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

# Some sites are spread over different two blocks so block-site structure doesn't currently reflect
# nested structure. Hard code correct block-site structure
alpha_dat_2 <- alpha_dat_1 %>%
  mutate(
    block = case_when(
      # still update all MUS113 blocks
      block == "MUS113" ~ "MUS124",
      # only this specific site gets moved to MUS061A
      block == "MUS124" & site == "MUS124-NT-C" ~ "MUS061A",
      TRUE ~ block
    )
  )

######### Fit GLMMs for alpha Hill numbers

### q = 0 (alpha richness)
alpha_q0_model_UD <- glmmTMB(
  alpha_q0 ~ treatment  + time_period + treatment*time_period + (1 | block/site),
  family = compois(link = "log"),
  data   = alpha_dat_2
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
  alpha_q1 ~ treatment  + time_period + treatment*time_period + (1 | block/site),
  family = Gamma(link = "log"),
  data   = alpha_dat_2
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
  alpha_q2 ~ treatment  + time_period + treatment*time_period + (1 | block/site),
  family = Gamma(link = "log"),
  data   = alpha_dat_2
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
##    on the RESPONSE (Hill-number) scale
## ==========================================

get_baci_UD <- function(emm_obj, order_q) {
  
  # Coefficient vectors defined for this row order:
  # 1 NT Before, 2 LR Before, 3 HR Before, 4 FR Before,
  # 5 NT After, 6 LR After, 7 HR After, 8 FR After
  L <- list(
    LR = c( 1, -1,  0,  0, -1,  1,  0,  0),
    HR = c( 1,  0, -1,  0, -1,  0,  1,  0),
    FR = c( 1,  0,  0, -1, -1,  0,  0,  1)
  )
  
  # 1) Regrid to RESPONSE scale (effective # of species)
  emm_resp <- regrid(emm_obj, transform = "response")
  
  # 2) BACI contrasts on response scale
  out <- contrast(emm_resp, method = L)
  
  # 3) Get estimates + tests + CIs (this includes p-values)
  out_df <- as.data.frame(summary(out, infer = TRUE))  # has p.value and CI cols
  names(out_df) <- sub("^p\\.value$", "p.value", names(out_df)) # safety
  
  # Robustly detect CI column names (emmeans version differences)
  lcl_name <- if ("lower.CL"  %in% names(out_df)) "lower.CL" else
    if ("asymp.LCL" %in% names(out_df)) "asymp.LCL" else
      if (".lower" %in% names(out_df)) ".lower" else NA_character_
  
  ucl_name <- if ("upper.CL"  %in% names(out_df)) "upper.CL" else
    if ("asymp.UCL" %in% names(out_df)) "asymp.UCL" else
      if (".upper" %in% names(out_df)) ".upper" else NA_character_
  
  if (is.na(lcl_name) || is.na(ucl_name)) {
    stop("Could not find LCL/UCL columns in summary(out, infer=TRUE). Names were: ",
         paste(names(out_df), collapse = ", "))
  }
  
  out_df %>%
    dplyr::transmute(
      Order.q   = order_q,
      Treatment = contrast,
      Estimate  = estimate,
      SE        = SE,
      df        = if ("df" %in% names(out_df)) df else NA_real_,
      stat      = if ("t.ratio" %in% names(out_df)) t.ratio else
        if ("z.ratio" %in% names(out_df)) z.ratio else NA_real_,
      p.value   = p.value,
      LCL       = .data[[lcl_name]],
      UCL       = .data[[ucl_name]]
    )
}


# Apply to each q
baci_q0_UD <- get_baci_UD(emm_q0_UD, 0)
baci_q1_UD <- get_baci_UD(emm_q1_UD, 1)
baci_q2_UD <- get_baci_UD(emm_q2_UD, 2)

# Combine and set factor orders for plotting
baci_results_alpha_UD <- bind_rows(baci_q0_UD, baci_q1_UD) %>%
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
  baci_results_alpha_UD$LCL,
  baci_results_alpha_UD$UCL
)), na.rm = TRUE)

p_baci_alpha_UD <- ggplot(
  baci_results_alpha_UD,
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
  scale_x_continuous(limits = c(-x_lim, x_lim)) +   
  scale_color_brewer(
    palette = "Dark2",
    breaks  = c("HR", "LR", "FR")
  ) +
  labs(
    title = "Unlimited Distance",
    x     = "Change in Effective Number of Species",
    y     = "Treatment",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position      = "right",
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_blank(),
    strip.text           = element_text(face = "bold"),
    plot.title           = element_text(hjust = 0.5),
    plot.title.position  = "plot"
  )

print(p_baci_alpha_UD)

# Save plot
ggsave(
  filename    = "Figures/Alpha diversity/alpha_div_UD.tiff",
  plot        = p_baci_alpha_UD,
  device      = "tiff",
  width       = 28,
  height      = 26,
  units       = "cm",
  dpi         = 600,
  bg          = "white",
  compression = "lzw")








## ============================================================
## Summarise alpha BACI contrasts into a table (Unlimited Distance)
## ============================================================

baci_table_alpha_UD <- baci_results_alpha_UD %>%
  mutate(
    CI_width    = UCL - LCL,
    CI_half     = CI_width / 2,
    Estimate_pm = paste0(
      round(Estimate, 2), " \u00B1 ", round(CI_half, 2)
    ),
    Interpretation = case_when(
      LCL > 0  ~ "Increase (CI > 0)",
      UCL < 0  ~ "Decrease (CI < 0)",
      TRUE     ~ "NS (includes 0)"
    )
  ) %>%
  arrange(
    Order.q,
    factor(Treatment, levels = c("FR", "LR", "HR"))
  ) %>%
  transmute(
    `Diversity order` = Order.q,
    Treatment,
    `Estimate ± half CI` = Estimate_pm,   # <-- NEW (reporting-ready)
    Estimate            = round(Estimate, 3),
    `LCL (95% CI)`      = round(LCL, 3),
    `UCL (95% CI)`      = round(UCL, 3),
    `p-value`           = round(p.value, 3),
    `CI width`          = round(CI_width, 2),
    `CI half-width`     = round(CI_half, 2),  # <-- NEW (numeric)
    `Interpretation`    = Interpretation
  )

# Look at the table in R
baci_table_alpha_UD

# Save to CSV
write.csv(
  baci_table_alpha_UD,
  "Output/Tables/AD_UD_table.csv",
  row.names = FALSE
)










## ============================================================
## Combine EMMs (q0 & q1) + build reporting table (After/Before/Grand Total)
## (Unlimited Distance)  -- SHORT VERSION
## ============================================================

emm_to_table <- function(emm, metric) {
  df <- as.data.frame(confint(regrid(emm, transform = "response")))
  lcl <- c("lower.CL","asymp.LCL",".lower")[c("lower.CL","asymp.LCL",".lower") %in% names(df)][1]
  ucl <- c("upper.CL","asymp.UCL",".upper")[c("upper.CL","asymp.UCL",".upper") %in% names(df)][1]
  if (is.na(lcl) || is.na(ucl)) stop("CI cols not found: ", paste(names(df), collapse=", "))
  df %>%
    transmute(
      Metric = metric,
      treatment,
      time_period,
      Mean = response,
      CI_half = (.data[[ucl]] - .data[[lcl]])/2,
      Cell = paste0(round(response, 2), " \u00B1 ", round((.data[[ucl]] - .data[[lcl]])/2, 2))
    )
}

emm_table_alpha_UD <- bind_rows(
  emm_to_table(emm_q0_UD, "Species richness (q=0)"),
  emm_to_table(emm_q1_UD, "Shannon diversity (q=1)")
) %>%
  mutate(
    Metric = factor(Metric, levels = c("Species richness (q=0)", "Shannon diversity (q=1)")),
    time_period = factor(time_period, levels = c("After","Before")),
    treatment = factor(treatment, levels = c("NT","HR","LR","FR"))
  )

# --- Wide table with After / Before / Grand Total (within each Metric x Treatment) ---
wide_cells <- emm_table_alpha_UD %>%
  select(Metric, treatment, time_period, Cell) %>%
  pivot_wider(names_from = time_period, values_from = Cell) %>%
  left_join(
    emm_table_alpha_UD %>%
      group_by(Metric, treatment) %>%
      summarise(Mean = mean(Mean), CI_half = mean(CI_half), .groups="drop") %>%
      mutate(`Grand Total` = paste0(round(Mean,2), " \u00B1 ", round(CI_half,2))) %>%
      select(Metric, treatment, `Grand Total`),
    by = c("Metric","treatment")
  ) %>%
  arrange(Metric, treatment)

# --- Add Metric header rows + final Grand Total row ---
emm_report_table_UD <- wide_cells %>%
  mutate(RowType = "treatment",
         Row     = as.character(treatment)) %>%
  select(Metric, RowType, Row, After, Before, `Grand Total`) %>%
  bind_rows(
    wide_cells %>%
      distinct(Metric) %>%
      transmute(
        Metric,
        RowType = "header",
        Row     = as.character(Metric),
        After = "", Before = "", `Grand Total` = ""
      )
  ) %>%
  arrange(
    Metric,
    factor(RowType, levels = c("header", "treatment")),
    factor(Row, levels = c("FR","LR","HR","NT"))
  ) %>%
  select(Row, Before, After, `Grand Total`) %>%
  bind_rows(tibble(Row="Grand Total", After="", Before="", `Grand Total`=""))

emm_report_table_UD

# Save
write.csv(emm_report_table_UD, "Output/Tables/AD_UD_EMM_summary_table.csv", row.names = FALSE)

