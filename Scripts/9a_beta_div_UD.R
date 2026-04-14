# ---
# title: "Beta diversity UD"
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
dat_div_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

# Create treatment column without remove original location column
dat1_div_UD <- dat_div_UD %>%
  mutate(site = location)

dat2_div_UD <- dat1_div_UD %>%
  separate_wider_delim(site, delim = "-", names = c("block", "treatment", "plot"))

# Identify species columns
species_cols <- setdiff(names(dat2_div_UD), c("location", "year", "block", "treatment", "plot", "recording_date_time"))

# Remove unneeded columns
dat3_div_UD <- dat2_div_UD %>%
  select(location, treatment, year, all_of(species_cols))

# Condense to one row per site-year by summing counts across surveys
dat4_div_UD <- dat3_div_UD %>%
  group_by(location, treatment, year) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

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

# Calculate 95% CIs based on the Standard Error of the Mean (SEM) ---
# Since N = 12 (6 sites * 2 contrasts) per treatment, df = 11.
t_score <- qt(0.975, 11) # Degrees of Freedom = N - 1 = 12 - 1 = 11

# 2. Summarise mean components and calculate CIs
beta_summary_FINAL_UD <- final_beta_data_UD %>%
  group_by(treatment) %>%
  summarise(
    # Turnover Component (beta.jtu)
    mean_turnover = mean(beta.jtu, na.rm = TRUE),
    sd_turnover   = sd(beta.jtu, na.rm = TRUE),
    
    # Nestedness Component (beta.jne)
    mean_nestedness = mean(beta.jne, na.rm = TRUE),
    sd_nestedness   = sd(beta.jne, na.rm = TRUE),
    
    # Total Beta (beta.jac)
    mean_total = mean(beta.jac, na.rm = TRUE),
    sd_total   = sd(beta.jac, na.rm = TRUE),
    
    n_contrasts = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Standard errors
    SE_turnover     = sd_turnover   / sqrt(n_contrasts),
    SE_nestedness   = sd_nestedness / sqrt(n_contrasts),
    SE_total        = sd_total      / sqrt(n_contrasts),
    
    # 95% confidence intervals for each component
    LCL_turnover   = mean_turnover   - t_score * SE_turnover,
    UCL_turnover   = mean_turnover   + t_score * SE_turnover,
    
    LCL_nestedness = mean_nestedness - t_score * SE_nestedness,
    UCL_nestedness = mean_nestedness + t_score * SE_nestedness,
    
    LCL_total      = mean_total      - t_score * SE_total,
    UCL_total      = mean_total      + t_score * SE_total
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
    treatment = factor(treatment, levels = c("NT", "HR", "LR", "FR"))
  )


# 4. Final Stacked Bar Plot with CIs on Total Mean
p_beta_stack_FINAL_UD <- ggplot(
  beta_long_FINAL_UD,
  aes(x = treatment, y = Mean_beta, fill = Component)
) +
  geom_col(color = "black") +
  geom_errorbar(
    aes(
      x = treatment,
      ymin = LCL_total,
      ymax = UCL_total,
      y = mean_total
    ),
    width = 0.2,
    linewidth = 0.8,
    color = "black",
    inherit.aes = FALSE,
    data = beta_summary_FINAL_UD
  ) +
  labs(
    x = "Treatment",
    y = "Mean β-diversity",
    fill = "",
    title = "      Unlimited Distance"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.title.position = "plot"
  )

print(p_beta_stack_FINAL_UD)

# Save plot
ggsave(
  filename    = "Figures/Beta diversity/beta_div_UD.tiff",
  plot        = p_beta_stack_FINAL_UD,
  device      = "tiff",
  width       = 28,
  height      = 26,
  units       = "cm",
  dpi         = 600,
  bg          = "white",
  compression = "lzw")





## ====================================
## 5. β-diversity results tables (UD)
## ====================================

# ---- 5a) LONG table ----
beta_table_FINAL_UD <- beta_summary_FINAL_UD %>%
  mutate(
    treatment = factor(treatment, levels = c("NT","FR","LR","HR")),
    
    CI_half_turnover   = (UCL_turnover   - LCL_turnover)   / 2,
    CI_half_nestedness = (UCL_nestedness - LCL_nestedness) / 2,
    CI_half_total      = (UCL_total      - LCL_total)      / 2,
    
    Turnover_pm   = paste0(round(mean_turnover,   3), " \u00B1 ", round(CI_half_turnover,   3)),
    Nestedness_pm = paste0(round(mean_nestedness, 3), " \u00B1 ", round(CI_half_nestedness, 3)),
    Total_pm      = paste0(round(mean_total,      3), " \u00B1 ", round(CI_half_total,      3))
  ) %>%
  transmute(
    Treatment = as.character(treatment),
    
    Metric = "Species replacement (β.jtu)",
    `Mean ± half CI` = Turnover_pm,
    Mean = round(mean_turnover, 3),
    `LCL (95% CI)` = round(LCL_turnover, 3),
    `UCL (95% CI)` = round(UCL_turnover, 3),
    `CI half-width` = round(CI_half_turnover, 3)
  ) %>%
  bind_rows(
    beta_summary_FINAL_UD %>%
      mutate(
        treatment = factor(treatment, levels = c("NT","FR","LR","HR")),
        CI_half_nestedness = (UCL_nestedness - LCL_nestedness) / 2,
        Nestedness_pm = paste0(round(mean_nestedness, 3), " \u00B1 ", round(CI_half_nestedness, 3))
      ) %>%
      transmute(
        Treatment = as.character(treatment),
        Metric = "Species loss (β.jne)",
        `Mean ± half CI` = Nestedness_pm,
        Mean = round(mean_nestedness, 3),
        `LCL (95% CI)` = round(LCL_nestedness, 3),
        `UCL (95% CI)` = round(UCL_nestedness, 3),
        `CI half-width` = round(CI_half_nestedness, 3)
      ),
    beta_summary_FINAL_UD %>%
      mutate(
        treatment = factor(treatment, levels = c("NT","FR","LR","HR")),
        CI_half_total = (UCL_total - LCL_total) / 2,
        Total_pm = paste0(round(mean_total, 3), " \u00B1 ", round(CI_half_total, 3))
      ) %>%
      transmute(
        Treatment = as.character(treatment),
        Metric = "Total β (β.jac)",
        `Mean ± half CI` = Total_pm,
        Mean = round(mean_total, 3),
        `LCL (95% CI)` = round(LCL_total, 3),
        `UCL (95% CI)` = round(UCL_total, 3),
        `CI half-width` = round(CI_half_total, 3)
      )
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Species replacement (β.jtu)",
                                       "Species loss (β.jne)",
                                       "Total β (β.jac)")),
    Treatment = factor(Treatment, levels = c("FR","LR","HR","NT"))
  ) %>%
  arrange(Metric, Treatment)

beta_table_FINAL_UD

write.csv(beta_table_FINAL_UD,
          "Output/Tables/Beta_diversity_UD_Table_long.csv",
          row.names = FALSE)


# ---- 5b) INTERLEAVED table (metric header row then treatments) ----
beta_table_FINAL_UD_interleaved <- beta_table_FINAL_UD %>%
  group_by(Metric) %>%
  group_modify(~ bind_rows(
    tibble(Treatment = NA_character_, `Mean ± half CI` = NA_character_,
           Mean = NA_real_, `LCL (95% CI)` = NA_real_, `UCL (95% CI)` = NA_real_,
           `CI half-width` = NA_real_),
    .x
  )) %>%
  ungroup() %>%
  mutate(
    Row = if_else(is.na(Treatment), as.character(Metric), as.character(Treatment))
  ) %>%
  select(Row, `Mean ± half CI`, Mean, `LCL (95% CI)`, `UCL (95% CI)`, `CI half-width`)

beta_table_FINAL_UD_interleaved

write.csv(beta_table_FINAL_UD_interleaved,
          "Output/Tables/Beta_diversity_UD_Table_interleaved.csv",
          row.names = FALSE)

