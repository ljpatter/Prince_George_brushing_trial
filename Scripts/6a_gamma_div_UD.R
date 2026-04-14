# ---
# title: "Gamma diversity (UD): standardized Before/After + BACI contrasts vs NT"
# author: "Leonard Patterson"
# created: "2025-12-15"
# ---

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(iNEXT)
})

## ========================
## Load + preprocess 
## ========================

abundance_all_years_UD <- read.csv("Output/Tabular Data/mean_count_all_years_UD.csv")

species_cols <- names(abundance_all_years_UD)[grepl("^[A-Z]{4}$", names(abundance_all_years_UD))]
stopifnot(length(species_cols) > 0)

abundance_all_years_UD <- abundance_all_years_UD %>%
  mutate(
    year = as.character(year),
    treatment = sub("^[^-]+-([^-]+)-.*$", "\\1", location)
  ) %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(species_cols), ~ replace_na(.x, 0))) %>%
  select(location, treatment, year, all_of(species_cols)) %>%
  group_by(location, treatment, year) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

## ============================================================
## 1) Helper: compute standardized gamma "cell" estimates + BACI
##    Returns a LIST: list(cell=..., baci=...)
## ============================================================

compute_baci_UD <- function(dat) {
  sp_cols <- names(dat)[grepl("^[A-Z]{4}$", names(dat))]
  if (length(sp_cols) == 0) stop("No species columns found (expected 4-letter codes).")
  
  dat <- dat %>%
    mutate(year = as.character(year)) %>%
    mutate(across(all_of(sp_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(sp_cols), ~ replace_na(.x, 0)))
  
  abund_before <- dat %>% filter(year == "2023")
  abund_after  <- dat %>% filter(year %in% c("2024", "2025"))
  
  empty_out <- list(cell = tibble(), baci = tibble())
  if (nrow(abund_before) == 0 || nrow(abund_after) == 0) return(empty_out)
  
  sum_before <- abund_before %>%
    group_by(treatment) %>%
    summarise(across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  sum_after <- abund_after %>%
    group_by(treatment) %>%
    summarise(across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  if (nrow(sum_before) == 0 || nrow(sum_after) == 0) return(empty_out)
  
  to_iNext_list <- function(df) {
    df %>%
      nest(.by = treatment) %>%
      mutate(vec = map(data, ~ .x %>% select(all_of(sp_cols)) %>% unlist(use.names = FALSE) %>% as.numeric())) %>%
      transmute(treatment, vec) %>%
      deframe() %>%
      lapply(\(x) x[x > 0]) # drop zero-abundance species
  }
  
  iNext_before <- to_iNext_list(sum_before)
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
  
  # ---- standardized gamma cell estimates (ALL treatments incl NT) ----
  cell <- equal_size %>%
    transmute(
      Treatment = Assemblage,
      Order.q   = as.integer(Order.q),
      Period,
      qD,
      LCL = qD.LCL,
      UCL = qD.UCL
    )
  
  # If NT missing, we cannot compute BACI contrasts, but we CAN return cell
  if (!("NT" %in% cell$Treatment)) return(list(cell = cell, baci = tibble()))
  
  # ---- BACI contrasts: (After-Before)_Treat - (After-Before)_NT ----
  df0 <- cell %>% mutate(SE = (UCL - LCL) / (2 * 1.96), Var = SE^2)
  
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
    filter(Treatment %in% c("FR", "LR", "HR")) %>%
    left_join(ctrl, by = "Order.q") %>%
    transmute(
      Order.q   = factor(Order.q, levels = c(0, 1, 2), labels = c("q = 0", "q = 1", "q = 2")),
      Treatment = factor(Treatment, levels = c("FR", "LR", "HR")),
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
## 3) Gamma EMM estimates table (q=0/1): Before/After/Total
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

## ========================
## 4) BACI contrasts table 
## ========================

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

## ===============================
## Gamma BACI contrasts WITH CIs
## ===============================

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
    # keep raw bootstrap draws as "x" conceptually (Estimate column untouched here)
    LCL = quantile(Estimate, 0.025, na.rm = TRUE),
    UCL = quantile(Estimate, 0.975, na.rm = TRUE),
    Est_mean = mean(Estimate, na.rm = TRUE),
    p.raw = 2 * pmin(
      mean(Estimate <= 0, na.rm = TRUE),
      mean(Estimate >= 0, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    CI_half = (UCL - LCL) / 2,
    `Estimate ± half CI` = paste0(round(Est_mean, 2), " \u00B1 ", round(CI_half, 2)),
    `p-value`  = signif(p.raw, 4),
    `p (Holm)` = signif(p.adjust(p.raw, method = "holm"), 3)
  ) %>%
  transmute(
    `Diversity order` = Metric,
    Treatment,
    `Estimate ± half CI`,
    Estimate = round(Est_mean, 3),
    `LCL (95% CI)` = round(LCL, 3),
    `UCL (95% CI)` = round(UCL, 3),
    `p-value`,
    `p (Holm)`,
    `CI half-width` = round(CI_half, 2)
  ) %>%
  arrange(`Diversity order`, Treatment)


write.csv(gamma_baci_table_UD, "Output/Tables/GD_UD_BACI_long.csv", row.names = FALSE)

## ========================
## 5) FINAL report table 
## ========================

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








## ========================
## 6. Plot BACI contrasts 
## ========================

# Build plotting df from your already-correct table
plot_df_gamma_UD <- gamma_baci_table_UD %>%
  transmute(
    Order.q = factor(`Diversity order`,
                     levels = c("Species richness", "Shannon diversity"),
                     labels = c("Species richness", "Shannon Diversity")),
    Treatment = factor(Treatment, levels = c("FR", "LR", "HR")),
    Estimate = Estimate,
    LCL = `LCL (95% CI)`,
    UCL = `UCL (95% CI)`
  )

# Symmetric x-axis limits around 0
x_lim <- max(abs(c(plot_df_gamma_UD$LCL, plot_df_gamma_UD$UCL)), na.rm = TRUE)

p_baci_gamma_UD <- ggplot(
  plot_df_gamma_UD,
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

print(p_baci_gamma_UD)

# Save (cowplot saves transparent by default; force white via ggsave)
ggsave(
  filename    = "Figures/Gamma diversity/gamma_div_UD.tiff",
  plot        = p_baci_gamma_UD,
  device      = "tiff",
  width       = 28,
  height      = 26,
  units       = "cm",
  dpi         = 600,
  bg          = "white",
  compression = "lzw")
