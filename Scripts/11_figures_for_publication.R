# ---
# title: "Figures for publication"
# author: "Leonard Patterson"
# created: "2025-12-12"
# description: 
# ---


# Load packages
library(cowplot)
library(ggplot2)


######### ALPHA DIV


# 1) Extract legend (from UD)
legend_alpha <- get_legend(
  p_baci_alpha_UD +
    theme(
      legend.position   = "right",
      legend.box.margin = margin(0, 0, 0, 12)
    )
)

# 2) Left plot (LA): remove legend + KILL inner right margin
p_alpha_LA_noleg <- p_baci_alpha_LA +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 1, b = 5, l = 5)   # <-- key
  )

# 3) Right plot (UD): hide y-axis but keep space + KILL inner left margin
p_alpha_UD_noleg <- p_baci_alpha_UD +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "transparent"),
    axis.text.y  = element_text(color = "transparent"),
    axis.ticks.y = element_line(color = "transparent"),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 1)  # <-- key
  )

# 4) Side-by-side layout (unchanged structure)
alpha_panels <- plot_grid(
  p_alpha_LA_noleg,
  p_alpha_UD_noleg,
  ncol  = 2,
  align = "h",
  axis  = "tb"
)

# 5) Add single legend
p_alpha_combined <- plot_grid(
  alpha_panels,
  legend_alpha,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

# 6) Draw
p_alpha_combined

# Save png
ggsave(
  filename = "Files/Manuscript/Figures/Figure_2.png",
  plot     = p_alpha_combined,
  width    = 10,
  height   = 6,
  dpi      = 600,
  bg       = "white"
)

# Save tiff
ggsave(
  filename = "Files/Manuscript/Figures/Figure_2.tiff",
  plot     = p_alpha_combined,
  device      = "tiff",
  width    = 10,
  height   = 6,
  dpi      = 600,
  bg       = "white"
)







######### GAMMA DIV 

# 1) Extract legend (from UD)
legend_gamma <- get_legend(
  p_baci_gamma_UD +
    theme(
      legend.position   = "right",
      legend.box.margin = margin(0, 0, 0, 12)
    )
)

# 2) Left plot (LA): remove legend + KILL inner right margin
p_gamma_LA_noleg <- p_baci_gamma_LA +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 1, b = 5, l = 5)   # <-- key
  )

# 3) Right plot (UD): hide y-axis but keep space + KILL inner left margin
p_gamma_UD_noleg <- p_baci_gamma_UD +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "transparent"),
    axis.text.y  = element_text(color = "transparent"),
    axis.ticks.y = element_line(color = "transparent"),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 1)  # <-- key
  )

# 4) Side-by-side layout (unchanged structure)
gamma_panels <- plot_grid(
  p_gamma_LA_noleg,
  p_gamma_UD_noleg,
  ncol  = 2,
  align = "h",
  axis  = "tb"
)

# 5) Add single legend
p_gamma_combined <- plot_grid(
  gamma_panels,
  legend_gamma,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

# 6) Draw
p_gamma_combined

# Save png
ggsave(
  filename = "Files/Manuscript/Figures/Figure_3.png",
  plot     = p_gamma_combined,
  width    = 10,
  height   = 6,
  dpi      = 300,
  bg       = "white"
)

# Save
ggsave(
  filename = "Files/Manuscript/Figures/Figure_3.tiff",
  plot     = p_gamma_combined,
  device      = "tiff",
  width    = 10,
  height   = 6,
  dpi      = 600,
  bg       = "white"
)











######### BETA DIV

# 1) Extract legend (from UD)
legend_beta <- get_legend(
  p_beta_stack_FINAL_UD +
    theme(
      legend.position   = "right",
      legend.box.margin = margin(0, 0, 0, 12)
    )
)

# 2) Left plot (LA): no legend + remove INNER right margin
p_beta_LA_noleg <- p_beta_stack_FINAL_LA +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 5, r = 1, b = 5, l = 5)
  )

# 3) Right plot (UD): hide y-axis but KEEP its space + remove INNER left margin
p_beta_UD_noleg <- p_beta_stack_FINAL_UD +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "transparent"),
    axis.text.y  = element_text(color = "transparent"),
    axis.ticks.y = element_line(color = "transparent"),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 1)
  )

# 4) Align + combine
aligned_beta <- align_plots(
  p_beta_LA_noleg,
  p_beta_UD_noleg,
  align = "h",
  axis  = "tb"
)

beta_panels <- plot_grid(
  aligned_beta[[1]],
  aligned_beta[[2]],
  ncol = 2,
  rel_widths = c(1, 1)
)

# 5) Add single legend
p_beta_combined <- plot_grid(
  beta_panels,
  legend_beta,
  ncol = 2,
  rel_widths = c(1.0, 0.25)
)

# 6) Draw
p_beta_combined

# Save png
ggsave(
  filename = "Files/Manuscript/Figures/Figure_4.png",
  plot     = p_beta_combined,
  width    = 10,
  height   = 6,
  dpi      = 600,
  bg       = "white"
)

# Save
ggsave(
  filename = "Files/Manuscript/Figures/Figure_4.tiff",
  plot     = p_beta_combined,
  device      = "tiff",
  width    = 10,
  height   = 6,
  dpi      = 600,
  bg       = "white"
)










### SSM

# 1) Stack with an explicit spacer row
p_contrast_stacked <- plot_grid(
  contrast_plot_log_treat_LA,
  NULL,                        # <-- spacer
  contrast_plot_log_treat_UD,
  ncol = 1,
  rel_heights = c(1, 0.03, 1)  # <-- controls gap size
)

# 2) Add white padding around the entire figure
p_contrast_final <- ggdraw(p_contrast_stacked) +
  theme(
    plot.margin = margin(12, 12, 12, 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Draw
p_contrast_final

# Save tiff
ggsave(
  filename = "Files/Manuscript/Figures/Figure_5.tiff",
  plot     = p_contrast_final,
  device      = "tiff",
  width    = 7,
  height   = 10,   # taller for vertical layout
  dpi      = 600,
  bg       = "white"
)

# Save png
ggsave(
  filename = "Files/Manuscript/Figures/Figure_5.png",
  plot     = p_contrast_final,
  width    = 7,
  height   = 10,   # taller for vertical layout
  dpi      = 600,
  bg       = "white")









### SSM EMMs



# 1) Stack with an explicit spacer row
p_EMM_stacked <- plot_grid(
  p_ssm_emm_LA,
  NULL,                        # <-- spacer
  p_ssm_emm_UD,
  ncol = 1,
  rel_heights = c(1, 0.03, 1)  # <-- controls gap size
)

# 2) Add white padding around the entire figure
p_EMM_final <- ggdraw(p_EMM_stacked) +
  theme(
    plot.margin = margin(12, 12, 12, 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Draw
p_EMM_final

# Save tiff
ggsave(
  filename = "Files/Manuscript/Figures/Figure_6.tiff",
  plot     = p_EMM_final,
  device      = "tiff",
  width    = 7,
  height   = 10,   # taller for vertical layout
  dpi      = 600,
  bg       = "white"
)

# Save png
ggsave(
  filename = "Files/Manuscript/Figures/Figure_6.png",
  plot     = p_EMM_final,
  width    = 8,
  height   = 10,   # taller for vertical layout
  dpi      = 600,
  bg       = "white")

