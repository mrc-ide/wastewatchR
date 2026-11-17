suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
})

output_dir <- "inst/analysis_scripts/outbreak_simulator_comparison"
replicates <- read.csv(file.path(output_dir, "replicate_detection_results.csv"))
summary <- read.csv(file.path(output_dir, "detection_frequency_summary.csv"))
replicates$parameter_set <- factor(replicates$parameter_set,
  levels = paste0("R0 = ", c(0.2, 0.4, 0.6, 0.8, 0.95)))
summary$parameter_set <- factor(summary$parameter_set,
  levels = levels(replicates$parameter_set))
method_levels <- c("Simple single outbreak", "Evolutionary single outbreak",
                   "Multi-outbreak simulator")
replicates$method <- factor(replicates$method, levels = method_levels)
summary$method <- factor(summary$method, levels = method_levels)

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())
p_size <- ggplot(replicates,
  aes(x = parameter_set, y = infections_in_window, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.82), outlier.alpha = 0.15,
               width = 0.72) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Infections during the three-year window",
       y = "Number of infections (log scale)", fill = NULL) +
  base_theme + theme(legend.position = "none")
p_clinical <- ggplot(summary,
  aes(x = parameter_set, y = clinical_rate, color = method, group = method)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = clinical_lower, ymax = clinical_upper), width = 0.1) +
  scale_color_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  labs(title = "Clinical detection frequency", y = "Replicates detected",
       color = NULL) + base_theme
legend <- get_legend(p_clinical + theme(legend.position = "bottom"))
figure <- plot_grid(p_size, p_clinical + theme(legend.position = "none"),
                    ncol = 1, align = "v", rel_heights = c(1.2, 1))
figure <- plot_grid(figure, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(file.path(output_dir, "outbreak_size_clinical_comparison.png"),
       figure, width = 11, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "outbreak_size_clinical_comparison.pdf"),
       figure, width = 11, height = 8, bg = "white")
