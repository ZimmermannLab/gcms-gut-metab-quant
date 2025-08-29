# ===============================
# Required libraries and theme
# ===============================
library(tidyverse)
library(ggpubr)

theme_set(
  theme_minimal(base_family = "Helvetica", base_size = 7) +
    theme(
      axis.text       = element_text(family = "Helvetica", size = 7),
      axis.title      = element_text(family = "Helvetica", size = 7),
      legend.text     = element_text(family = "Helvetica", size = 7),
      legend.title    = element_text(family = "Helvetica", size = 7),
      plot.title      = element_text(family = "Helvetica", size = 7, hjust = 0.5),
      strip.text      = element_text(family = "Helvetica", size = 7)
    )
)

# ===============================
# Load and preprocess data
# ===============================
data <- read_csv2("data_results.csv") %>%   # <-- Replace with relative/public path
  separate(data_filename, into = c("mice", "tissue"), sep = "_") %>% 
  separate(tissue, into = c("tissue"), sep = "\\.") %>%
  filter(!tissue %in% c("liver", "plasma"))

# Define tissue order
data$tissue <- factor(
  data$tissue, 
  levels = c("duod", "jeju", "ileum", "cecum", "colon", "feces")
)

# ===============================
# Function: Plot metabolite boxplots
# ===============================
plot_metabolite <- function(df, metabolite, show_labels = TRUE, point_size = 3) {
  
  # Filter for metabolite of interest
  metabolite_data <- df %>% filter(name_original == metabolite)
  
  # Compute per-tissue t-tests
  p_values <- metabolite_data %>%
    group_by(tissue) %>%
    summarise(p_value = t.test(concentration_calculated_uM ~ mice)$p.value, .groups = "drop") %>%
    filter(p_value < 0.05)
  
  # Brackets for significance
  signif_comparisons <- p_values %>%
    mutate(
      x1 = as.numeric(tissue) - 0.25,
      x2 = as.numeric(tissue) + 0.25,
      y  = max(metabolite_data$concentration_calculated_uM) * 1.05,
      label = if (show_labels) sprintf("p = %.3f", p_value) else "*"
    )
  
  # Base plot
  p <- ggplot(metabolite_data, aes(x = tissue, y = concentration_calculated_uM, fill = mice)) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
    geom_jitter(
      position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),
      alpha = 0.7, size = point_size, shape = 21, stroke = 0.5, color = "black"
    ) +
    stat_summary(
      fun.data = mean_cl_normal, geom = "errorbar", width = 0.2,
      position = position_dodge(width = 0.75), color = "black"
    ) +
    labs(
      x = "Tissue",
      y = "Concentration (nmol per mg)",
      fill = "Group",
      title = metabolite
    ) +
    scale_fill_manual(values = c("GF" = "#FF9999", "SPF" = "#9999FF")) +
    scale_x_discrete(labels = c(
      "duod"  = "Duodenum",
      "jeju"  = "Jejunum",
      "ileum" = "Ileum",
      "cecum" = "Cecum",
      "colon" = "Colon",
      "feces" = "Feces"
    )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add significance brackets (if any)
  if (nrow(signif_comparisons) > 0) {
    p <- p + geom_bracket(
      data = signif_comparisons,
      aes(xmin = x1, xmax = x2, y.position = y, label = label),
      inherit.aes = FALSE, tip.length = 0.03, size = 0.5, color = "black"
    )
  }
  
  return(p)
}

# ===============================
# Example usage
# ===============================

# Hydroxyhexanoic acid plot
p1 <- plot_metabolite(data, "Hydroxyhexanoic acid", show_labels = FALSE, point_size = 3)
print(p1)

# Serine plot (with labels)
p2 <- plot_metabolite(data, "Serine", show_labels = TRUE, point_size = 3)
print(p2)

# Propionic acid plot (with labels)
p2 <- plot_metabolite(data, "Propionic acid", show_labels = TRUE, point_size = 3)
print(p2)