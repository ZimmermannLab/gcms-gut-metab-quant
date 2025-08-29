# ================================
# Libraries
# ================================
library(tidyverse)
library(viridis)
library(plotly)
library(ggrepel)

# ================================
# Step 1: Prepare scatter Plot Data
# ================================

data <- data %>%
  mutate(
    MouseType = str_split(data_filename, "_", simplify = TRUE)[, 1],
    Tissue    = str_split(data_filename, "_", simplify = TRUE)[, 2],
    Replicate = str_split(data_filename, "_", simplify = TRUE)[, 3]
  )

filtered_data <- data %>%
  filter(Tissue == "liver") %>% 
  mutate(l10_concentration_calculated_uM = log10(concentration_calculated_uM + 1e-4))

stat_data <- filtered_data %>%
  group_by(name_original, MouseType) %>%
  summarize(
    Mean_Concentration = mean(concentration_calculated_uM, na.rm = TRUE),
    SD_Concentration = sd(concentration_calculated_uM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = MouseType, values_from = c(Mean_Concentration, SD_Concentration))

p_values <- filtered_data %>%
  group_by(name_original) %>%
  summarize(
    p_value = tryCatch({
      values_SPF <- l10_concentration_calculated_uM[MouseType == "SPF"]
      values_GF  <- l10_concentration_calculated_uM[MouseType == "GF"]
      
      if (length(values_SPF) > 1 && length(values_GF) > 1) {
        wilcox.test(values_SPF, values_GF, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH"))

wide_concentrations <- stat_data %>%
  left_join(p_values, by = "name_original") %>%
  mutate(Significant = !is.na(p_value_adj) & p_value_adj < 0.1)

# ================================
# Step 2: Merge Grouping Information
# ================================
df <- read.csv2('/g/zimmermann/Members/Nikita/For_sunburst_4.csv')

group_mapping <- df %>%
  mutate(labels = as.character(labels), parents = as.character(parents)) %>%
  dplyr::select(labels, parents) %>%
  rename(name_original = labels, group = parents)

wide_concentrations <- wide_concentrations %>%
  left_join(group_mapping, by = "name_original")

# ================================
# Step 3: Define Custom Group Colors
# ================================
group_color_map <- c(
  "Benzenoids" = "darkviolet",
  "Organic acids" = "deepskyblue3",
  "Lipids" = "orange",           
  "Amino acids" = "chartreuse3",       
  "Amino acids derivatives" = "red",  
  "Imidazopyrimidines" = "coral4", 
  "Pyrimidines" = "bisque2",
  "Furans" = "chartreuse3",
  "Pyridines" = "darkolivegreen", 
  "Flavonoids, isoflavonoids" = "darkgrey",
  "Phenylpropanoids" = "deepskyblue3",
  "Alkaloids" = "coral",
  "Furanoid lignans" = "aquamarine4",
  "Organooxygens" = "darkturquoise",
  "Lignans" ="bisque1",
  "Organonitrogens" = "darksalmon"
)

wide_concentrations <- wide_concentrations %>%
  mutate(group_color = group_color_map[group]) %>%
  mutate(group_color = ifelse(is.na(group_color), "black", group_color))

# ================================
# Step 4: Prepare log10 means
# ================================
pseudo <- 1e-4
wide_concentrations <- wide_concentrations %>%
  mutate(
    log10_SPF = log10(Mean_Concentration_SPF + pseudo),
    log10_GF  = log10(Mean_Concentration_GF + pseudo)
  )

min_val <- min(c(wide_concentrations$log10_SPF, wide_concentrations$log10_GF), na.rm = TRUE) - 0.5
max_val <- max(c(wide_concentrations$log10_SPF, wide_concentrations$log10_GF), na.rm = TRUE) + 0.5

# =============================
# Manually specify plasma metabolites to label
# =============================
manually_label_all <- c(
  "Aspartic acid", "Acetic acid", "Indole-propionic acid", "3-Indoleacetic acid", "Suberic acid", "Valeric acid"
)


# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Aspartic acid" ~ -0.3,
      name_original == "Acetic acid" ~ 0.3,
      name_original == "Indole-propionic acid" ~ 0.5,
      name_original == "3-Indoleacetic acid" ~ 0.2,
      name_original == "Suberic acid" ~ -0.3,
      name_original == "Valeric acid" ~ -0.5,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Aspartic acid" ~ 0.4,
      name_original == "Acetic acid" ~ 0.1,
      name_original == "Indole-propionic acid" ~ 0.5,
      name_original == "3-Indoleacetic acid" ~ 0.4,
      name_original == "Suberic acid" ~ -0.4,
      name_original == "Valeric acid" ~ -0.3,
      TRUE ~ 0
    )
  )

# ================================
# Step 6: Plot
# ================================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 6, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 11,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 60,
    force_pull = 4.2,
    segment.size = 0.8,
    segment.color = "black",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(plot)

# ================================
# Step 7: Save Plot
# ================================
ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_plasma_scatter.pdf", 
       plot = plot, device = "pdf", units = "in", width = 16, height = 16, dpi = 600)




# =============================
# Manually specify liver metabolites to label
# =============================
manually_label_all <- c(
  "12-Hydroxystearic acid", "Acetic acid", "Tryptamine", "Suberic acid", "Hippuric acid"
)

# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "12-Hydroxystearic acid" ~ -0.7,
      name_original == "Acetic acid" ~ 0.1,
      name_original == "Tryptamine" ~ 0.1,
      name_original == "Suberic acid" ~ 0.4,
      name_original == "Hippuric acid" ~ -0.6,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Aspartic acid" ~ -0.7,
      name_original == "Acetic acid" ~ 0.1,
      name_original == "Tryptamine" ~ 0.1,
      name_original == "Suberic acid" ~ -1.8,
      name_original == "Hippuric acid" ~ -0.6,
      TRUE ~ 0
    )
  )

# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 6, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 11,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 60,
    force_pull = 4.2,
    segment.size = 0.8,
    segment.color = "black",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )



# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_liver_scatter.pdf", 
       device = "pdf", units = "in", width = 16, height = 16, dpi = 600)











# =============================
# Ensure group_color is a valid character and no NA values
# =============================
wide_concentrations <- wide_concentrations %>%
  mutate(group_color = ifelse(is.na(group_color), "black", group_color))

# =============================
# Manually specify duodenum metabolites to label
# =============================
manually_label_all <- c(
  "12-Hydroxystearic acid", "3-hydroxybutyric acid", "Hydroxyoctanoic acid", "Acetic acid",
  "2-Ketoglutaric acid", "Tyrosol", "2-Hydroxyisovaleric acid", "2-Ketobutyric acid",
  "4-Hydroxycinnamic acid", "Alanine", "Cinnamoylglycine", "Cystine", "Dihydroferulic acid",
  "Dimethylmalonic acid", "Ectoine", "GABA", "Glycine", "Isocaproic acid",
  "Isovaleric acid", "Leucine", "N-methylnicotin amide", "Nutriacholic acid", "Phenylvaleric acid",
  "Threonine", "Tricarballylic acid", "Uric acid", "Valine", "5-Keto-D-gluconate",
  "Malonic acid", "3-Hydroxyhippuric acid", "Indole-3-lactic acid", "Phthalic acid"
)

# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "12-Hydroxystearic acid" ~ -0.7,
      name_original == "3-hydroxybutyric acid" ~ 1.4,
      name_original == "Hydroxyoctanoic acid" ~ 0.5,
      name_original == "Acetic acid" ~ 0.1,
      name_original == "2-Ketoglutaric acid" ~ -1.8,
      name_original == "Tyrosol" ~ -0.6,
      name_original == "2-Hydroxyisovaleric acid" ~ 5.6,
      name_original == "2-Ketobutyric acid" ~ -0.5,
      name_original == "4-Hydroxycinnamic acid" ~ -0.6,
      name_original == "Alanine" ~ 0.2,
      name_original == "Cinnamoylglycine" ~ -0.12,
      name_original == "Cystine" ~ -0.4,
      name_original == "Dihydroferulic acid" ~ 0.2,
      name_original == "Dimethylmalonic acid" ~ 0.9,
      name_original == "Ectoine" ~ 0.6,
      name_original == "GABA" ~ 0.1,
      name_original == "Glycine" ~ 0.15,
      name_original == "Isocaproic acid" ~ -0.3,
      name_original == "Isovaleric acid" ~ 0.6,
      name_original == "Leucine" ~ 0.3,
      name_original == "N-methylnicotin amide" ~ 4.0,
      name_original == "Nutriacholic acid" ~ 0.6,
      name_original == "Phenylvaleric acid" ~ 0.5,
      name_original == "Threonine" ~ -0.2,
      name_original == "Tricarballylic acid" ~ 0.9,
      name_original == "Uric acid" ~ -0.3,
      name_original == "Valine" ~ -0.2,
      name_original == "5-Keto-D-gluconate" ~ -0.2,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "3-Hydroxyhippuric acid" ~ 0.2,
      name_original == "Indole-3-lactic acid" ~ -0.5,
      name_original == "Phthalic acid" ~ -0.3,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "12-Hydroxystearic acid" ~ 0.1,
      name_original == "3-hydroxybutyric acid" ~ -0.9,
      name_original == "Hydroxyoctanoic acid" ~ -0.4,
      name_original == "Acetic acid" ~ 0.1,
      name_original == "2-Ketoglutaric acid" ~ 0.3,
      name_original == "Tyrosol" ~ -0.6,
      name_original == "2-Hydroxyisovaleric acid" ~ 0.8,
      name_original == "2-Ketobutyric acid" ~ 0.5,
      name_original == "4-Hydroxycinnamic acid" ~ -0.9,
      name_original == "Alanine" ~ 0.7,
      name_original == "Cinnamoylglycine" ~ -0.3,
      name_original == "Cystine" ~ 0.4,
      name_original == "Dihydroferulic acid" ~ -0.9,
      name_original == "Dimethylmalonic acid" ~ -0.4,
      name_original == "Ectoine" ~ -0.01,
      name_original == "GABA" ~ -0.2,
      name_original == "Glycine" ~ 0.4,
      name_original == "Isocaproic acid" ~ 0.3,
      name_original == "Isovaleric acid" ~ 0.1,
      name_original == "Leucine" ~ 0.7,
      name_original == "N-methylnicotin amide" ~ 0.2,
      name_original == "Nutriacholic acid" ~ 0.4,
      name_original == "Phenylvaleric acid" ~ -0.3,
      name_original == "Threonine" ~ 0.7,
      name_original == "Tricarballylic acid" ~ -0.6,
      name_original == "Uric acid" ~ -0.25,
      name_original == "Valine" ~ 1.9,
      name_original == "5-Keto-D-gluconate" ~ 0.2,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "3-Hydroxyhippuric acid" ~ -0.2,
      name_original == "Indole-3-lactic acid" ~ -0.1,
      name_original == "Phthalic acid" ~ -0.3,
      TRUE ~ 0
    )
  )

# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 60,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol/mg, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol/mg, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_duod_scatter.pdf", 
       device = "pdf", units = "in", width = 15, height = 15, dpi = 600)



# =============================
# Ensure group_color is a valid character and no NA values
# =============================
wide_concentrations <- wide_concentrations %>%
  mutate(group_color = ifelse(is.na(group_color), "black", group_color))

# =============================
# Manually specify jejunum metabolites to label
# =============================
manually_label_all <- c(
  "Alanine", "Glutamic acid", "Leucine", "Glycine", "Phenylalanine", "Malic acid", "Acetic acid", "Lactic acid", "Succinic acid",
  "Citramalic acid", "12-Hydroxystearic acid", "Malonic acid", "2-Ketobutyric acid", "Proline", "Aspartic acid",
  "Oxalic acid", "3-hydroxybutyric acid", "Arachidonic acid", "3-Hydroxyanthranilic acid",
  "4-Aminobenzoic acid", "Kynurenine", "Fumaric acid", "2-Ketoisocaproic acid"
)

# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Alanine" ~ 0.8,
      name_original == "Glutamic acid" ~ -0.5,
      name_original == "Leucine" ~ 0.5,
      name_original == "Glycine" ~ -0.7,
      name_original == "Phenylalanine" ~ -1.3,
      name_original == "Malic acid" ~ 0.7,
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Lactic acid" ~ -0.1,
      name_original == "Succinic acid" ~ 0.8,
      name_original == "Citramalic acid" ~ 0.7,
      name_original == "12-Hydroxystearic acid" ~ -0.7,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "2-Ketobutyric acid" ~ -0.2,
      name_original == "Proline" ~ 0.6,
      name_original == "Aspartic acid" ~ 0.9,
      name_original == "Oxalic acid" ~ 0.5,
      name_original == "3-hydroxybutyric acid" ~ 1.4,
      name_original == "Arachidonic acid" ~ 0.4,
      name_original == "3-Hydroxyanthranilic acid" ~ 0.3,
      name_original == "4-Aminobenzoic acid" ~ 0.4,
      name_original == "Kynurenine" ~ -0.6,
      name_original == "Fumaric acid" ~ -0.6,
      name_original == "2-Ketoisocaproic acid" ~ -0.2,
      TRUE ~ 0  # fallback if any missed
    ),
    nudge_y = case_when(
      name_original == "Alanine" ~ 0.3,
      name_original == "Glutamic acid" ~ 0.5,
      name_original == "Leucine" ~ 0.3,
      name_original == "Glycine" ~ 0.3,
      name_original == "Phenylalanine" ~ -0.2,
      name_original == "Malic acid" ~ 0.2,
      name_original == "Acetic acid" ~ 0.3,
      name_original == "Lactic acid" ~ -0.3,
      name_original == "Succinic acid" ~ 0.1,
      name_original == "Citramalic acid" ~ -0.2,
      name_original == "12-Hydroxystearic acid" ~ 0.1,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "2-Ketobutyric acid" ~ 0.2,
      name_original == "Proline" ~ 0.2,
      name_original == "Aspartic acid" ~ 0.5,
      name_original == "Oxalic acid" ~ -0.4,
      name_original == "3-hydroxybutyric acid" ~ -0.9,
      name_original == "Arachidonic acid" ~ 0.4,
      name_original == "3-Hydroxyanthranilic acid" ~ -0.3,
      name_original == "4-Aminobenzoic acid" ~ -0.5,
      name_original == "Kynurenine" ~ 0.3,
      name_original == "Fumaric acid" ~ 0.5,
      name_original == "2-Ketoisocaproic acid" ~ -0.6,
      TRUE ~ 0
    )
  )


# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 50,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol/mg, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol/mg, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_jejunum_scatter.pdf", 
       device = "pdf", units = "in", width = 15, height = 15, dpi = 600)





# Manually specify feces metabolites to label
# =============================
manually_label_all <- c("Acetic acid", "Dimethylsuccinic acid", "Succinic acid", "Valeric acid", "Isovaleric acid",
                        "Butyric acid", "Propionic acid", "Lactic acid", "Enterolactone", "Hydroxyhexanoic acid",
                        "3-Indoleacetic acid", "3-Indoleacetonitrile", "Equol", "Isocaproic acid",
                        "Vanillylmandelic acid", "Tyrosol", "Dimethylmalonic acid", "Cystine", "4-Hydroxyhippuric acid",
                        "Malonic acid", "Isoleucine", "Glycine", "Leucine",
                        "Threonine", "Proline", "Serine", "Proline", "2-Ketobutyric acid",
                        "Methylmalonic acid", "6-Hydroxynicotinic acid", "Dihydroferulic acid"
                        
)

# =============================
# Generate label data with nudges
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Dihydroferulic acid" ~ 0.1,
      name_original == "6-Hydroxynicotinic acid" ~ 0.2,
      name_original == "Dimethylsuccinic acid" ~ 0.45,
      name_original == "Succinic acid" ~ 0.3,
      name_original == "Valeric acid" ~ 0.4,
      name_original == "Isovaleric acid" ~ 0.65,
      name_original == "Butyric acid" ~ 0.4,
      name_original == "Propionic acid" ~ 0.2,
      name_original == "Lactic acid" ~ 0.25,
      name_original == "Enterolactone" ~ 0.2,
      name_original == "Hydroxyhexanoic acid" ~ 0.25,
      name_original == "3-Indoleacetic acid" ~ -0.25,
      name_original == "3-Indoleacetonitrile" ~ 0.2,
      name_original == "Equol" ~ 0.45,
      name_original == "Isocaproic acid" ~ -0.3,
      name_original == "Vanillylmandelic acid" ~ 0.5,
      name_original == "Tyrosol" ~ -0.6,
      name_original == "Dimethylmalonic acid" ~ 0.9,
      name_original == "Cystine" ~ -0.4,
      name_original == "4-Hydroxyhippuric acid" ~ 0.25,
      name_original == "Malonic acid" ~ -0.25,
      name_original == "Isoleucine" ~ -0.4,     # doubled
      name_original == "Glycine" ~ -0.8,        # doubled      # doubled
      name_original == "Leucine" ~ -1.6,        # doubled
      name_original == "Threonine" ~ -0.5,     # doubled
      name_original == "Proline" ~ 0.8,        # doubled
      name_original == "Serine" ~ 1.9,         # doubled
      name_original == "2-Ketobutyric acid" ~ 2.4,
      name_original == "Methylmalonic acid" ~ -0.65,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Dihydroferulic acid" ~ -0.1,
      name_original == "6-Hydroxynicotinic acid" ~ 0.1,
      name_original == "Dimethylsuccinic acid" ~ -0.4,
      name_original == "Succinic acid" ~ -0.4,
      name_original == "Valeric acid" ~ -0.45,
      name_original == "Isovaleric acid" ~ 0.2,
      name_original == "Butyric acid" ~ 0.3,
      name_original == "Propionic acid" ~ 0.4,
      name_original == "Lactic acid" ~ 0.3,
      name_original == "Enterolactone" ~ 0.2,
      name_original == "Hydroxyhexanoic acid" ~ -0.25,
      name_original == "3-Indoleacetic acid" ~ -0.4,
      name_original == "3-Indoleacetonitrile" ~ -0.4,
      name_original == "Equol" ~ 0.15,
      name_original == "Isocaproic acid" ~ 0.3,
      name_original == "Vanillylmandelic acid" ~ -0.2,
      name_original == "Tyrosol" ~ -0.2,
      name_original == "Dimethylmalonic acid" ~ -0.4,
      name_original == "Cystine" ~ 0.4,
      name_original == "4-Hydroxyhippuric acid" ~ 0.4,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "Isoleucine" ~ -1.3,       # doubled
      name_original == "Glycine" ~ 0.6,          # doubled       # doubled
      name_original == "Leucine" ~ 1.4,          # doubled
      name_original == "Threonine" ~ 0.85,        # doubled
      name_original == "Proline" ~ 0.7,          # doubled
      name_original == "Serine" ~ -1.4,          # doubled
      name_original == "2-Ketobutyric acid" ~ 1.5,
      name_original == "Methylmalonic acid" ~ 1.3,
      TRUE ~ 0
    )
  )



# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.55,
    point.padding = 0.55,
    force = 85,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol/mg, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol/mg, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_feces_scatter.pdf", 
       device = "pdf", units = "in", width = 15, height = 15, dpi = 600)



# =============================
# Ensure group_color is a valid character and no NA values
# =============================
wide_concentrations <- wide_concentrations %>%
  mutate(group_color = ifelse(is.na(group_color), "black", group_color))

# =============================
# Manually specify colon metabolites to label
# =============================
manually_label_all <- c("Acetic acid", "Glutamic acid", "Proline", "Threonine", "Lysine",
                        "Tyrosine", "Oxalic acid", "Creatinine", "Methionine", "Phenylalanine",
                        "Isoleucine", "2-ketoglutaric acid", "Malonic acid", "Cystine", "Maleic acid", "Tyrosol",
                        "Aconitic acid", "3-Indoleacetonitrile", "3-Indoleacetic acid", "Indole-propionic acid",
                        "Hydroxyhlutaric acid", "2-3-Dihydroxybenzoic acid", "Heptanoic acid", "Hydroxyhexanoic acid",
                        "Isobutyric acid", "Citramalic acid", "Isovaleric acid", "Enterolactone", "Propionic acid",
                        "Butyric acid", "Succinic acid", "Suberic acid", "Dimethylsuccinic acid", "Fumaric acid",
                        "Vanillylmandelic acid", "6-Hydroxynicotinic acid", "Isocaproic acid", "Citric acid", "Dimethylmalonic acid",
                        "Dimethylsuccinic acid", "Adipic acid", "Valeric acid", "Hydroxyoctanoic acid", "Equol"
)

# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Glutamic acid" ~ 0.2,
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Proline" ~ 0.2,
      name_original == "Threonine" ~ -0.8,
      name_original == "Lysine" ~ 0.2,
      name_original == "Tyrosine" ~ -1.6,
      name_original == "Oxalic acid" ~ 0.4,
      name_original == "Creatinine" ~ -0.5,
      name_original == "Methionine" ~ -0.6,
      name_original == "Phenylalanine" ~ -1.6,
      name_original == "Isoleucine" ~ -0.33,
      name_original == "2-ketoglutaric acid" ~ -1.8,
      name_original == "Malonic acid" ~ -0.3,
      name_original == "Cystine" ~ -0.4,
      name_original == "Maleic acid" ~ 0.2,
      name_original == "Tyrosol" ~ -0.2,
      name_original == "Aconitic acid" ~ -0.3,
      name_original == "3-Indoleacetonitrile" ~ -0.3,
      name_original == "3-Indoleacetic acid" ~ 0.25,
      name_original == "Indole-propionic acid" ~ 0.2,
      name_original == "Hydroxyhlutaric acid" ~ 0.2,
      name_original == "2-3-Dihydroxybenzoic acid" ~ -0.25,
      name_original == "Heptanoic acid" ~ 0.2,
      name_original == "Hydroxyhexanoic acid" ~ 0.6,
      name_original == "Isobutyric acid" ~ 0.2,
      name_original == "Citramalic acid" ~ 0.1,
      name_original == "Isovaleric acid" ~ -0.3,
      name_original == "Enterolactone" ~ -0.2,
      name_original == "Propionic acid" ~ 0.2,
      name_original == "Butyric acid" ~ 0.2,
      name_original == "Succinic acid" ~ 0.4,
      name_original == "Suberic acid" ~ 0.2,
      name_original == "Dimethylsuccinic acid" ~ 0.2,
      name_original == "Fumaric acid" ~ 1.6,
      name_original == "Vanillylmandelic acid" ~ 0.6,
      name_original == "6-Hydroxynicotinic acid" ~ 0.3,
      name_original == "Isocaproic acid" ~ 0.6,
      name_original == "Citric acid" ~ -0.2,
      name_original == "Dimethylmalonic acid" ~ -1.2,
      name_original == "Adipic acid" ~ -0.4,
      name_original == "Valeric acid" ~ 1.5,
      name_original == "Equol" ~ 0.2,
      name_original == "Hydroxyoctanoic acid" ~ 0.5,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Glutamic acid" ~ -0.2,
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Proline" ~ 0.2,
      name_original == "Threonine" ~ 0.7,
      name_original == "Lysine" ~ 0.1,
      name_original == "Tyrosine" ~ 0.5,
      name_original == "Oxalic acid" ~ -0.6,
      name_original == "Creatinine" ~ -0.5,
      name_original == "Methionine" ~ 0.7,
      name_original == "Phenylalanine" ~ -0.2,
      name_original == "Isoleucine" ~ 1.8,
      name_original == "2-ketoglutaric acid" ~ 0.3,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "Cystine" ~ 0.4,
      name_original == "Maleic acid" ~ -0.4,
      name_original == "Tyrosol" ~ -0.2,
      name_original == "Aconitic acid" ~ -0.3,
      name_original == "3-Indoleacetonitrile" ~ 0.2,
      name_original == "3-Indoleacetic acid" ~ 0.4,
      name_original == "Indole-propionic acid" ~ -0.2,
      name_original == "Hydroxyhlutaric acid" ~ -0.2,
      name_original == "2-3-Dihydroxybenzoic acid" ~ -0.2,
      name_original == "Heptanoic acid" ~ -0.4,
      name_original == "Hydroxyhexanoic acid" ~ 1.1,
      name_original == "Isobutyric acid" ~ -0.1,
      name_original == "Citramalic acid" ~ 0.3,
      name_original == "Isovaleric acid" ~ 0.1,
      name_original == "Enterolactone" ~ -0.1,
      name_original == "Propionic acid" ~ -0.2,
      name_original == "Butyric acid" ~ -0.2,
      name_original == "Succinic acid" ~ -0.2,
      name_original == "Suberic acid" ~ 0.4,
      name_original == "Dimethylsuccinic acid" ~ -0.2,
      name_original == "Fumaric acid" ~ -0.25,
      name_original == "Vanillylmandelic acid" ~ -0.3,
      name_original == "6-Hydroxynicotinic acid" ~ -0.4,
      name_original == "Isocaproic acid" ~ -0.5,
      name_original == "Citric acid" ~ 0.2,
      name_original == "Dimethylmalonic acid" ~ -0.4,
      name_original == "Adipic acid" ~ -0.1,
      name_original == "Valeric acid" ~ 0.5,
      name_original == "Equol" ~ 0.2,
      name_original == "Hydroxyoctanoic acid" ~ -0.2,
      TRUE ~ 0
    )
  )


# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 60,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol/mg, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol/mg, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_colon_scatter.pdf", 
       device = "pdf", units = "in", width = 15, height = 15, dpi = 600)





# =============================
# Ensure group_color is a valid character and no NA values
# =============================
wide_concentrations <- wide_concentrations %>%
  mutate(group_color = ifelse(is.na(group_color), "black", group_color))

# =============================
# Manually specify ileum metabolites to label
# =============================
manually_label_all <- c(
  "Alanine", "Glutamic acid", "Leucine", "Glycine", "Phenylalanine", "Malic acid", "Acetic acid", "Lactic acid", "Succinic acid",
  "Citramalic acid", "12-Hydroxystearic acid", "Malonic acid", "2-Ketobutyric acid", "Aspartic acid",
  "Oxalic acid", "3-hydroxybutyric acid",
  "4-Aminobenzoic acid", "Kynurenine", "Fumaric acid", "2-Ketoisocaproic acid", "Citric acid", "Dihydroferulic acid",
  "3-Indoleacetic acid", "Indole-3-lactic acid", "Enterolactone", "Propionic acid"
)

# =============================
# Generate label data with nudges
# =============================
label_data_all <- wide_concentrations %>%
  filter(name_original %in% manually_label_all) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Alanine" ~ 1.2,
      name_original == "Glutamic acid" ~ -1.4,
      name_original == "Leucine" ~ 0.5,
      name_original == "Glycine" ~ -0.7,
      name_original == "Phenylalanine" ~ 1.3,
      name_original == "Malic acid" ~ 0.74,
      name_original == "Acetic acid" ~ 0.2,
      name_original == "Lactic acid" ~ 0.1,
      name_original == "Succinic acid" ~ 0.8,
      name_original == "Citramalic acid" ~ 0.7,
      name_original == "12-Hydroxystearic acid" ~ -0.7,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "2-Ketobutyric acid" ~ -0.2,
      name_original == "Aspartic acid" ~ 1.5,
      name_original == "Oxalic acid" ~ 0.7,
      name_original == "3-hydroxybutyric acid" ~ 1.4,
      name_original == "4-Aminobenzoic acid" ~ -0.4,
      name_original == "Kynurenine" ~ -0.6,
      name_original == "Fumaric acid" ~ -0.6,
      name_original == "2-Ketoisocaproic acid" ~ -0.5,
      name_original == "Citric acid" ~ -0.2,
      name_original == "Dihydroferulic acid" ~ -0.2,
      name_original == "Propionic acid" ~ -0.2,
      name_original == "Enterolactone" ~ -0.2,
      name_original == "Indole-3-lactic acid" ~ -0.2,
      name_original == "3-Indoleacetic acid" ~ -0.2,
      
      TRUE ~ 0  # fallback if any missed
    ),
    nudge_y = case_when(
      name_original == "Alanine" ~ 0.3,
      name_original == "Glutamic acid" ~ 0.6,
      name_original == "Leucine" ~ 0.4,
      name_original == "Glycine" ~ 0.3,
      name_original == "Phenylalanine" ~ 0.3,
      name_original == "Malic acid" ~ -0.2,
      name_original == "Acetic acid" ~ 0.5,
      name_original == "Lactic acid" ~ 0.01,
      name_original == "Succinic acid" ~ -0.4,
      name_original == "Citramalic acid" ~ -0.6,
      name_original == "12-Hydroxystearic acid" ~ 0.1,
      name_original == "Malonic acid" ~ -0.2,
      name_original == "2-Ketobutyric acid" ~ 0.2,
      name_original == "Aspartic acid" ~ 0.5,
      name_original == "Oxalic acid" ~ 0.1,
      name_original == "3-hydroxybutyric acid" ~ -0.9,
      name_original == "4-Aminobenzoic acid" ~ -0.7,
      name_original == "Kynurenine" ~ 0.45,
      name_original == "Fumaric acid" ~ 0.5,
      name_original == "2-Ketoisocaproic acid" ~ -0.6,
      name_original == "Citric acid" ~ -0.35,
      name_original == "Dihydroferulic acid" ~ -0.5,
      name_original == "Propionic acid" ~ -0.2,
      name_original == "Enterolactone" ~ -0.2,
      name_original == "Indole-3-lactic acid" ~ -0.2,
      name_original == "3-Indoleacetic acid" ~ 0.3,
      TRUE ~ 0
    )
  )


# =============================
# Plot
# =============================
plot <- ggplot(wide_concentrations, aes(x = log10_SPF, y = log10_GF)) +
  geom_point(aes(fill = group_color), shape = 21, size = 4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkblue") +
  
  geom_errorbar(data = filter(wide_concentrations, Significant),
                aes(
                  ymin = log10_GF - (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10))),
                  ymax = log10_GF + (SD_Concentration_GF / (ifelse(Mean_Concentration_GF == 0, pseudo, Mean_Concentration_GF) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = filter(wide_concentrations, Significant),
                 aes(
                   xmin = log10_SPF - (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10))),
                   xmax = log10_SPF + (SD_Concentration_SPF / (ifelse(Mean_Concentration_SPF == 0, pseudo, Mean_Concentration_SPF) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data_all,
    aes(label = name_original, color = group_color),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.25,
    force = 50,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data_all$nudge_x,
    nudge_y = label_data_all$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (SPF plasma) (nmol/mg, log10 scale)",
    y = "Mean concentration (GF plasma) (nmol/mg, log10 scale)",
    color = "Metabolite Group"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  
  scale_x_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  scale_y_continuous(
    breaks = seq(floor(min_val), ceiling(max_val), 1),
    minor_breaks = seq(floor(min_val), ceiling(max_val), 0.2),
    labels = function(x) round(10^x, digits = 4)
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(color = "black", hjust = 0.99, vjust = 0.3),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =============================
# Show plot and save
# =============================
#print(plot)

ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/GF_SPF_ileum_scatter.pdf", 
       device = "pdf", units = "in", width = 15, height = 15, dpi = 600)
















# ================================
# Libraries
# ================================
library(tidyverse)
library(ggrepel)
library(viridis)

# ================================
# Step 1: Load and Prepare Data
# ================================
data <- read_csv2("/g/zimmermann/Members/Nikita/GC-MSMS_application/Buffers_results/_calculated_concentrations_ACN.csv")

data <- data %>%
  mutate(
    Buffer = str_split(data_filename, "_", simplify = TRUE)[, 1],
    Replicate = str_split(data_filename, "_", simplify = TRUE)[, 2],
    SampleType = str_split(data_filename, "_", simplify = TRUE)[, 3]
  )

filtered_data <- data %>%
  filter(SampleType == "ACN", Buffer %in% c("FreshlyFrozen", "OMNIgutFrozen"))

# ================================
# Step 2: Summary + Stats with FDR
# ================================
stat_data <- filtered_data %>%
  group_by(name_original, Buffer) %>%
  summarize(
    Mean_Concentration = mean(concentration_calculated_uM, na.rm = TRUE),
    SD_Concentration = sd(concentration_calculated_uM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Buffer, values_from = c(Mean_Concentration, SD_Concentration))

p_values <- filtered_data %>%
  group_by(name_original) %>%
  summarize(
    p_value = tryCatch({
      v1 <- concentration_calculated_uM[Buffer == "FreshlyFrozen"]
      v2 <- concentration_calculated_uM[Buffer == "OMNIgutFrozen"]
      if (length(v1) > 1 && length(v2) > 1) {
        wilcox.test(v1, v2, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_value_fdr = p.adjust(p_value, method = "fdr"))

wide_concentrations <- stat_data %>%
  left_join(p_values, by = "name_original") %>%
  mutate(
    Significant = !is.na(p_value) & p_value < 0.05,
    Significant_FDR = !is.na(p_value_fdr) & p_value_fdr < 0.1
  )

# ================================
# Step 3: Transformations
# ================================
pseudo <- 1e-4

wide_concentrations <- wide_concentrations %>%
  mutate(
    Residual = Mean_Concentration_FreshlyFrozen - Mean_Concentration_OMNIgutFrozen,
    log10_Fresh = log10(ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen)),
    log10_OMNI = log10(ifelse(Mean_Concentration_OMNIgutFrozen == 0, pseudo, Mean_Concentration_OMNIgutFrozen)),
    log10_Residual = log10(abs(Residual) + pseudo) * sign(Residual)
  )

# ================================
# Step 4: Add Group Colors
# ================================
group_df <- read.csv2('/g/zimmermann/Members/Nikita/For_sunburst_4.csv')

group_mapping <- group_df %>%
  mutate(labels = as.character(labels), parents = as.character(parents)) %>%
  dplyr::select(labels, parents) %>%
  rename(name_original = labels, group = parents)

group_color_map <- c(
  "Benzenoids" = "darkviolet",
  "Organic acids" = "deepskyblue3",
  "Lipids" = "orange",           
  "Amino acids" = "chartreuse3",       
  "Amino acids derivatives" = "red",  
  "Imidazopyrimidines" = "coral4", 
  "Pyrimidines" = "bisque2",
  "Furans" = "chartreuse3",
  "Pyridines" = "darkolivegreen", 
  "Flavonoids, isoflavonoids" = "darkgrey",
  "Phenylpropanoids" = "deepskyblue3",
  "Alkaloids" = "coral",
  "Furanoid lignans" = "aquamarine4",
  "Organooxygens" = "darkturquoise",
  "Lignans" ="bisque1",
  "Organonitrogens" = "darksalmon"
)

wide_concentrations <- wide_concentrations %>%
  left_join(group_mapping, by = "name_original") %>%
  mutate(group_color = group_color_map[group],
         group_color = ifelse(is.na(group_color), "black", group_color))

# ================================
# Step 5: Manual Labels Setup
# ================================
manually_label <- c("Syringic acid", "Cystine", "Tyrosol", "Alanine", "Isoleucine", "Isovaleric acid", "Glutamic acid")

label_data <- wide_concentrations %>%
  filter(name_original %in% manually_label) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Tyrosol" ~ 0.2,
      name_original == "Syringic acid" ~ 0.2,
      name_original == "Cystine" ~ -0.4,
      name_original == "Isovaleric acid" ~ -0.3,
      name_original == "Isoleucine" ~ 0.65,
      name_original == "Alanine" ~ 0.65,
      name_original == "Glutamic acid" ~ 0.65,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Tyrosol" ~ -0.5,
      name_original == "Syringic acid" ~ 0.2,
      name_original == "Cystine" ~ -0.4,
      name_original == "Isovaleric acid" ~ -0.6,
      name_original == "Isoleucine" ~ -1.65,
      name_original == "Alanine" ~ -1.65,
      name_original == "Glutamic acid" ~ -1.25,
      TRUE ~ 0
    )
  )

# ================================
# Step 6: Plot
# ================================
min_val <- min(c(wide_concentrations$log10_Fresh, wide_concentrations$log10_OMNI), na.rm = TRUE) - 0.5
max_val <- max(c(wide_concentrations$log10_Fresh, wide_concentrations$log10_OMNI), na.rm = TRUE) + 0.5

plot <- ggplot(wide_concentrations, aes(x = log10_Fresh, y = log10_OMNI)) +
  geom_point(aes(fill = group_color, shape = Significant_FDR), size = 3, color = "black") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +  # triangle for FDR-significant
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  
  # Error bars
  geom_errorbar(data = label_data,
                aes(
                  ymin = log10_OMNI - (SD_Concentration_OMNIgutFrozen / (ifelse(Mean_Concentration_OMNIgutFrozen == 0, pseudo, Mean_Concentration_OMNIgutFrozen) * log(10))),
                  ymax = log10_OMNI + (SD_Concentration_OMNIgutFrozen / (ifelse(Mean_Concentration_OMNIgutFrozen == 0, pseudo, Mean_Concentration_OMNIgutFrozen) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = label_data,
                 aes(
                   xmin = log10_Fresh - (SD_Concentration_FreshlyFrozen / (ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen) * log(10))),
                   xmax = log10_Fresh + (SD_Concentration_FreshlyFrozen / (ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  # Labels
  geom_text_repel(
    data = label_data,
    aes(label = name_original, color = group_color),
    size = 4,
    box.padding = 0.1,
    point.padding = 0.1,
    force = 825,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data$nudge_x,
    nudge_y = label_data$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (Freshly Frozen, log10 scale)",
    y = "Mean concentration (OMNIgut Frozen, log10 scale)"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    breaks = seq(0, max_val, 1),
    labels = parse(text = paste0("10^", seq(0, max_val, 1)))
  ) +
  scale_y_continuous(
    breaks = seq(0, max_val, 1),
    labels = parse(text = paste0("10^", seq(0, max_val, 1)))
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 0.99, vjust = 0.99),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(plot)

# ================================
# Step 7: Save Plot
# ================================
ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/OMNIgut_scatter.pdf", 
       plot = plot, device = "pdf", units = "in", width = 5, height = 5, dpi = 600)





# ================================
# Libraries
# ================================
library(tidyverse)
library(ggrepel)
library(viridis)

# ================================
# Step 1: Load and Prepare Data
# ================================
data <- read_csv2("/g/zimmermann/Members/Nikita/GC-MSMS_application/Buffers_results/_calculated_concentrations_ACN.csv")

data <- data %>%
  mutate(
    Buffer = str_split(data_filename, "_", simplify = TRUE)[, 1],
    Replicate = str_split(data_filename, "_", simplify = TRUE)[, 2],
    SampleType = str_split(data_filename, "_", simplify = TRUE)[, 3]
  )

filtered_data <- data %>%
  filter(SampleType == "ACN", Buffer %in% c("FreshlyFrozen", "InvitekFrozen"))

# ================================
# Step 2: Summary + Stats
# ================================

stat_data <- filtered_data %>%
  group_by(name_original, Buffer) %>%
  summarize(
    Mean_Concentration = mean(concentration_calculated_uM, na.rm = TRUE),
    SD_Concentration = sd(concentration_calculated_uM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Buffer, values_from = c(Mean_Concentration, SD_Concentration))


p_values <- filtered_data %>%
  group_by(name_original) %>%
  summarize(
    p_value = tryCatch({
      v1 <- concentration_calculated_uM[Buffer == "FreshlyFrozen"]
      v2 <- concentration_calculated_uM[Buffer == "InvitekFrozen"]
      if (length(v1) > 1 && length(v2) > 1) {
        wilcox.test(v1, v2, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_value_fdr = p.adjust(p_value, method = "fdr"))

wide_concentrations <- stat_data %>%
  left_join(p_values, by = "name_original") %>%
  mutate(
    Significant = !is.na(p_value) & p_value < 0.031,
    Significant_FDR = !is.na(p_value_fdr) & p_value_fdr < 0.1
  )

# ================================
# Step 3: Transformations
# ================================
pseudo <- 1e-4

wide_concentrations <- wide_concentrations %>%
  mutate(
    Residual = Mean_Concentration_FreshlyFrozen - Mean_Concentration_InvitekFrozen,
    log10_Fresh = log10(ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen)),
    log10_Invitek = log10(ifelse(Mean_Concentration_InvitekFrozen == 0, pseudo, Mean_Concentration_InvitekFrozen)),
    log10_Residual = log10(abs(Residual) + pseudo) * sign(Residual)
  )

# ================================
# Step 4: Add Group Colors
# ================================
group_df <- read.csv2('/g/zimmermann/Members/Nikita/For_sunburst_4.csv')

group_mapping <- group_df %>%
  mutate(labels = as.character(labels), parents = as.character(parents)) %>%
  dplyr::select(labels, parents) %>%
  rename(name_original = labels, group = parents)

group_color_map <- c(
  "Benzenoids" = "darkviolet",
  "Organic acids" = "deepskyblue3",
  "Lipids" = "orange",           
  "Amino acids" = "chartreuse3",       
  "Amino acids derivatives" = "red",  
  "Imidazopyrimidines" = "coral4", 
  "Pyrimidines" = "bisque2",
  "Furans" = "chartreuse3",
  "Pyridines" = "darkolivegreen", 
  "Flavonoids, isoflavonoids" = "darkgrey",
  "Phenylpropanoids" = "deepskyblue3",
  "Alkaloids" = "coral",
  "Furanoid lignans" = "aquamarine4",
  "Organooxygens" = "darkturquoise",
  "Lignans" ="bisque1",
  "Organonitrogens" = "darksalmon"
)

wide_concentrations <- wide_concentrations %>%
  left_join(group_mapping, by = "name_original") %>%
  mutate(group_color = ifelse(is.na(group_color_map[group]), "black", group_color_map[group]))

# ================================
# Step 5: Manual Labels Setup
# ================================
manually_label <- c("Syringic acid", "Isovaleric acid", "Tyrosol", "Isoleucine", "Alanine", "Glutamic acid", "Cystine")

label_data <- wide_concentrations %>%
  filter(name_original %in% manually_label) %>%
  mutate(
    nudge_x = case_when(
      name_original == "Tyrosol" ~ 0.2,
      name_original == "Syringic acid" ~ 0.2,
      name_original == "Cystine" ~ 0.65,
      name_original == "Isovaleric acid" ~ -0.3,
      name_original == "Isoleucine" ~ 0.65,
      name_original == "Alanine" ~ 0.65,
      name_original == "Glutamic acid" ~ 0.65,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      name_original == "Tyrosol" ~ -0.5,
      name_original == "Syringic acid" ~ 0.2,
      name_original == "Cystine" ~ -1.65,
      name_original == "Isovaleric acid" ~ -0.6,
      name_original == "Isoleucine" ~ -1.65,
      name_original == "Alanine" ~ -1.65,
      name_original == "Glutamic acid" ~ -1.25,
      TRUE ~ 0
    )
  )

# ================================
# Step 6: Plot
# ================================
min_val <- min(c(wide_concentrations$log10_Fresh, wide_concentrations$log10_Invitek), na.rm = TRUE) - 0.5
max_val <- max(c(wide_concentrations$log10_Fresh, wide_concentrations$log10_Invitek), na.rm = TRUE) + 0.5

plot <- ggplot(wide_concentrations, aes(x = log10_Fresh, y = log10_Invitek)) +
  geom_point(aes(fill = group_color), shape = 21, size = 3, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  
  geom_errorbar(data = label_data,
                aes(
                  ymin = log10_Invitek - (SD_Concentration_InvitekFrozen / (ifelse(Mean_Concentration_InvitekFrozen == 0, pseudo, Mean_Concentration_InvitekFrozen) * log(10))),
                  ymax = log10_Invitek + (SD_Concentration_InvitekFrozen / (ifelse(Mean_Concentration_InvitekFrozen == 0, pseudo, Mean_Concentration_InvitekFrozen) * log(10)))
                ),
                width = 0.1, color = "grey30") +
  
  geom_errorbarh(data = label_data,
                 aes(
                   xmin = log10_Fresh - (SD_Concentration_FreshlyFrozen / (ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen) * log(10))),
                   xmax = log10_Fresh + (SD_Concentration_FreshlyFrozen / (ifelse(Mean_Concentration_FreshlyFrozen == 0, pseudo, Mean_Concentration_FreshlyFrozen) * log(10)))
                 ),
                 height = 0.1, color = "grey30") +
  
  geom_text_repel(
    data = label_data,
    aes(label = name_original, color = group_color),
    size = 6,
    box.padding = 0.1,
    point.padding = 0.1,
    force = 825,
    force_pull = 4.2,
    segment.size = 0.5,
    segment.color = "black",
    segment.curvature = 0,
    segment.linetype = "solid",
    min.segment.length = 0,
    nudge_x = label_data$nudge_x,
    nudge_y = label_data$nudge_y
  ) +
  
  labs(
    x = "Mean concentration (Freshly Frozen, log10 scale)",
    y = "Mean concentration (Invitek Frozen, log10 scale)"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    breaks = seq(0, max_val, 1),
    labels = parse(text = paste0("10^", seq(0, max_val, 1)))
  ) +
  scale_y_continuous(
    breaks = seq(0, max_val, 1),
    labels = parse(text = paste0("10^", seq(0, max_val, 1)))
  ) +
  coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 0.99, vjust = 0.99),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ================================
# Step 7: Save Plot
# ================================
ggsave("/g/zimmermann/Members/Nikita/GC-MSMS_application/Figure Panels/Illustrator/Invitek_scatter.pdf", 
       plot = plot, device = "pdf", units = "in", width = 5, height = 5, dpi = 600)


