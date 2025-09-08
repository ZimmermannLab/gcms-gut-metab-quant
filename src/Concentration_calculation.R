# GC-MS/MS Data Processing Pipeline
# Author: Nikita Denisov
# Description: Computes calibration curves and concentrations for labelled and unlabelled metabolites.

# Libraries ----
library(tidyverse)
library(progress)

# Source helper functions
# Place utils.R in the same folder or adjust path
source("utils.R")

#* 01: Data import ----
input_file_path <- "data/input_file.tsv"   # <-- adjust as needed
df <- f_import_from_filepath(input_file_path)

#* 02: Compute Area Ratios ----

## Labelled metabolites
labelled_metabolites_vec <- df %>%
  filter(type == "labelled") %>%
  pull(name_clean) %>%
  unique() %>%
  sort()

labelled_metabolites_df <- df %>%
  filter(name_clean %in% labelled_metabolites_vec)

processed_labelled_metabolites_df <- f_compute_area_ratios_labelledMetabolites(
  labelled_metabolites_df
)

## Unlabelled metabolites (map to internal standards)
metabolite_to_IS_map <- list(
  hexanoic_acid         = "heptanoic_acid",
  l_valine              = "n_methylnicotin_amide",
  malonic_acid          = "dimethylmalonic",
  indole_propionic_acid = c("ectoine","indole_6_carboxaldehyde","indole_3_acetonitrile","propylparaben",
                            "sumikis_acid","tyrosol","harman","6_hydroxynicotinic_acid",
                            "ureidopropionic_acid","quinolinic_acid","3_indoleacetic_acid",
                            "tryptamine","palmitoleic_acid","3_methyluric_acid","kynurenine",
                            "indole_3_carbinol","indole_3_lactic_acid","palmitoyl_ethanolamide",
                            "3_indoleacrylic_acid","phenylacetil_glutaric_acid","stearoylethanolamide",
                            "piperine","indole_pyruvic_acid","3_hydroxy_dl_kynurenine","spermine",
                            "daidzen"),
  succinic_acid         = c("dimethylsuccinic","2_aminoheptanoic_acid","2_hydroxyoctanoic_acid","5_hydroxyhexanoic_acid"),
  gaba                  = "putrescine",
  fumaric_acid          = c("2_aminoheptanoic_acid","2_hydroxyoctanoic_acid","5_hydroxyhexanoic_acid"),
  l_proline             = "methylcysteine",
  valeric_acid          = "phenylvaleric_acid",
  adipic_acid           = c("levulinic_acid","3_hydroxypalmitic_acid","2_hydroxyhippuric_acid",
                            "4_hydroxyhippuric_acid","2_hydroxystearic_acid","hexadecanedioic_acid",
                            "cholestan_3_ol"),
  l_threonine           = "agmatine",
  ketoglutaric_acid     = "5_keto_d_gluconate",
  l_phenylalanine       = "acetylspermidine",
  malic_acid            = "citramalic_acid",
  phthalic_acid         = c("hydroxyphenylacetic_acid","dihydroxybenzoic_acid"),
  aconitic_acid         = "erythronic_acid",
  uric_acid             = c("arachidonic_acid","cholestenone","enterolactone","stigmasterol","nutriacholic_acid"),
  orotic_acid           = c("syringic_acid","dihidroferulic_acid"),
  vanillylmandelic_acid = c("tricarballylic_acid","undecanedioic_acid"),
  l_lysine              = "cinnamoyglycine"
)

all_processed_unlabelled_df <- purrr::map_dfr(
  names(metabolite_to_IS_map),
  function(IS) {
    group <- metabolite_to_IS_map[[IS]]
    group_df <- df %>%
      filter(name_clean %in% group | (name_clean == IS & type == "labelled"))
    f_compute_area_ratios_unlabelledMetabolites(group_df, IS)
  }
)

#* 03: Generate calibration curves ----
output_folder <- "results/calibration_plots/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

f_generate_calibration_curves(
  ratios_calibration_curve_df = processed_labelled_metabolites_df,
  save_plot_folder            = output_folder,
  coefficients_file_name      = file.path(output_folder, "coefficients_labelled.tsv"),
  points_to_exclude           = list()
)

f_generate_calibration_curves(
  ratios_calibration_curve_df = all_processed_unlabelled_df,
  save_plot_folder            = output_folder,
  coefficients_file_name      = file.path(output_folder, "coefficients_unlabelled_multiple_IS.tsv")
)

#* 04: Combine coefficient files ----
coef_labelled_df   <- read_tsv(file.path(output_folder, "coefficients_labelled.tsv"))
coef_unlabelled_df <- read_tsv(file.path(output_folder, "coefficients_unlabelled_multiple_IS.tsv"))
coef_combined_df   <- bind_rows(coef_labelled_df, coef_unlabelled_df)

if (any(duplicated(coef_combined_df$name_clean))) {
  warning("Duplicated metabolites detected:")
  print(coef_combined_df %>%
          filter(duplicated(name_clean) | duplicated(name_clean, fromLast = TRUE)))
}
stopifnot(!any(duplicated(coef_combined_df$name_clean)))

write_tsv(coef_combined_df, file.path(output_folder, "coefficients_combined.tsv"))

#* 05: Process unknown metabolites ----
df_unknown_metabolites <- f_import_from_filepath("data/results_unknown.tsv")
coefficients_file_path <- file.path(output_folder, "coefficients_combined.tsv")

## Labelled
labelled_metabolites_vec <- df_unknown_metabolites %>%
  filter(type == "labelled") %>%
  pull(name_clean) %>%
  unique()

processed_labelled_unknown_metabolites_df <- f_compute_area_ratios_labelledMetabolites(
  df_unknown_metabolites %>% filter(name_clean %in% labelled_metabolites_vec)
)

## Unlabelled
processed_unlabelled_unknown_metabolites_df <- purrr::map_dfr(
  names(metabolite_to_IS_map),
  function(IS) {
    group <- metabolite_to_IS_map[[IS]]
    group_df <- df_unknown_metabolites %>%
      filter(name_clean %in% group | (name_clean == IS & type == "labelled"))
    f_compute_area_ratios_unlabelledMetabolites(group_df, IS)
  }
)

## Combine all
area_ratios_unknown_metabolites_combined <- bind_rows(
  processed_labelled_unknown_metabolites_df,
  processed_unlabelled_unknown_metabolites_df
)

#* 06: Calculate concentrations ----
calculated_concentrations_unknown_metabolites_df <- f_compute_sample_concentration(
  unknown_metabolites_df   = area_ratios_unknown_metabolites_combined,
  coefficients_file_path   = coefficients_file_path
)

write_tsv(
  calculated_concentrations_unknown_metabolites_df,
  "results/calculated_concentrations.tsv"
)
