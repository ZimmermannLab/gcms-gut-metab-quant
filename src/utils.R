#' ---
#' title: "GC–MS utilities"
#' author: "Fabian Springer"
#' manuscript: "Denisov et al., Development of a GC-MS/MS method to quantify 120 gut microbiota-derived metabolites "
#' date: "2025-09-08"
#' license: "MIT"
#' description: >
#'   Helper functions to import instrument exports, compute area ratios
#'   for labelled/unlabelled compounds, fit calibration curves, infer
#'   concentrations for unknowns, plot QC figures, and export results.
#' notes: >
#'   - This script assumes a specific tab-delimited export format
#'     (see f_import_from_filepath()).
#' ---

require(parallel)
require(tidyverse)
require(janitor)
require(progress)
require(openxlsx)

#' Theme for ggplot figures
#'
#' @param font.size Base font size.
#' @param line.size Base line width for lines and borders.
#' @param panel.grid Grid visibility: 'none', 'major', 'full', 'major_x', 'major_y'.
#'
#' @return A ggplot2 theme object.
#' @examples
#' # ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_embl()
theme_embl <- function(font.size = 12, line.size=0.45, panel.grid='none') {
  # Taken from ggembl (internal EMBL repo), lightly adapted for this project.
  if (!panel.grid %in% c('none', 'major', 'full', 'major_x', 'major_y')){
    warning("Unrecognized panel.grid info... Reverting to 'none'")
    panel.grid <- 'none'
  }
  thm <- theme_bw() +
    theme(line = element_line(colour = "black",
                              lineend = "square",
                              linetype = "solid",
                              size=line.size),
          rect = element_rect(fill = NA,
                              colour = "black",
                              linetype = "solid",
                              size=line.size),
          text = element_text(colour = "black",
                              face = "plain",
                              size = font.size,
                              vjust = 0.5,
                              hjust = 0.5,
                              lineheight = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA, colour='black', size=line.size),
          strip.background = element_blank(),
          legend.key = element_blank(),
          title = element_text(size = rel(1)),
          plot.title = element_text(size = rel(1.2), face = "bold"),
          strip.text = element_text(colour='black'),
          axis.ticks.length = unit(0.3, "lines"),
          axis.text = element_text(size=rel(0.8), colour='black'),
          axis.ticks = element_line(colour='black'))
  if (panel.grid=='none'){
    thm <- thm + theme(panel.grid = element_blank())
  } else if (stringr::str_detect(panel.grid,'major')){
    thm <- thm + theme(panel.grid = element_blank())
    if (panel.grid=='major') {
      thm <- thm + theme(panel.grid.major = element_line(colour='lightgrey',
                                                         size=0.6*line.size))
    } else if (panel.grid=='major_x'){
      thm <- thm + theme(panel.grid.major.x = element_line(colour='lightgrey',
                                                           size=0.6*line.size))
    } else if (panel.grid =='major_y'){
      thm <- thm + theme(panel.grid.major.y = element_line(colour='lightgrey',
                                                           size=0.6*line.size))
    }
  } else if (panel.grid=='full'){
    thm <- thm + theme(panel.grid.minor = element_line(colour='lightgrey',
                                                       size=0.3*line.size),
                       panel.grid.major = element_line(colour='lightgrey',
                                                       size=0.6*line.size))
  }
  return(thm)
}

# Project default theme with bolder axes/legend for figures in the manuscript
theme_fs <- theme_embl(font.size = 12, line.size = 0.7, panel.grid = "major") +
  theme(
    axis.title = element_text(face = "bold", size = 12), # Bold axis titles
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

#' Import instrument export file (tab-delimited)
#'
#' Parses the instrument's standard export into a tidy tibble, detects
#' labelled vs unlabelled entries, cleans names, and extracts known
#' concentrations from file names.
#'
#' @param input_file_path Path to the raw export (.txt/.tsv) file.
#'
#' @return Tibble with columns:
#'   name_original, name_clean, data_filename, type,
#'   area, (other numeric QC columns), concentration_known (numeric)
#'
#' @details
#' Expects blocks delimited by rows where column 1 == "Name".
#' The function temporarily trims to the first 7 columns.
#'
#' @examples
#' # df <- f_import_from_filepath("data/raw_export.txt")
f_import_from_filepath <- function(input_file_path) {
  # Basic file check and read
  if (file.exists(input_file_path)) {
    data_raw <- read.delim(input_file_path, header = F, sep = "\t", comment.char = "@") %>%
      as.data.frame()
  } else {
    stop("The file does not exist")
  }
  
  # Manually remove column beyond the first 7 (temporary workaround for current export)
  data_raw <- data_raw[,1:7]
  if (ncol(data_raw) != 7) {
    stop("The number of columns is not 7")
  }

  # Remove empty rows (rows where all entries are "")
  dim(data_raw) # keeps a brief console output (intentional)
  data_raw <- data_raw[!apply(data_raw == "", 1, all), ]

  # Identify the start of each metabolite block via "Name" in column 1
  row_indices_NAME <- which(str_detect(data_raw[, 1], "Name"))
  if (length(row_indices_NAME) == 0) {
    stop("No 'Name' entries found in column 1 of the provided file - check the format (or ask Fabi)")
  }

  # Accumulate parsed blocks here
  metabolite_df <- data.frame(matrix(ncol = 7, nrow = 0))

  # Iterate over blocks and collect rows
  pb <- progress_bar$new(total = length(row_indices_NAME))
  message("Processing ", length(row_indices_NAME), " metabolites")
  for (i in seq(1, length(row_indices_NAME))) {
    entry_start <- row_indices_NAME[i] + 2 # data start (2 rows after "Name")

    # End of block: 2 rows before the next "Name" or end of file
    if (i == length(row_indices_NAME)) {
      entry_end <- nrow(data_raw)
    } else {
      entry_end <- row_indices_NAME[i + 1] - 2
    }

    # Extract fields for this metabolite and tag with its name
    block_entries <- data_raw[entry_start:entry_end, ]
    metabolite_name <- data_raw[entry_start - 2, 2]
    block_entries$V1 <- metabolite_name

    metabolite_df <- rbind(metabolite_df, block_entries)
    pb$tick()
  }

  # Assign column names based on the header row following the first "Name"
  col_names <- c("Name", as.character(data_raw[row_indices_NAME[1] + 1, 2:7]))
  colnames(metabolite_df) <- col_names

  # Tidy: to tibble, clean names
  metabolite_df <- metabolite_df %>%
    as_tibble() %>%
    janitor::clean_names()

  # Convert numeric-like columns (with commas as decimal separators) to numeric
  metabolite_df <- metabolite_df %>%
    mutate(across(area:percent_rsd_area, ~ str_replace_all(.x, ",", ".") %>% as.numeric()))

  # Detect labelled/unlabelled, standardize names, and extract known concentration
  mod_metabolite_df <- metabolite_df %>%
    mutate(
      # Labeled molecules have ", " in the original name
      type = ifelse(str_detect(name, ", "), "labelled", "unlabelled"),
      # Cleaned lower-case name without trailing label info (after comma)
      name_clean = str_to_lower(str_remove(name, ", .*")),
      # Replace spaces/dashes with underscores for safe identifiers
      name_clean = str_replace_all(name_clean, pattern = " |\\-", replacement = "_")
    ) %>% 
    mutate(
      # Known concentration parsed from data_filename (e.g., "0_6859" -> "0.6859")
      concentration_known = case_when(
        type %in% c("labelled","unlabelled") ~ str_replace(data_filename,pattern = "_",replacement = "."),
        TRUE ~ NA_character_
      ),
      concentration_known = as.numeric(concentration_known)
    ) %>% 
    glimpse() %>%  # intentional console glimpse for QC
    dplyr::rename(name_original = name) %>% 
    dplyr::relocate(name_original, name_clean, data_filename, type)

  return(mod_metabolite_df)
}

#' Compute area ratios using the labelled counterpart as reference
#'
#' For each (name_clean, data_filename) group, divides the area of each row
#' by the area of the corresponding labelled entry.
#'
#' @param labelled_metabolites_df Tibble from importer, containing both labelled and unlabelled rows.
#'
#' @return Input tibble with an added `area_ratio` column and `reference_metabolite = "labelled"`.
#' @examples
#' # ratios <- f_compute_area_ratios_labelledMetabolites(df)
f_compute_area_ratios_labelledMetabolites <- function(labelled_metabolites_df) {
  N_labelled_metabs <- unique(labelled_metabolites_df$name_clean) %>% length()
  message("Number of labelled metabolites: ",N_labelled_metabs)

  # Identify presence/absence of labelled/unlabelled counterparts per data_filename
  metab_labelled <- labelled_metabolites_df %>%
    filter(type == "labelled") %>%
    dplyr::select(name_clean,data_filename) %>%
    distinct() %>% 
    arrange(data_filename,name_clean)

  metab_unlabelled <- labelled_metabolites_df %>%
    filter(type == "unlabelled") %>%
    dplyr::select(name_clean,data_filename) %>%
    distinct() %>% 
    arrange(data_filename,name_clean)
  
  m_labelled_vec <- metab_labelled$name_clean %>% unique()
  m_unlabelle_vec <- metab_unlabelled$name_clean %>% unique()

  # Report missing counterparts
  missing_metabolites_unlabelled <- c()
  for(m in m_labelled_vec){
    if(!(m %in% m_unlabelle_vec)){
      missing_metabolites_unlabelled <- c(missing_metabolites_unlabelled,m)
    }
  }

  missing_metabolites_labelled <- c()
  for(m in m_unlabelle_vec){
    if(!(m %in% m_labelled_vec)){
      missing_metabolites_labelled <- c(missing_metabolites_labelled,m)
    }
  }
  
  if (length(missing_metabolites_unlabelled) > 0) {
    message(
      "Attention: The following LABELLED metabolites do not have an UNLABELLED counterpart:\n - ",
      paste0(missing_metabolites_unlabelled, collapse = "\n - ")
    )
  }
  if (length(missing_metabolites_labelled) > 0) {
    message(
      "Attention: The following UNLABELLED metabolites do not have a LABELLED counterpart:\n",
      paste0(missing_metabolites_labelled, collapse = "\n")
    )
  }
  
  # Compute unlabelled/labelled ratios within each (name_clean, data_filename)
  res_df <- labelled_metabolites_df %>%
    inner_join(.,inner_join(metab_labelled,metab_unlabelled)) %>% 
    group_by(name_clean, data_filename) %>%
    mutate(area_ratio = area / area[type == "labelled"]) %>%  # assumes exactly one labelled per group
    dplyr::relocate(area_ratio) %>%
    ungroup() %>%
    mutate(reference_metabolite = "labelled") %>% 
    dplyr::relocate(name_original, name_clean, type, data_filename, concentration_known, area, area_ratio,reference_metabolite)
  
  return(res_df)
}

#' Compute area ratios for unlabelled metabolites vs a chosen reference metabolite
#'
#' Keeps the labelled entry of the reference metabolite and divides each metabolite's
#' area by the reference area when present for the same `data_filename`.
#'
#' @param unlabelled_metabolites_df Tibble from importer containing both types.
#' @param reference_molecule_name_clean Clean name (lowercase, underscores) of the reference metabolite.
#'
#' @return Tibble with `area_ratio` and `reference_metabolite` set to the provided reference.
#' @examples
#' # ratios <- f_compute_area_ratios_unlabelledMetabolites(df, "glucose")
f_compute_area_ratios_unlabelledMetabolites <- function(unlabelled_metabolites_df, reference_molecule_name_clean){
  message("Computing area ratios of unlabelled metabolites based on reference molecule: ",reference_molecule_name_clean)
  N_unlab_metabs <- unlabelled_metabolites_df %>% filter(name_clean != reference_molecule_name_clean) %>% pull(name_clean) %>% unique() %>% length()
  message("Number of unlabelled metabolites: ",N_unlab_metabs)

  # Validate the presence of the labelled reference
  labelled_ref_molecule <- unlabelled_metabolites_df %>% filter(type == "labelled") %>% pull(name_clean) %>% unique()
  if(!(reference_molecule_name_clean %in% labelled_ref_molecule)){
    stop("The reference molecule in its labelled form is not present in the provided dataframe")
  }

  # Warn if some metabolites do not share the same data_filename set as the reference
  reference_data_filenames <- unlabelled_metabolites_df %>% filter(name_clean == reference_molecule_name_clean & type == "labelled") %>% pull(data_filename) %>% unique()
  unmatched_filenames_df <- unlabelled_metabolites_df %>%
    filter(name_clean != reference_molecule_name_clean) %>%
    filter(!(data_filename %in% reference_data_filenames)) %>% 
    dplyr::select(name_clean,data_filename) %>% distinct()

  if (nrow(unmatched_filenames_df) > 0) {
    message("Following metabolites have different data_filenames than the selected reference:")
    for(m in unique(unmatched_filenames_df$name_clean)){
      if(!(m %in% labelled_ref_molecule)){
        message(paste0(m,": ",paste0(unmatched_filenames_df %>% filter(name_clean == m) %>% pull(data_filename) %>% unique(),collapse = "; ")))
      }
    }        
  }

  # Compute ratios vs reference (keep labelled ref row; drop ref at the end)
  res_df <-
    unlabelled_metabolites_df %>%
    filter(name_clean != reference_molecule_name_clean | type == "labelled") %>% 
    group_by(data_filename) %>%     # data_filename acts as concentration identifier
    mutate(area_ratio = if (any(name_clean == reference_molecule_name_clean)) area / area[name_clean == reference_molecule_name_clean] else NA_real_) %>%    
    dplyr::relocate(area_ratio) %>%
    ungroup() %>%
    mutate(reference_metabolite = reference_molecule_name_clean) %>%
    dplyr::relocate(name_original, name_clean, type, data_filename, concentration_known, area, area_ratio, reference_metabolite) %>%
    filter(name_clean != reference_molecule_name_clean) # remove the reference molecule
  
  res_df %>% filter(data_filename == "0_6859")  # intentional side-effect line as in original
  
  return(res_df)
}

#' Plot a calibration curve and fit a linear model (log10-log10)
#'
#' Takes a single-molecule tibble (multiple concentrations) and fits a linear
#' model in log10 space, returning the ggplot and model coefficients.
#'
#' @param molecule_df Tibble for one molecule (mixture of include/exclude points).
#' @param measured_values (unused placeholder for future extension)
#' @param points_to_exclude_vec Numeric vector of concentrations (rounded via signif(.,2)) to exclude.
#'
#' @return List with elements `plot` (ggplot) and `model_df` (slope, intercept, r_squared).
#' @examples
#' # res <- f_plot_calibration_curve(df_mol, points_to_exclude_vec = c(0.1))
f_plot_calibration_curve <- function(molecule_df,measured_values=NULL,points_to_exclude_vec=NULL){
  # Normalize exclude vector
  if(is.null(points_to_exclude_vec)){
    points_to_exclude_vec <- c()
  }

  # Add log10 columns and mark points to include/exclude
  molecule_df <- molecule_df %>% mutate(
    l10_conc = log10(concentration_known),
    l10_area_ratio = log10(area_ratio),
    include_point = ifelse(signif(concentration_known,2) %in% points_to_exclude_vec, "exclude", "include")
    )
  
  molecule_name <- unique(molecule_df$name_original)
  stopifnot("Provided dataframe contains more than one molecule!" = length(molecule_name) == 1)
  
  # Basic ranges for annotations/axes
  max_y <- max(molecule_df$l10_area_ratio,na.rm = T)
  min_y <- min(molecule_df$l10_area_ratio,na.rm = T)
  min_x <- min(molecule_df$l10_conc,na.rm = T)
  max_x <- max(molecule_df$l10_conc,na.rm = T)
  
  # Fit a simple linear model in log10 space using included points
  mod <- lm(l10_area_ratio ~ l10_conc,molecule_df %>% filter(include_point == "include"))
  coefficients <- coef(mod)
  intercept <- signif(coefficients[1], 2)
  slope <- signif(coefficients[2], 2)
  r_squared <- signif(summary(mod)$r.squared,2)
  
  # Prediction line (for visual guide)
  pred_df <- tibble(l10_conc = seq(min_x,max_x,length.out = 10))
  prediction <- predict(mod,pred_df)
  pred_df$l10_area_ratio <- prediction
  
  # Build plot
  xBreaks <- molecule_df$l10_conc
  xLabs <- 10^xBreaks
  p <-    
    molecule_df %>% 
    ggplot(aes(y = l10_area_ratio, x = l10_conc)) +
    geom_line(data = pred_df) +
    # included points
    geom_point(
      data = molecule_df %>% filter(!is.na(l10_area_ratio),include_point == "include"),
      pch = 21, stroke = 1, fill = "lightgrey", size = 2.5) +
    # excluded points
    geom_point(  
      data = molecule_df %>% filter(!is.na(l10_area_ratio),include_point == "exclude"),
      pch = 21, fill = "white", size = 2.5)+
    ylab("Area ratio [log10]") +
    xlab("Concentration [uM]") +
    ggtitle(str_to_title(molecule_name)) +
    scale_x_continuous(breaks = xBreaks, labels = signif(xLabs, 2)) +
    theme_fs
  
  # Add R^2 annotation
  text <- paste("R^2 =", r_squared)
  p <- p +
    geom_text(
      data = data.frame(x = -Inf, y = max_y, label = text),
      x = -Inf, y = max_y, label = text, hjust = -0.2, vjust = 0.5)
  
  # Return model coefficients
  mod_df <- tibble(slope = slope, intercept = intercept, r_squared = r_squared)
  res_list <- list(plot = p, model_df = mod_df)
  return(res_list)
}

#' Generate calibration curves for all molecules and save figures + coefficients
#'
#' Iterates over molecules in `ratios_calibration_curve_df`, fits curves via
#' f_plot_calibration_curve(), saves per-molecule PDFs, and writes a TSV of
#' coefficients (slope, intercept, R^2).
#'
#' @param ratios_calibration_curve_df Tibble with area ratios for unlabelled molecules across concentrations.
#' @param save_plot_folder Output folder for per-molecule plots (PDF).
#' @param coefficients_file_name Optional path for TSV; defaults to <save_plot_folder>/_model_coefficients.tsv.
#' @param points_to_exclude Named list: molecule_name_clean -> numeric vector of concentrations to exclude.
#'
#' @return Invisibly returns the path to the coefficients file (side-effect: writes files).
#' @examples
#' # f_generate_calibration_curves(ratios_df, "results/calib_plots")
f_generate_calibration_curves <- function(ratios_calibration_curve_df, save_plot_folder, coefficients_file_name=NULL,points_to_exclude = NULL) {
  c <- 1 # dummy (kept intentionally)

  if(!("area_ratio" %in% colnames(ratios_calibration_curve_df))){
    stop("Provided dataframe does not contain the column 'area_ratio'")
  }
  
  out_file_name <- if (is.null(coefficients_file_name)) {
    file.path(save_plot_folder,"_model_coefficients.tsv")
  } else {
    coefficients_file_name
  }
  
  if(!(dir.exists(save_plot_folder))){dir.create(save_plot_folder,recursive = T)}
  
  all_molecules <- unique(ratios_calibration_curve_df$name_clean)
  message("Processing ",length(all_molecules)," molecules")

  pb <- progress::progress_bar$new(total = length(all_molecules))
  model_coefficients_df <- tibble()
  mol <- all_molecules[2] # pre-set for debug (kept intentionally)

  for (mol in all_molecules) {
    message("\n",mol)

    # Keep only unlabelled rows for this molecule
    c_df <- ratios_calibration_curve_df %>%
      filter(name_clean == mol) %>%
      filter(type == "unlabelled")

    if (nrow(c_df) == 0 | !any(!is.na(c_df$area_ratio))) {
      message("No non-NA area-ratios for ", mol, " - skipping")
      pb$tick()
      next
    }

    # Per-molecule exclusions
    if (!is.null(points_to_exclude) & mol %in% names(points_to_exclude)) {
      points_to_exclude_vec <- points_to_exclude[[mol]]
    }else{
      points_to_exclude_vec <- c()
    }

    # Fit and collect coefficients
    c_res <- f_plot_calibration_curve(molecule_df = c_df,points_to_exclude_vec = points_to_exclude_vec)
    c_model_df <- c_res$model_df

    model_coefficients_df <- bind_rows(model_coefficients_df, c_model_df %>%
      mutate(
        name_clean = mol,
        name_original = unique(c_df$name_original)
      ))

    # Save plot
    c_plot <- c_res$plot + theme(aspect.ratio = 1)
    ggsave(c_plot, filename = file.path(save_plot_folder, paste0(mol, ".pdf")), width = 5, height = 5)
    pb$tick()
  }

  # Write coefficients table
  write_tsv(model_coefficients_df %>% dplyr::relocate(name_original,name_clean), file = out_file_name)
  
  message("-------------- Finished --------------")
  message("Plots_saved_to ", save_plot_folder)
  message("Model coefficients saved to file: ", out_file_name)
}

#' Compute concentrations for unknown samples using calibration coefficients
#'
#' Joins unknown sample area ratios to fitted model coefficients and back-solves
#' the linear model (in log10 space) to obtain concentration estimates.
#'
#' @param unknown_metabolites_df Tibble with columns including area_ratio (unlabelled).
#' @param coefficients_file_path Path to TSV produced by f_generate_calibration_curves().
#'
#' @return Tibble with calculated concentrations (linear and log10 scales) and model diagnostics.
#' @examples
#' # conc_df <- f_compute_sample_concentration(unknown_df, "results/_model_coefficients.tsv")
f_compute_sample_concentration <- function(unknown_metabolites_df, coefficients_file_path) {
  if(!file.exists(coefficients_file_path)){
    stop("The provided coefficients file does not exist")
  }else{
    coefficients_df <- read_tsv(coefficients_file_path) %>% dplyr::select(-name_original)
  }

  molecule_vec <- unique(unknown_metabolites_df$name_clean)
  message("Processing ",length(molecule_vec)," molecules")

  # Only consider molecules present in the coefficients file
  molecule_vec <- molecule_vec[molecule_vec %in% coefficients_df$name_clean]
  message(length(molecule_vec)," molecules present in the coefficients file - skipping the others")

  calculated_concentrations_df <- tibble()
  mol <- molecule_vec[2] # pre-set for debug (kept intentionally)

  for (mol in molecule_vec) {
    message(mol)

    # Keep unlabelled rows for this molecule
    c_df <- unknown_metabolites_df %>%
      filter(name_clean == mol) %>%
      filter(type == "unlabelled")

    if (nrow(c_df) == 0 | !any(!is.na(c_df$area_ratio))) {
      message("No non-NA area-ratios for ", mol, " - skipping")
      next
    }
    
    # Join coefficients and compute log10(area_ratio)
    c_df <- suppressMessages(c_df %>%
      inner_join(., coefficients_df) %>% 
      mutate(l10_area_ratio = log10(area_ratio)))

    # Back-calculate concentration in log10 space; revert to linear
    c_conc_df <- c_df %>%
      dplyr::select(name_original, name_clean, data_filename, area, area_ratio, l10_area_ratio, slope:r_squared,reference_metabolite) %>%
      mutate(
        concentration_calculated_log10_uM = (l10_area_ratio - intercept) / slope,
        concentration_calculated_uM = 10^concentration_calculated_log10_uM
      ) %>%
      dplyr::relocate(concentration_calculated_uM, concentration_calculated_log10_uM, .after = area_ratio)
    
    calculated_concentrations_df <- bind_rows(calculated_concentrations_df,c_conc_df)    
  }
  return(calculated_concentrations_df %>% dplyr::select(name_original:concentration_calculated_uM,reference_metabolite))
}

#' Export a data frame to Excel (.xlsx)
#'
#' Writes a single-sheet workbook with auto-width columns and a header filter.
#'
#' @param dataframe Tibble/data.frame to export.
#' @param out_file_name Target path (".xlsx" added if missing).
#'
#' @return Invisibly returns the written path (side-effect: writes a file).
#' @examples
#' # f_export_to_excel(results_df, "results/unknown_concentrations")
f_export_to_excel <- function(dataframe,out_file_name){
  # Add .xlsx suffix if not present
  if(!str_detect(out_file_name,".xlsx$")){
    out_file_name <- paste0(out_file_name,".xlsx")
  }

  require(openxlsx)

  # Build workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Sheet1")
  writeData(wb, "Sheet1", dataframe)
  
  # Add column filter on the header row
  addFilter(wb, "Sheet1", cols = 1:ncol(dataframe),row = 1)

  # Auto column widths
  setColWidths(wb, "Sheet1", cols = 1:ncol(dataframe), widths = "auto")

  saveWorkbook(wb, file = out_file_name, overwrite = TRUE)
  message("Dataframe saved to: ",out_file_name)
}
