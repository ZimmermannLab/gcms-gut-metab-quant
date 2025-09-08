require(parallel)

theme_embl <- function(font.size = 12, line.size=0.45, panel.grid='none') {
  #* Theme for plotting, taken from ggembl: ----
  # https://git.embl.de/grp-zeller/ggembl/-/blob/master/R/themes.R
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
theme_fs <- theme_embl(font.size = 12, line.size = 0.7, panel.grid = "major") +
  theme(
    axis.title = element_text(face = "bold", size = 12), # Bold axis titles
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    # axis.line = element_line(size = 1.5), # Change line thickness of axes

    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


f_import_from_filepath <- function(input_file_path) {
  #* this function imports the standard output files provided by the XZY software ----
  # loop over every row and check whether "Name" is in column 1.
  # If so, this characterizes a metabolite block.
  # The matbolite block ends at the next Name index or at the end of the file
  # Generate a clean tibble with the metabolites and check for labelled/unlabelled status/ concentrations etc

  # Do a few checks
  if (file.exists(input_file_path)) {
    data_raw <- read.delim(input_file_path, header = F, sep = "\t", comment.char = "@") %>%
      as.data.frame()
  } else {
    stop("The file does not exist")
  }
  
  #just for now: remove Conc. ppm column manually (wont be in the real export according to Nikita)
  data_raw <- data_raw[,1:7]
  if (ncol(data_raw) != 7) {
    stop("The number of columns is not 7")
  }

  # # get index of empty rows
  # empty_rows <- which(apply(data_raw == "", 1, all))

  # Remove empty rows
  dim(data_raw)
  data_raw <- data_raw[!apply(data_raw == "", 1, all), ]


  #* Part1: Import the raw file and prepare the table ----
  row_indices_NAME <- which(str_detect(data_raw[, 1], "Name"))
  if (length(row_indices_NAME) == 0) {
    stop("No 'Name' entries found in column 1 of the provided file - check the format (or ask Fabi)")
  }

  metabolite_df <- data.frame(matrix(ncol = 7, nrow = 0))
  # extract and prepare the blocks: loop over row indices, extrac the rows, generate new DF with the respective fields
  pb <- progress_bar$new(total = length(row_indices_NAME))
  message("Processing ", length(row_indices_NAME), " metabolites")
  for (i in seq(1, length(row_indices_NAME))) {
    entry_start <- row_indices_NAME[i] + 2 # here starts the entry

    if (i == length(row_indices_NAME)) { # define the end of the entry: either 2rows before the next Name index or EOF
      entry_end <- nrow(data_raw)
    } else {
      entry_end <- row_indices_NAME[i + 1] - 2
    }

    # extract the fields
    block_entries <- data_raw[entry_start:entry_end, ]
    metabolite_name <- data_raw[entry_start - 2, 2]
    block_entries$V1 <- metabolite_name


    # Add to metabolite_df
    metabolite_df <- rbind(metabolite_df, block_entries)
    pb$tick()
  }

  # Extract the column names from the raw df and assign
  col_names <- c("Name", as.character(data_raw[row_indices_NAME[1] + 1, 2:7]))
  colnames(metabolite_df) <- col_names

  # Convert to tibble and make colnames clean
  metabolite_df <- metabolite_df %>%
    as_tibble() %>%
    janitor::clean_names()

  # Convert all , to . and convert to numeric
  metabolite_df <- metabolite_df %>%
    mutate(across(area:percent_rsd_area, ~ str_replace_all(.x, ",", ".") %>% as.numeric()))


  #* Part 2: Detect whether it is a labelled molecule and highlight accordingly ----
  # ! Important: The detection of labelled molecules works based on searech for "," followed by a space in the name
  mod_metabolite_df <- metabolite_df %>%
    mutate(
      type = ifelse(str_detect(name, ", "), "labelled", "unlabelled"),
      name_clean = str_to_lower(str_remove(name, ", .*")),
      name_clean = str_replace_all(name_clean, pattern = " |\\-", replacement = "_") # replace all spaces and dashes with underscores
    ) %>% 
    # For the labelled metabolites: extract concentration from the data_filename
    mutate(
      concentration_known = case_when(
        type %in% c("labelled","unlabelled") ~ str_replace(data_filename,pattern = "_",replacement = "."),
        TRUE ~ NA_character_
      ),
      concentration_known = as.numeric(concentration_known)
    ) %>% 
    glimpse() %>% 
    dplyr::rename(name_original = name) %>% 
    dplyr::relocate(name_original, name_clean, data_filename, type)

  return(mod_metabolite_df)
}


f_compute_area_ratios_labelledMetabolites <- function(labelled_metabolites_df) {
  #* Computes peak area ratios for metabolites compared to the labelled counterpart ----
  
  N_labelled_metabs <- unique(labelled_metabolites_df$name_clean) %>% length()
  message("Number of labelled metabolites: ",N_labelled_metabs)

  # Check if only labelled metabolites are provided
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
  # Check whether labelled metabolotes have a counterpart
  missing_metabolites_unlabelled <- c()
  for(m in m_labelled_vec){
    if(!(m %in% m_unlabelle_vec)){
      missing_metabolites_unlabelled <- c(missing_metabolites_unlabelled,m)
    }
  }
  # check whether unlabelled metabolites have a counterpart
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
  
  
  # generate the ratios of unlabelled/labelled
  res_df <- labelled_metabolites_df %>%
    inner_join(.,inner_join(metab_labelled,metab_unlabelled)) %>% 
    group_by(name_clean, data_filename) %>%
    mutate(area_ratio = area / area[type == "labelled"]) %>%
    dplyr::relocate(area_ratio) %>%
    ungroup() %>%
    mutate(reference_metabolite = "labelled") %>% 
    # re-arrange columns
    dplyr::relocate(name_original, name_clean, type, data_filename, concentration_known, area, area_ratio,reference_metabolite)
  
  return(res_df)
}

f_compute_area_ratios_unlabelledMetabolites <- function(unlabelled_metabolites_df, reference_molecule_name_clean){
  #* Computes area ratios of unlabelled metabolites based on reference molecule ----
  # Reference molecule name_clean must be present in the name_clean column

  message("Computing area ratios of unlabelled metabolites based on reference molecule: ",reference_molecule_name_clean)
  N_unlab_metabs <- unlabelled_metabolites_df %>% filter(name_clean != reference_molecule_name_clean) %>% pull(name_clean) %>% unique() %>% length()
  message("Number of unlabelled metabolites: ",N_unlab_metabs)

  # Check whether the reference molecule is present
  labelled_ref_molecule <- unlabelled_metabolites_df %>% filter(type == "labelled") %>% pull(name_clean) %>% unique()
  if(!(reference_molecule_name_clean %in% labelled_ref_molecule)){
    stop("The reference molecule in its labelled form is not present in the provided dataframe")
  }

  #@ Check if there are metabolites that have different data_filename than the selected reference
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

  # Compute the area ratios compared to the reference molecule  
  res_df <-
    unlabelled_metabolites_df %>%
    filter(name_clean != reference_molecule_name_clean | type == "labelled") %>% # make sure to keep only the labelled version of the reference molecule (if unlabelled was also provided)
    
    group_by(data_filename) %>% # only group by filename (which represents concentrations)     
     #! if the reference molecule is present with the same data_filename, compute the area ratio, if not, then return NA
    mutate(area_ratio = if (any(name_clean == reference_molecule_name_clean)) area / area[name_clean == reference_molecule_name_clean] else NA_real_) %>%    
    
    dplyr::relocate(area_ratio) %>%
    ungroup() %>%
    mutate(reference_metabolite = reference_molecule_name_clean) %>%
    # re-arrange columns
    dplyr::relocate(name_original, name_clean, type, data_filename, concentration_known, area, area_ratio, reference_metabolite) %>%
    filter(name_clean != reference_molecule_name_clean) # remove reference molecule
  
  res_df %>% filter(data_filename == "0_6859")
  
  return(res_df)

}

f_plot_calibration_curve <- function(molecule_df,measured_values=NULL,points_to_exclude_vec=NULL){
  ## Function to fit a linear model for the calibration curve
  # generates a plot of area ratio vs concentration
  # Returns ggplot and linear model
  #tmp_molDF
  if(is.null(points_to_exclude_vec)){
    points_to_exclude_vec <- c()
  }
  molecule_df <- molecule_df %>% mutate(
    l10_conc = log10(concentration_known),
    l10_area_ratio = log10(area_ratio),
    include_point = ifelse(signif(concentration_known,2) %in% points_to_exclude_vec, "exclude", "include")
    )
  
  molecule_name <- unique(molecule_df$name_original)
  stopifnot("Provided dataframe contains more than one molecule!" = length(molecule_name) == 1)
  
  # get maximum area_ratio
  max_y <- max(molecule_df$l10_area_ratio,na.rm = T)
  min_y <- min(molecule_df$l10_area_ratio,na.rm = T)
  
  # get min and max of concentration
  min_x <- min(molecule_df$l10_conc,na.rm = T)
  max_x <- max(molecule_df$l10_conc,na.rm = T)
  
  # fit a linear model
  #mod <- lm(concentration ~ area_ratio,molecule_df)
  mod <- lm(l10_area_ratio ~ l10_conc,molecule_df %>% filter(include_point == "include"))
  # Extract coefficients
  coefficients <- coef(mod)
  intercept <- signif(coefficients[1], 2)
  slope <- signif(coefficients[2], 2)
  
  # Extract R-squared value
  r_squared <- signif(summary(mod)$r.squared,2)
  
  # generate calibration curve for the provided concentrations
  pred_df <- tibble(l10_conc = seq(min_x,max_x,length.out = 10))
  prediction <- predict(mod,pred_df)
  pred_df$l10_area_ratio <- prediction
  
  # generate the calibration curve
  xBreaks <- molecule_df$l10_conc
  xLabs <- 10^xBreaks
  p <-    
    molecule_df %>% 
    ggplot(aes(y = l10_area_ratio, x = l10_conc)) +
    geom_line(data = pred_df) +
    # plot included points
    geom_point(
      data = molecule_df %>% filter(!is.na(l10_area_ratio),include_point == "include"),
      pch = 21, stroke = 1, fill = "lightgrey", size = 2.5) +
    # plot excluded points
    geom_point(  
      data = molecule_df %>% filter(!is.na(l10_area_ratio),include_point == "exclude"),
      pch = 21, fill = "white", size = 2.5)+
    # add the linear model
    ylab("Area ratio [log10]") +
    xlab("Concentration [uM]") +
    ggtitle(str_to_title(molecule_name)) +
    scale_x_continuous(breaks = xBreaks, labels = signif(xLabs, 2)) +
    theme_fs
  
  # Add the formula of the calibration curve and R2
  # text <- paste("y =", intercept, "+", slope, "* x","; R2 =", r_squared)
  text <- paste("R^2 =", r_squared)
  p <- p +
    geom_text(
      data = data.frame(x = -Inf, y = max_y, label = text),
      x = -Inf, y = max_y, label = text, hjust = -0.2, vjust = 0.5)
  
  # for the output of the model: generate a tibble with slope, intercept and r_squared
  mod_df <- tibble(slope = slope, intercept = intercept, r_squared = r_squared)
  
  res_list <- list(plot = p, model_df = mod_df)
  return(res_list)
}

f_generate_calibration_curves <- function(ratios_calibration_curve_df, save_plot_folder, coefficients_file_name=NULL,points_to_exclude = NULL) {
  # Computes calibration curves of all UNLABELLED molecules based on their area ratios computed previously
  # the save_plot_folder
  c <- 1

  if(!("area_ratio" %in% colnames(ratios_calibration_curve_df))){
    stop("Provided dataframe does not contain the column 'area_ratio'")
  }
  
  if(is.null(coefficients_file_name)){
    out_file_name <- file.path(save_plot_folder,"_model_coefficients.tsv")
  }else{
    out_file_name <- coefficients_file_name
  }
  
  if(!(dir.exists(save_plot_folder))){dir.create(save_plot_folder,recursive = T)}
  
  all_molecules <- unique(ratios_calibration_curve_df$name_clean)
  message("Processing ",length(all_molecules)," molecules")

  pb <- progress::progress_bar$new(total = length(all_molecules))
  model_coefficients_df <- tibble()
  mol <- all_molecules[2]
  for (mol in all_molecules) {
    message("\n",mol)
    c_df <- ratios_calibration_curve_df %>%
      filter(name_clean == mol) %>%
      filter(type == "unlabelled")    
    if (nrow(c_df) == 0 | !any(!is.na(c_df$area_ratio))) {
      message("No non-NA area-ratios for ", mol, " - skipping")
      # skip if dataframe is empty
      pb$tick()
      next
    }

    # forward points to exclude (if any)
    if (!is.null(points_to_exclude) & mol %in% names(points_to_exclude)) {
      # check whether the points to exclude are in the dataframe
      points_to_exclude_vec <- points_to_exclude[[mol]]
    }else{
      points_to_exclude_vec <- c()
    }

    # Generate the calibration curve and the linear model
    c_res <- f_plot_calibration_curve(molecule_df = c_df,points_to_exclude_vec = points_to_exclude_vec)
    # extract the model
    c_model_df <- c_res$model_df
    # combine the model df
    model_coefficients_df <- bind_rows(model_coefficients_df, c_model_df %>%
      mutate(
        name_clean = mol,
        name_original = unique(c_df$name_original)
      ))
    # save the plot
    # extract the calibration curve
    c_plot <- c_res$plot + theme(aspect.ratio = 1)
    ggsave(c_plot, filename = file.path(save_plot_folder, paste0(mol, ".pdf")), width = 5, height = 5)
    pb$tick()
  }
  # save the dataframe with all model coefficients and R2 values:
  write_tsv(model_coefficients_df %>% dplyr::relocate(name_original,name_clean), file = out_file_name)
  
  message("-------------- Finished --------------")
  message("Plots_saved_to ", save_plot_folder)
  
  message("Model coefficients saved to file: ", out_file_name)


}

f_compute_sample_concentration <- function(unknown_metabolites_df, coefficients_file_path) {
  #* Computes concentrations of the provided samples based on the saved model coefficients in the file ----# 
  # returns dataframe with calculated concentrations

  if(!file.exists(coefficients_file_path)){
    stop("The provided coefficients file does not exist")
  }else{
    coefficients_df <- read_tsv(coefficients_file_path) %>% dplyr::select(-name_original)
  }

  molecule_vec <- unique(unknown_metabolites_df$name_clean)

  message("Processing ",length(molecule_vec)," molecules")

  # keep only molecules presented in the coefficients file
  molecule_vec <- molecule_vec[molecule_vec %in% coefficients_df$name_clean]
  message(length(molecule_vec)," molecules present in the coefficients file - skipping the others")

  calculated_concentrations_df <- tibble()
  mol <- molecule_vec[2]
  for (mol in molecule_vec) {
    message(mol)
    c_df <- unknown_metabolites_df %>%
      filter(name_clean == mol) %>%
      filter(type == "unlabelled")    

    if (nrow(c_df) == 0 | !any(!is.na(c_df$area_ratio))) {
      message("No non-NA area-ratios for ", mol, " - skipping")
      # skip if dataframe is empty      
      next
    }
    
    # Join the coefficients to the dataframe
    c_df <- suppressMessages(c_df %>%
      inner_join(., coefficients_df) %>% 
      mutate(l10_area_ratio = log10(area_ratio)))

    # The formula from the calibration curve is:
    # y = mx + b with y being the area_ratio and x the concentration, m the slope and b the intercept
    # The concentration can be computed as: x = (y-b)/m

    # Compute concentrations based on the area_ratios from the given sample
    #! Important: Compute log10 area ratios since the linear-model parameters were also computed in log-10-space
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

f_export_to_excel <- function(dataframe,out_file_name){
  #takes a dataframe and exports it to .xlsx

  #add .xlsx suffix if not installed
  if(!str_detect(out_file_name,".xlsx$")){
    out_file_name <- paste0(out_file_name,".xlsx")
  }

  require(openxlsx)

  wb <- createWorkbook()
  addWorksheet(wb, "Sheet1")
  writeData(wb, "Sheet1", dataframe)
  
  # add column filter
  addFilter(wb, "Sheet1", cols = 1:ncol(dataframe),row = 1)

  # make columns as wide as header
  setColWidths(wb, "Sheet1", cols = 1:ncol(dataframe), widths = "auto")
  saveWorkbook(wb, file = out_file_name, overwrite = TRUE)
  message("Dataframe saved to: ",out_file_name)

}
