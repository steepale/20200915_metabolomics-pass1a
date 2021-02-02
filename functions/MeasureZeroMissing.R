MeasureZeroMissing <- function(RUN_N,
                       tissue,
                       shared_runs,
                       shared_runs_out,
                       shared_mets,
                       shared_sams,
                       nxp_out,
                       plot_df,
                       plot_join_df,
                       norm_df, 
                       norm_join_df,
                       zero_rm = 0.95, na_rm = 0.6, save_file = FALSE){
  # Identify Zero/Missing Values in Metabolites
  ################################################################################
  na_df <- plot_df %>%
    ungroup() %>% 
    group_by(METABOLITE_NAME, STUDY_INSTITUTE_METAB_FAMILY) %>%
    mutate(ZERO_LOG = ifelse(VALUE == 0, 1, 0)) %>%
    mutate(ZERO_LOG = ifelse(is.na(VALUE), 0, ZERO_LOG)) %>%
    mutate(ZERO_N = sum(ZERO_LOG)) %>%
    mutate(ZERO_FREQ = round(ZERO_N/length(shared_sams), digits = 2)) %>%
    mutate(NA_LOG = ifelse(is.na(VALUE), 1, 0)) %>%
    mutate(NA_N = sum(NA_LOG)) %>%
    mutate(NA_FREQ = round(NA_N/length(shared_sams), digits = 2)) %>%
    mutate(MIN = round(min(VALUE, na.rm = T)), digits = 2) %>%
    mutate(MAX = round(max(VALUE, na.rm = T)), digits = 2) %>%
    mutate(MEAN = round(mean(VALUE, na.rm = T)), digits = 2) %>%
    mutate(MEDIAN = round(median(VALUE, na.rm = T)), digits = 2) %>%
    mutate(Q1 = quantile(VALUE, na.rm = TRUE, probs = 0.25) %>% round(digits = 2)) %>%
    mutate(Q3 = quantile(VALUE, na.rm = TRUE, probs = 0.75) %>% round(digits = 2)) %>%
    select(TISSUE,STUDY_INSTITUTE_METAB_FAMILY, METABOLITE_NAME, ZERO_N, ZERO_FREQ, NA_N, NA_FREQ, 
           MEAN, MIN,Q1, MEDIAN, Q3, MAX) %>%
    unique() %>%
    arrange(STUDY_INSTITUTE_METAB_FAMILY, desc(NA_FREQ), desc(ZERO_FREQ)) %>%
    ungroup()

# Remove Metabolites with high NAor Zero Frequencies
################################################################################
met_rm <- na_df %>%
  filter((ZERO_FREQ >= zero_rm) | (NA_FREQ >= na_rm)) %>%
  select(METABOLITE_NAME) %>% unlist() %>% unique()
shared_mets <- shared_mets[shared_mets %!in% met_rm]
plot_df <- plot_df %>%
  mutate(METABOLITE_NAME = as.character(METABOLITE_NAME)) %>%
  filter(METABOLITE_NAME %in% shared_mets)
plot_join_df <- plot_join_df %>%
  mutate(METABOLITE_NAME = as.character(METABOLITE_NAME)) %>%
  filter(METABOLITE_NAME %in% shared_mets)
norm_df <- norm_df %>%
  mutate(METABOLITE_NAME = as.character(METABOLITE_NAME)) %>%
  filter(METABOLITE_NAME %in% shared_mets)
norm_join_df <- norm_join_df %>%
  mutate(METABOLITE_NAME = as.character(METABOLITE_NAME)) %>%
  filter(METABOLITE_NAME %in% shared_mets)

# Label metabolite that require imputation
met_imp <- na_df %>%
  filter(NA_FREQ > 0) %>%
  select(METABOLITE_NAME) %>% unlist() %>% as.character()
  
print(paste0("Metabolites Removed: ", paste(met_rm, collapse = '; ')))
print(paste0("Metabolites to Impute: ", paste(met_imp, collapse = '; ')))

if(save_file == TRUE){
  # Save the dataframe as an R object
  ####################################
  data_dir <- paste0(WD,"/data/site_comparisons/",shared_runs_out)
  cmd <- paste0("mkdir -p ",data_dir)
  system(cmd)
  na_file <- paste0(data_dir,"/20201104_pass1a-metabolomics-",tissue,"-zero-na-df_steep.txt")
  write.table(na_df, file = na_file, quote = F, row.names = F, sep = '\t')
  print(na_file)
}

  # Load the out list
      out_list <- list()
      out_list[['RUN_N']] <- RUN_N
      out_list[['tissue']] <- tissue
      out_list[['shared_runs']] <- shared_runs 
      out_list[['shared_runs_out']] <- shared_runs_out
      out_list[['shared_mets']] <- shared_mets
      out_list[['shared_sams']] <- shared_sams
      out_list[['nxp_out']] <- nxp_out
      out_list[['plot_df']] <- plot_df
      out_list[['plot_join_df']] <- plot_join_df
      out_list[['norm_df']] <- norm_df
      out_list[['norm_join_df']] <- norm_join_df
      out_list[['met_rm']] <- met_rm
      out_list[['na_df']] <- na_df
      out_list[['met_imp']] <- met_imp
      # Return the out list
      out_list
}
