PrepareData2Sites <- function(comp_id, run_df, count_df, named = "named", 
                              feature_scale = 'linear_mean', feature_center  = 'linear_sd',
                              render = FALSE){

  # Ensure counts are in dataframe form
  if(!is.data.frame(count_df)){
    stop("Counts must be nested dataframe")
  }
  if(!is.character(comp_id)){
    stop("Comparison identifier must be a character string")
  }
  # Variables
  RUN_N <- run_df %>%
    filter(COMP_ID == comp_id) %>%
    select(RUN_N) %>% unlist() %>% as.numeric()
  # Ensure run comparison is for 2 sites
  if(RUN_N != 2){
      stop("Comparison must be 2 site comparison")
  }
  tissue <- run_df %>%
    filter(COMP_ID == comp_id) %>%
    select(TISSUE) %>% unlist() %>% as.character()
  shared_runs <- run_df %>%
    filter(COMP_ID == comp_id) %>%
    select(RUNS) %>% unlist() %>% as.character() %>% 
    str_split(pattern = ';') %>% unlist()
  shared_runs_out <- shared_runs %>% 
      str_replace_all('_','-') %>% str_replace_all(' ','-') %>% tolower() %>%
      paste(collapse = '_')
  shared_mets <- run_df %>%
    filter(COMP_ID == comp_id) %>%
    select(METABOLITES_RUN) %>% unlist() %>% as.character() %>% 
    str_split(pattern = ';') %>% unlist()
  shared_sams <- run_df %>%
    filter(COMP_ID == comp_id) %>%
    select(SAMPLES_RUN) %>% unlist() %>% as.character() %>% 
    str_split(pattern = ';') %>% unlist()
  nxp_out <- paste0('m',as.character(length(shared_mets)),'-s',as.character(length(shared_sams)))
   # Subset Data
  ################################################################################
  # Collect abundances in long format
  plot_df <- count_df %>%
      ungroup() %>%
      filter(TISSUE == tissue) %>%
      filter(NAMED == named) %>%
      unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
      filter(STUDY_INSTITUTE_METAB_FAMILY %in% shared_runs) %>%
      unnest(COUNT_DATA) %>%
      filter(METABOLITE_NAME %in% shared_mets) %>%
      filter(labelid %in% shared_sams)
  
      # Collect abundances in long "joined" format
      plot_join_df <- count_df %>%
        ungroup() %>%
        unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
        filter(TISSUE == tissue) %>%
        filter(NAMED == named) %>%
        filter(STUDY_INSTITUTE_METAB_FAMILY %in% shared_runs) %>%
        unnest(COUNT_DATA) %>%
        filter(METABOLITE_NAME %in% shared_mets) %>%
        filter(labelid %in% shared_sams) %>%
        select(STUDY_INSTITUTE_METAB_FAMILY,METABOLITE_NAME,labelid,VALUE) %>%
        mutate(STUDY_INSTITUTE_METAB_FAMILY = as.character(STUDY_INSTITUTE_METAB_FAMILY)) %>%
        split(.$STUDY_INSTITUTE_METAB_FAMILY) %>%
        map(~dplyr::rename(., !!sym(unique(.$STUDY_INSTITUTE_METAB_FAMILY)) := "VALUE")) %>%
        map(~select(., -STUDY_INSTITUTE_METAB_FAMILY)) %>%
        purrr::reduce(left_join, by = c("labelid","METABOLITE_NAME"))
      
      if(feature_scale == 'linear_mean' & feature_center  == 'linear_sd'){
        # The normalized counts
      norm_df <- plot_df %>%
        unite(labelid, STUDY_INSTITUTE_METAB_FAMILY, 
              col = "labelid_STUDY_INSTITUTE_METAB_FAMILY", remove = F,sep =';') %>%
        split(.$STUDY_INSTITUTE_METAB_FAMILY) %>%
        map(select, labelid_STUDY_INSTITUTE_METAB_FAMILY, METABOLITE_NAME, VALUE) %>%
        map(pivot_wider, names_from = METABOLITE_NAME, values_from = VALUE) %>%
        map(column_to_rownames, var = "labelid_STUDY_INSTITUTE_METAB_FAMILY") %>%
        map(as.matrix) %>%
        map(AutoScaleMatrix) %>%
        map(as.data.frame) %>%
        map(rownames_to_column, var = "labelid_STUDY_INSTITUTE_METAB_FAMILY") %>%
        map(pivot_longer, names_to = 'METABOLITE_NAME', values_to = 'VALUE', cols = all_of(shared_mets)) %>%
        map(separate, col = "labelid_STUDY_INSTITUTE_METAB_FAMILY", 
            into = c("labelid", "STUDY_INSTITUTE_METAB_FAMILY"), sep =';') %>%
        bind_rows()
      # Collect normalized abundances in long "joined" format
      norm_join_df <- plot_df %>%
        unite(labelid, STUDY_INSTITUTE_METAB_FAMILY, 
              col = "labelid_STUDY_INSTITUTE_METAB_FAMILY", remove = F,sep =';') %>%
        select(STUDY_INSTITUTE_METAB_FAMILY, labelid_STUDY_INSTITUTE_METAB_FAMILY, 
               METABOLITE_NAME, VALUE) %>%
        mutate(STUDY_INSTITUTE_METAB_FAMILY = as.character(STUDY_INSTITUTE_METAB_FAMILY)) %>%
        split(.$STUDY_INSTITUTE_METAB_FAMILY) %>%
        map(select, -STUDY_INSTITUTE_METAB_FAMILY) %>%
        map(pivot_wider, names_from = METABOLITE_NAME, values_from = VALUE) %>%
        map(column_to_rownames, var = "labelid_STUDY_INSTITUTE_METAB_FAMILY") %>%
        map(as.matrix) %>%
        map(AutoScaleMatrix) %>%
        map(as.data.frame) %>%
        map(rownames_to_column, var = "labelid_STUDY_INSTITUTE_METAB_FAMILY") %>%
        map(pivot_longer, names_to = 'METABOLITE_NAME', values_to = 'VALUE', cols = all_of(shared_mets)) %>%
        map(separate, col = "labelid_STUDY_INSTITUTE_METAB_FAMILY", 
            into = c("labelid", "STUDY_INSTITUTE_METAB_FAMILY"), sep =';') %>%
        map(~dplyr::rename(., !!sym(unique(.$STUDY_INSTITUTE_METAB_FAMILY)) := "VALUE")) %>%
        map(~select(., -STUDY_INSTITUTE_METAB_FAMILY)) %>%
        purrr::reduce(left_join, by = c("labelid","METABOLITE_NAME"))
      }
      
      if(render == TRUE){
        # Render the Primary Analysis
        ##########################################################################
        html_dir = '/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/site_comparisons'
        html_file = paste0("20201104_named-site-comparison-",comp_id,'_',tissue,'_',nxp_out,'_',
                         shared_runs_out,'_steep.html')
        rmarkdown::render('/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/scripts/20201104_pass1a-metabolomics-site-specific-eda-analysis-file_steep.Rmd',
        output_format = "html_document",
        output_file = html_file,
        output_dir = html_dir, 
        quiet = TRUE,
        clean = TRUE) 
        #print(paste0('Output File: ', html_file))
        print(paste0('Output Files PATH: ', html_dir))
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
      # Return the out list
      out_list
}
