# DEV
################
# imp_met <- 'Biliverdin'
# na_df <- out_list[['na_df']]
# plot_df <- out_list[['plot_df']]
#site <- shared_runs[1]

ImputeMissing <- function(RUN_N,
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
                       met_rm,
                       na_df,
                       met_imp,
                       imputation_method = 'half-min'){
  # Imputation
  #########################
  if(imputation_method == 'half-min'){
    for(imp_met in met_imp){
      for(site in shared_runs){
        # Collect minimum value
        min_val <- na_df %>%
          filter(METABOLITE_NAME == imp_met) %>%
          filter(STUDY_INSTITUTE_METAB_FAMILY == site) %>%
          select(MIN) %>% unlist() %>% as.numeric()
        # Simple Imputation
        plot_df <- plot_df %>%
          mutate(VALUE = ifelse(METABOLITE_NAME == imp_met & STUDY_INSTITUTE_METAB_FAMILY == site & is.na(VALUE), 
                                min_val/2, VALUE))
      } # Iterate over sites
    } # Iterate over metabolites that require imuptation
  } # Boolean, imputation decision
  
  # Data Output
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
