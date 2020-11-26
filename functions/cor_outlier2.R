# Variable correlate with 'suspected' outlier status
################################################################################
cor_outlier2 <- function(rld = rld, ntop = 20000, intgroups = names(colData(rld))){
        # Create out dataframe
        out_df <- data.frame(matrix(ncol = 2, nrow = 0))
        names(out_df) <- c('Adjusted_R_Sq','Condition')
        for(intgroup in intgroups) {
                #intgroup <- 'animal.key.anirandgroup'
                group <- colData(rld)[[intgroup]]
                d <- colData(rld) %>% as_tibble()
                d <- d %>% mutate(sample_outlier = ifelse(sample_key %in% MIS_SUS, 1, 0))
                # Perform linear model on PC vs
                if(length(unique(group[!is.na(group)])) >=2){
                        int_r2 <- (lm(d[['sample_outlier']] ~ d[[intgroup]]) %>% summary())$adj.r.squared
                        int_df <- data.frame(Adjusted_R_Sq = int_r2,
                                             Condition = intgroup) %>% as_tibble()
                        out_df <- rbind(out_df, int_df)
                }else{
                        NA
                }
        }
        out_df %>%
                arrange(desc(Adjusted_R_Sq)) %>%
                filter(Adjusted_R_Sq >= 0.1) %>%
                filter(Condition  != 'sample_outlier')
}
################################################################################
