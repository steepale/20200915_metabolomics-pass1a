# Function to grab topp 500 variance genes and find correlated variables
################################################################################
cor_PC_1_6 <- function(rld = rld, ntop = 20000, intgroups){
        rv <- rowVars(assay(rld))
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(assay(rld)[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        out_df <- data.frame(matrix(ncol = 3, nrow = 0))
        names(out_df) <- c('PC','Adjusted_R_Sq','Condition')
        for(intgroup in intgroups) {
                intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
                group <- colData(rld)[[intgroup]]
                d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],PC5 = pca$x[, 5],PC6 = pca$x[, 6],
                                group = group, 
                                intgroup.df, name = colnames(rld))
                # Perform linear model on PC vs
                if(length(unique(group[!is.na(group)])) >=2){
                        r2_pc1 <- (lm(PC1 ~ group, data = d) %>% summary())$adj.r.squared
                        r2_pc2 <- (lm(PC2 ~ group, data = d) %>% summary())$adj.r.squared
                        r2_pc3 <- (lm(PC3 ~ group, data = d) %>% summary())$adj.r.squared
                        r2_pc4 <- (lm(PC4 ~ group, data = d) %>% summary())$adj.r.squared
                        r2_pc5 <- (lm(PC5 ~ group, data = d) %>% summary())$adj.r.squared
                        r2_pc6 <- (lm(PC6 ~ group, data = d) %>% summary())$adj.r.squared
                        int_df <- data.frame(PC = 1:6, 
                                             Adjusted_R_Sq = c(r2_pc1, r2_pc2,r2_pc3, r2_pc4,r2_pc5,r2_pc6),
                                             Condition = intgroup) %>% as_tibble()
                        out_df <- rbind(out_df, int_df)
                }else{
                        NA
                }
        }
        out_df %>%
                arrange(desc(Adjusted_R_Sq))
        
}
################################################################################
