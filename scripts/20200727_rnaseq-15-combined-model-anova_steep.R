
#' ## Goals of Analysis:
#' 
#' ## Setup the Environment

#+ Setup Environment, message=FALSE, results='hide', warning = FALSE
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("gam")
#install.packages("XML")
# for(p in packages){
#         if(!require(p, character.only = T)){
#                 install.packages(p)
#         }
#         library(p, character.only = T)
# }

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","reshape2","xtable","kohonen","caret","enrichR","gplots","tiff","splines","gam","DESeq2","car","KEGGREST")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
slice <- dplyr::slice
filter <- dplyr::filter
counts <- DESeq2::counts
map <- purrr::map
seq_range <- modelr::seq_range

# Global options
options(dplyr.print_max = 100)

# Source the functions
source(paste0(WD,'/functions/not_in.R'))
source(paste0(WD,'/functions/rat_mouse_ortho.R'))
source(paste0(WD,'/functions/mouse2rat_ortho.R'))
source(paste0(WD,'/functions/lmp.R'))
source(paste0(WD,'/functions/cor_PC_1_6.R'))
source(paste0(WD,'/functions/elbow_finder.R'))
source(paste0(WD,'/functions/cor_outlier2.R'))
source(paste0(WD,'/functions/sin.R'))
source(paste0(WD,'/functions/cos.R'))
source(paste0(WD,'/functions/circleFun.R'))

# A function to collect the non-parametric p-value
################################################################################
p.test.cos <- function(data,TOD,iter,every=10) {
        pval <- rep(0,dim(data)[1])
        cat(iter, "permutations in progress\n")
        PVE<-pve.FR(data,TOD)
        for(i in 1:iter) {
                # Shuffles the samples
                new.data<-data[,sample(1:dim(data)[2])]
                # Calculate the proportion of variance explained
                tmp.PVE <- pve.FR(new.data,TOD)
                # Generates a running count of simulated PVE 
                # that are larger or equal to the actual PVE
                pval <- pval + as.numeric(tmp.PVE >= PVE)
                if(i %% every == 0) cat(every)
                if(i %% (10*every) == 0) cat("\n")
        }
        # Calculate the p value by dividing the running count by the number of iterations
        pval <- pval / iter
}
################################################################################

# A function to collect the non-parametric p-value
################################################################################
p.test.ce <- function(data,TOD,iter,every=10) {
        pval <- rep(0,dim(data)[1])
        cat(iter, "permutations in progress\n")
        PVE<-pve.FR(data,TOD)
        for(i in 1:iter) {
                # Shuffles the samples
                new.data<-data[,sample(1:dim(data)[2])]
                # Calculate the proportion of variance explained
                tmp.PVE <- pve.FR(new.data,TOD)
                # Generates a running count of simulated PVE 
                # that are larger or equal to the actual PVE
                pval <- pval + as.numeric(tmp.PVE >= PVE)
                if(i %% every == 0) cat(every)
                if(i %% (10*every) == 0) cat("\n")
        }
        # Calculate the p value by dividing the running count by the number of iterations
        pval <- pval / iter
}
################################################################################


# Collect the percent variance explained 
################################################################################
PVE_e <- function(df){
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output pve exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained 
################################################################################
PVE_e2 <- function(df){
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output pve exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained 
################################################################################
PVE_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_res <- ss_res/ss_total
        # Output all 3 PVE (circ, exer, res)
        pve_circ
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb_e <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb2_e <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_circ
}
################################################################################


# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb2_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_circ
}
################################################################################

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
# Declare Tissue

# SCN: Hypothalamus (Suprachiasmatic nucleus)
# LIV: Liver
# KID: Kidney
# AOR: Aorta
# SKM: Gastrocnemius
# HAT: Heart
# ADG: Adrenal gland
# BAT: Brown adipose tissue
# WAT: White adipose tissue
# COR: Cortex
# HIP: Hippocampus
# LUNG: Lung
# OVR: Ovaries
# SPL: Spleen
# TES: Testes

# Load the decision table
table_file <- paste0(WD,'/data/20200603_rnaseq-tissue-data-assambly-table_steep.txt')
df_tbl <- read.table(file = table_file,sep = '\t', header = T, check.names = F)

by_gene_by_tissue_df <- tibble()
models_df <- data.frame()
#TISSUE <- "Kidney"
# for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
for(TISSUE in c('White Adipose')){
        print(TISSUE)
        TISSUE1 <- TISSUE
        # # Collect the formula
        # design <- df_tbl %>%
        # filter(Tissue == TISSUE) %>%
        # select(Formula) %>% unique() %>% 
        # unlist() %>% as.character() %>% as.formula()
        # Collect the Outliers
        OUTLIERS <- df_tbl %>%
                filter(Tissue == TISSUE) %>%
                select(Outliers) %>% unique() %>% 
                unlist() %>% as.character()
        # Collect Adjusted_Variance
        ADJ_VAR <- df_tbl %>%
                filter(Tissue == TISSUE) %>%
                select(Adjusted_Variance) %>% unique() %>% 
                unlist() %>% as.character()
        # Collect the TIS Symbol
        TIS <- df_tbl %>%
                filter(Tissue == TISSUE) %>%
                select(Tis) %>% unique() %>% 
                unlist() %>% as.character()
        # Collect the Formula
        FORMULA <- df_tbl %>%
                filter(Tissue == TISSUE) %>%
                select(Formula) %>% unique() %>% 
                unlist() %>% as.character() %>% as.formula()
        
        #' ## Load & Clean Data
        #' ##### Data files to load:
        #' * Count Matrix and Metadata Table from:
        #'     * RNA-Seq from Mt. Sinai
        #'         * 3 sequencing batches & metadata
        #'     * RNA-Seq from Stanford
        #'         * 2 sequencing batches & metadata
        
        #+ Load the Data
        ################################################################################
        #####     Load & Clean Data      ###############################################
        ################################################################################
        
        # TODO: Take this if statement from this script and incorporate it into the bollinger script
        # Files last saved in: 20200603_rnaseq-tissue-data-assembly_steep.R
        
        # Set a vector for Exercise/Control Levels and Colors
        ec_levels <- c('Exercise - IPE',
                       'Exercise - 0.5 hr',
                       'Exercise - 1 hr',
                       'Exercise - 4 hr',
                       'Exercise - 7 hr',
                       'Exercise - 24 hr',
                       'Exercise - 48 hr',
                       'Control - IPE',
                       'Control - 7 hr')
        ec_colors <- c('gold',
                       'darkgoldenrod1',
                       'orange',
                       'darkorange',
                       'darkorange2',
                       'darkorange3',
                       'darkorange4',
                       'steelblue1',
                       'steelblue4')
        
        # Load Metadata and count data as R objects
        ################################################################################
        # Restore the metadata object
        meta_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
        col_data <- readRDS(file = meta_file)
        
        # Restore the count object
        count_file <- paste0(WD,'/data/jiayu_subset/all_metabolites_normalized.txt')
        count_file1 <- paste0(WD,'/data/jiayu_subset/single_metabolte_2-hydroxybutyrate.csv')
        count_file2 <- paste0(WD,'/data/jiayu_subset/single_metabolte_N-acetylphenylalanine.csv')
        count_data <- read.table(count_file, header = T, sep = '\t', row.names = 1, 
                                 check.names = F)
        
        # Determine the median values by which to center
        median_c0 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Control - IPE') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e24 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 24 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e48 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 48 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e0.5 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 0.5 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e0 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - IPE') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e1 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 1 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e4 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 4 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_c7 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Control - 7 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        median_e7 <- col_data %>%
                filter(Tissue == TISSUE1) %>%
                filter(animal.key.anirandgroup == 'Exercise - 7 hr') %>%
                select(specimen.collection.t_death_hour) %>%
                unlist() %>% as.numeric() %>% median()
        
        # Apply the median centered groups to time of death (modeling purposes)
        col_data <- col_data %>%
                mutate(specimen.collection.t_death_hour_mc = 
                               case_when(animal.key.anirandgroup == 'Control - IPE' ~
                                                 median_c0,
                                         animal.key.anirandgroup == 'Control - 7 hr' ~
                                                 median_c7,
                                         animal.key.anirandgroup == 'Exercise - IPE' ~
                                                 median_e0,
                                         animal.key.anirandgroup == 'Exercise - 0.5 hr' ~
                                                 median_e0.5,
                                         animal.key.anirandgroup == 'Exercise - 1 hr' ~
                                                 median_e1,
                                         animal.key.anirandgroup == 'Exercise - 4 hr' ~
                                                 median_e4,
                                         animal.key.anirandgroup == 'Exercise - 7 hr' ~
                                                 median_e7,
                                         animal.key.anirandgroup == 'Exercise - 24 hr' ~
                                                 median_e24,
                                         animal.key.anirandgroup == 'Exercise - 48 hr' ~
                                                 median_e48))
        
        #' ## Collect Samples of Interest and Normalize
        
        #+ Collect Samples of Interest and Normalize
        ################################################################################
        #####     Collect Samples of Interest and Normalize      #######################
        ################################################################################
        if(TISSUE == c('Gastrocnemius')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Gastrocnemius') %>%
                        filter(Seq_batch == 'MSSM_1') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS)
        }else if(TISSUE == c('Lung')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Lung') %>%
                        filter(Seq_batch == 'MSSM_3') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS)
        }else{
                # Filter Samples (meta)
                tod_cols <- col_data %>%
                        filter(Tissue == TISSUE) %>%
                        filter(!is.na(animal.registration.sex))# %>%
                        #filter(sample_key %!in% OUTLIERS)
        }
        rownames(tod_cols) <- tod_cols$sample_key
        
        
        # Time post exercise
        tod_cols <- tod_cols %>%
                mutate(specimen.collection.t_exercise_hour = case_when(
                        animal.key.anirandgroup == 'Control - IPE' ~ -1,
                        animal.key.anirandgroup == 'Control - 7 hr' ~ 7,
                        animal.key.anirandgroup == 'Exercise - IPE' ~ 0,
                        animal.key.anirandgroup == 'Exercise - 0.5 hr' ~ 0.5,
                        animal.key.anirandgroup == 'Exercise - 1 hr' ~ 1,
                        animal.key.anirandgroup == 'Exercise - 4 hr' ~ 4,
                        animal.key.anirandgroup == 'Exercise - 7 hr' ~ 7,
                        animal.key.anirandgroup == 'Exercise - 24 hr' ~ 24,
                        animal.key.anirandgroup == 'Exercise - 48 hr' ~ 48))
        
        # Take the absolute value of the square root of seconds post exercise (consider negative numbers)
        # Make sure to Subtract 1 hour (3600s) from "Control - IPE" groups to account for exercise effect
        # tod_cols <- tod_cols %>%
        #         mutate(calculated.variables.deathtime_after_acute =
        #                        ifelse(animal.key.anirandgroup == 'Control - IPE', 
        #                               calculated.variables.deathtime_after_acute - 3600,
        #                               calculated.variables.deathtime_after_acute))
        # tod_cols <- tod_cols %>%
        #         mutate(specimen.collection.t_exercise_hour_sqrt = ifelse(
        #                 calculated.variables.deathtime_after_acute < 0, 
        #                 (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)*(-1), 
        #                 (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)))
        
        row.names(tod_cols) <- tod_cols$sample_key
        # # Examine histograms
        # tod_cols %>%
        #         filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        #         ggplot(aes(x=calculated.variables.deathtime_after_acute)) +
        #         geom_histogram(bins = 68)
        # tod_cols %>%
        #         filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        #         ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt)) +
        #         geom_histogram(bins = 68)
        
        # Collect samples without NA values in TOD
        nona_sams <- tod_cols %>%
                filter(!is.na(specimen.collection.t_death_hour)) %>%
                #filter(sample_key %!in% OUTLIERS) %>%
                filter(!is.na(animal.registration.sex)) %>%
                select(BID) %>% unlist() %>% as.character()
        nona_sams <- paste0(nona_sams, '015906')
        # Collect tissue specific counts
        nona_sams <- nona_sams[nona_sams %in% colnames(count_data)]
        tod_counts <- count_data[,nona_sams]
        # Subset the metadata as well
        rownames(tod_cols) <- paste0(tod_cols$BID, '015906')
        tod_cols <- tod_cols[nona_sams,]
        
        #' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
        all(rownames(tod_cols) == colnames(tod_counts))
        tod_cols$sample_key <- row.names(tod_cols)
        tod_counts$metabolite <- row.names(tod_counts)
        
        # Specify id.vars: the variables to keep but not split apart on
        t_counts <- reshape2::melt(tod_counts, id.vars=c("metabolite"))
        names(t_counts) <- c('metabolite','sample_key','count')
        tod_counts_t <- tod_counts %>%
                select(-metabolite) %>% t() %>% as.data.frame()
        tod_counts_t$sample_key <- row.names(tod_counts_t)
        
        #' #### Annotate Data for Modeling By Cluster
        
        #+ Annotate Data for Modeling By Cluster
        ################################################################################
        ################ Annotate Data for Modeling By Cluster #########################
        ################################################################################
        
        # Join the dataframes and nest
        by_gene_df <- left_join(tod_cols, t_counts, by = "sample_key") %>%
                filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
                group_by(metabolite) %>%
                arrange(sample_key) %>%
                nest()
        
        # Join the dataframes and nest
        by_gene_df7 <- left_join(tod_cols, t_counts, by = "sample_key") %>%
                filter(animal.key.anirandgroup %in% c('Control - 7 hr')) %>%
                group_by(metabolite) %>%
                arrange(sample_key) %>%
                nest()
        
        # Must be true
        all(by_gene_df$metabolite == by_gene_df7$metabolite)
        
        # Generate model functions for the dataframes
        gam_mod <- function(df) {
                lm(count ~ ns(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        poly4_mod <- function(df) {
                lm(count ~ poly(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        
        # Generate a model function for the dataframes
        sin_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc),
                   data = df)
        }
        
        # Generate a model that combines circadian with exercise (circadian first)
        ce_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc) +
                           ns(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        ce2_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc) +
                           poly(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        # Generate a model that combines circadian with exercise (exercise first)
        ec_mod <- function(df) {
                lm(count ~ ns(specimen.collection.t_exercise_hour_sqrt, df = 4) +
                           SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc), data = df)
        }
        ec2_mod <- function(df) {
                lm(count ~ poly(specimen.collection.t_exercise_hour_sqrt, df = 4) +
                           SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc), data = df)
        }
        
        
        
        # In case you'd like to subset the data
        #by_gene_df_bk <- by_gene_df
        #by_gene_df <- by_gene_df_bk
        
        by_gene_df <- by_gene_df %>%
                ungroup()
        
        # Capture the rejects mdl <- 'ce2'
        rejects_df <- data.frame()
        for(mdl in c('sin','poly4','ce2')){
                # circadian with exercise Model (circadian first)
                ################################################################################
                # Collect variables
                ( model_col <- as.symbol(paste0(mdl,'_model')) )
                anova_int <- as.symbol(paste0(mdl,'_anova_int'))
                anova_col <- as.symbol(paste0(mdl,'_anova'))
                resid_col <- as.symbol(paste0(mdl,'_resid'))
                metrics_col <- as.symbol(paste0(mdl,'_metrics'))
                summary_col <- as.symbol(paste0(mdl,'_summary'))
                deviance_col <- as.symbol(paste0(mdl,'_deviance'))
                MODEL <- match.fun(paste0(mdl,'_mod'))
                # Run models and save as a column
                by_gene_df <- by_gene_df %>%
                        mutate(!!model_col := map(data, MODEL))
                # Examine the ANOVA report on models
                #car::Anova defaults to type 2
                # Filter out genes with deviance of zero
                
                by_gene_df <- by_gene_df %>%
                        mutate(!!deviance_col := map_dbl(!!model_col, deviance))
                assign(paste0(mdl,'_rejects_df'), (by_gene_df %>%
                                                           filter(!!deviance_col <= sqrt(.Machine$double.eps))))
                
                rejects_df <- get(paste0(mdl,'_rejects_df'))
        }
        
        # Perform a series of analyses for each model
        for(mdl in c('sin','poly4','ce2')){
                # circadian with exercise Model (circadian first)
                ################################################################################
                # Collect variables
                ( model_col <- as.symbol(paste0(mdl,'_model')) )
                anova_int <- as.symbol(paste0(mdl,'_anova_int'))
                anova_col <- as.symbol(paste0(mdl,'_anova'))
                resid_col <- as.symbol(paste0(mdl,'_resid'))
                metrics_col <- as.symbol(paste0(mdl,'_metrics'))
                summary_col <- as.symbol(paste0(mdl,'_summary'))
                deviance_col <- as.symbol(paste0(mdl,'_deviance'))
                MODEL <- match.fun(paste0(mdl,'_mod'))
                # Run models and save as a column
                by_gene_df <- by_gene_df %>%
                        mutate(!!model_col := map(data, MODEL))
                # Examine the ANOVA report on models
                #car::Anova defaults to type 2
                # Filter out genes with deviance of zero
                
                by_gene_df <- by_gene_df %>%
                        mutate(!!deviance_col := map_dbl(!!model_col, deviance))
                by_gene_df <- by_gene_df %>%
                        filter(!!deviance_col > sqrt(.Machine$double.eps))
                by_gene_df <- by_gene_df %>%
                        mutate(!!anova_int := map(!!model_col, car::Anova)) %>%
                        mutate(!!anova_col := map(!!anova_int, broom::tidy)) %>%
                        select(-all_of(anova_int))
                # Add the residuals
                # by_gene_df <- by_gene_df %>%
                #   mutate(!!resid_col := map2(data, !!model_col, modelr::add_residuals))
                # # Examine the model metrics
                by_gene_df <- by_gene_df %>%
                   mutate(!!metrics_col := map(!!model_col, broom::glance))
                # # Examine some model summaries
                by_gene_df <- by_gene_df %>%
                   mutate(!!summary_col := map(!!model_col, summary))
        }
        
        # Collect all the PVEs
        pve_df <- by_gene_df %>%
                # select(ENSEMBL_RAT, SYMBOL_RAT,gam_anova,sin_anova,ce_anova,ec_anova,poly4_anova,ec2_anova,ce2_anova)
                select(metabolite,sin_anova,poly4_anova,ce2_anova)
        # Collect the PVE from models and different aspects of combined models
        pve_df <- pve_df %>%
                # mutate(pve_e = map_dbl(gam_anova, PVE_e)) %>%
                mutate(pve_e2 = map_dbl(poly4_anova, PVE_e2)) %>%
                mutate(pve_c = map_dbl(sin_anova, PVE_c)) %>%
                # mutate(pve_ec_c = map_dbl(ec_anova, PVE_comb_c)) %>%
                # mutate(pve_ec_e = map_dbl(ec_anova, PVE_comb_e)) %>%
                # mutate(pve_ce_c = map_dbl(ce_anova, PVE_comb_c)) %>%
                # mutate(pve_ce_e = map_dbl(ce_anova, PVE_comb_e)) %>%
                # mutate(pve_ec2_c = map_dbl(ec2_anova, PVE_comb2_c)) %>%
                # mutate(pve_ec2_e = map_dbl(ec2_anova, PVE_comb2_e)) %>%
                mutate(pve_ce2_c = map_dbl(ce2_anova, PVE_comb2_c)) %>%
                mutate(pve_ce2_e = map_dbl(ce2_anova, PVE_comb2_e)) %>%
                # select(-gam_anova,-sin_anova,-ce_anova,-ec_anova,
                #        -poly4_anova,-ce2_anova,-ec2_anova) %>%
                select(-sin_anova,-poly4_anova,-ce2_anova) %>%
                mutate(TISSUE = TISSUE1)
        pve_rejects_df <- data.frame(metabolite = rejects_df$metabolite,
                                     pve_e2 = 0,
                                     pve_c = 0,
                                     pve_ce2_c = 0,
                                     pve_ce2_e = 0,
                                     TISSUE = TISSUE1)
        pve_df <- rbind(pve_df, pve_rejects_df)
        
        # Concatenate the dataframes
        #models_df <- pve_df
        models_df <- rbind(models_df, pve_df)
        # Save dfs
        by_gene_df$TISSUE <- TISSUE1
        # by_gene_by_tissue_df <- rbind(by_gene_by_tissue_df, by_gene_df)
        by_gene_by_tissue_df_file <- paste0(WD,'/data/20200727_rnaseq-',TISSUE1,'-models-data_steep.rds')
        # saveRDS(by_gene_df, file = by_gene_by_tissue_df_file)
        # Save the final output table
        models_file <- paste0(WD,'/data/20200603_rnaseq-',TISSUE1,'-models-pve-table2_steep.txt')
        #write.table(models_df, file = models_file,sep = '\t',row.names = F,quote = F)
}
