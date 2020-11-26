'---
#' title: "PASS1A Rat Tissue: -- ANOVA on Combined Model"
#' author: "Alec Steep" 
#' date: "20200727"
#' output:
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'     
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

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
pacs...man <- c("tidyverse","GenomicRanges","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","caret","enrichR","gplots","tiff","splines","gam","DESeq2","car","KEGGREST")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
slice <- dplyr::slice
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


by_gene_by_tissue_df <- tibble()
models_df <- data.frame()
#TISSUE <- "Kidney"
for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
        
        print(TISSUE)
        TISSUE1 <- TISSUE
        
        #' ## Load & Clean Data
        #' ##### Data files to load:
        
        
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
                        filter(!is.na(animal.registration.sex)) #%>%
                        #filter(sample_key %!in% OUTLIERS)
        }
        
        
        
        rownames(tod_cols) <- paste0(tod_cols$BID, '015906')
        
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


# Perform a series of analyses for each model
i <- 1
#mdl <- 'poly4'
for(mdl in c('sin','poly4','ce2')){
        # circadian with exercise Model (circadian first)
        ################################################################################
        # Collect variables
        print(i)
        i <- i + 1
        ( model_col <- as.symbol(paste0(mdl,'_model')) )
        anova_int <- as.symbol(paste0(mdl,'_anova_int'))
        anova_col <- as.symbol(paste0(mdl,'_anova'))
        resid_col <- as.symbol(paste0(mdl,'_resid'))
        metrics_col <- as.symbol(paste0(mdl,'_metrics'))
        summary_col <- as.symbol(paste0(mdl,'_summary'))
        MODEL <- match.fun(paste0(mdl,'_mod'))
        # Run models and save as a column
        by_gene_df <- by_gene_df %>%
                dplyr::mutate(!!model_col := map(data, MODEL))
        # Examine the ANOVA report on models
        #car::Anova defaults to type 2
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
        # select(metabolite, SYMBOL_RAT,gam_anova,sin_anova,ce_anova,ec_anova,poly4_anova,ec2_anova,ce2_anova)
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


# Visualize models on single gene
################################################################################
#################   Visualize models on single gene  ###########################
################################################################################

# hr7_df2 <- by_gene_df7 %>%
#         filter(metabolite == g) %>%
#         select(metabolite, SYMBOL_RAT, data) %>%
#         ungroup() %>%
#         unnest(data) %>%
#         filter(animal.key.anirandgroup == "Control - 7 hr") %>%
#         select(metabolite, SYMBOL_RAT, specimen.collection.t_exercise_hour, 
#                specimen.collection.t_death_hour_mc, count)
# sub_data2 <- sub_data %>%
#         mutate(specimen.collection.t_exercise_hour = 
#         ifelse(animal.key.anirandgroup == "Control - IPE", 0, specimen.collection.t_exercise_hour))
# sub_data2$specimen.collection.t_exercise_hour %>% table()

# ggplot() +
#         geom_point(data = sub_data2,
#                    aes(x = specimen.collection.t_exercise_hour, y = count,
#                        color = animal.key.anirandgroup)) +
#         geom_point(data = hr7_df2,
#                    mapping = aes(specimen.collection.t_exercise_hour, count),
#                    color = ec_colors[9]) +
#         geom_line(data = sub_data2, 
#                   aes(x = specimen.collection.t_exercise_hour, count), 
#                   size = 1, alpha = 0.6, color = 'blue',
#                   stat = "smooth", method = "lm", formula = y ~ ns(x, 5), se = F) +
#         scale_color_manual(values=ec_colors, drop = F) +
#         theme(legend.title = element_blank()) +
#         ggtitle(
#                 paste0(unique(model_pred_df$SYMBOL_RAT),
#                        " in ",TISSUE)) +
#         ylab("Counts (Normalized)") +
#         xlab("Hours Post Exercise")
# 
# ggplot() +
#         geom_point(data = sub_data2,
#                    aes(x = specimen.collection.t_death_hour_mc, y = count,
#                        color = animal.key.anirandgroup)) +
#         geom_point(data = hr7_df2,
#                    mapping = aes(specimen.collection.t_death_hour_mc, count),
#                    color = ec_colors[9]) +
#         scale_color_manual(values=ec_colors, drop = F) +
#         theme(legend.title = element_blank()) +
#         ggtitle(
#                 paste0(unique(model_pred_df$SYMBOL_RAT),
#                        " in ",TISSUE)) +
#         ylab("Counts (Normalized)") +
#         xlab("Hour of Death")


bk <- by_gene_df
#by_gene_df <- bk

models_df %>%
        arrange(pve_e2)

genes <- '2-hydroxybutyrate'
genes <- 'N-acetylphenylalanine'
genes <- 'hydroxyphenyllactic acid'

by_gene_df <- bk %>%
        filter(metabolite == genes)

# Visualize Genes in Single models
poly4_title <- "Counts ~ B0 + (B1...B4)*poly(sqrt(HPE),4)"
sin_title <- "Counts ~ B0 + B1*SIN(TOD) + B2*COS(TOD)"
ce2_title <- "Counts ~ B0 + (B1...B4)*poly(sqrt(HPE),4) +\nB5*SIN(TOD) + B6*COS(TOD)"
# mdl <- 'poly4'
#summary(mod)
by_gene_df$sin_anova
g <- genes
# bin_df %>%
#         filter(TISSUE == 'Kidney') %>%
#         filter(SYMBOL_RAT == 'Arntl')
mdl <- 'sin'
mdl <- 'ce2'
mdl <- 'poly4'

for(mdl in c('sin','poly4','ce2')){
        print(mdl)
        # Collect variables
        ###########################
        # 1. Add the predictions
        ######################################################
        genes <- by_gene_df$metabolite %>% as.character() %>% unique()
        # Must be true
        all(genes == by_gene_df$metabolite)
        ( model_pred <- as.symbol(paste0(mdl,'_pred')) )
        ( model_col <- as.symbol(paste0(mdl,'_model')) )
        ( model_pred_df <- as.symbol(paste0(mdl,'_pred_df')) )
        pred_list <- list()
        grid <- data.frame()
        
        for( g in genes){
                # Subset the original counts from the gene (Y original)
                sub_data <- (by_gene_df %>% 
                                     filter(metabolite == g) %>%
                                     ungroup() %>%
                                     select(data))[[1]] %>% as.data.frame() %>% 
                        arrange(specimen.collection.t_death_hour_mc)
                # Generate a grid of x values (grid X values) that are within the same range as actual x values
                if( mdl %in% c('gam', 'poly4')){
                        grid <- data.frame(
                                specimen.collection.t_exercise_hour_sqrt = 
                                        seq_range(sub_data$specimen.collection.t_exercise_hour_sqrt, 
                                                  n = length(sub_data$specimen.collection.t_exercise_hour_sqrt)))
                }else if( mdl == 'sin'){
                        grid <- data.frame(specimen.collection.t_death_hour_mc = 
                                                   seq_range(sub_data$specimen.collection.t_death_hour_mc, 
                                                             n = length(sub_data$specimen.collection.t_death_hour_mc)))
                }else if( mdl == 'ce2'){
                        # sort(sub_data$specimen.collection.t_death_hour_mc) %>% table()
                        grid <- data.frame()
                        for (tod_int in c('10.00-10.75','10.75-11.6','11.60-12.10','12.10-13.50',
                                          '13.50-14.10','14.10-14.80','14.80-17.80')){
                                #table(tod_cols$specimen.collection.t_death_hour_mc)
                                if(tod_int %in% c('10.00-10.75')){
                                        n_int = 12
                                }else if(tod_int %in% c('10.75-11.6')){
                                        n_int = 5
                                }else if(tod_int %in% c('11.60-12.10')){
                                        n_int = 7
                                }else if(tod_int %in% c('12.10-13.50')){
                                        n_int = 8
                                }else if(tod_int %in% c('13.50-14.10')){
                                        n_int = 7
                                }else if(tod_int %in% c('14.10-14.80')){
                                        n_int = 7
                                }else if(tod_int %in% c('14.80-17.80')){
                                        n_int = 4
                                }
                                tod_int_s <- (strsplit(tod_int,'-') %>% unlist())[1] %>% as.numeric()
                                tod_int_e <- (strsplit(tod_int,'-') %>% unlist())[2] %>% as.numeric()
                                tod_int_m <- median(c(tod_int_s,tod_int_e))
                                data_int <- sub_data %>%
                                        filter(specimen.collection.t_death_hour_mc >= tod_int_s &
                                                       specimen.collection.t_death_hour_mc < tod_int_e)
                                grid_int <- data.frame(specimen.collection.t_exercise_hour_sqrt = 
                                                               rep(median(data_int$specimen.collection.t_exercise_hour_sqrt), 
                                                                   n_int),
                                                       specimen.collection.t_death_hour_mc = 
                                                               rep(median(data_int$specimen.collection.t_death_hour_mc), 
                                                                   n_int))
                                grid <- rbind(grid, grid_int)
                        }
                }
                # Collect the actual model
                mod <- (by_gene_df %>% 
                                filter(metabolite == g) %>% ungroup() %>% 
                                select(model_col))[[1]][[1]]
                if( mdl == 'ce2'){
                        # Collect the coefficients
                        res_co <- (mod$coefficients %>%
                                           as.data.frame())['(Intercept)',]
                        sin_co <- (mod$coefficients %>%
                                           as.data.frame())['SIN(specimen.collection.t_death_hour_mc)',]
                        cos_co <- (mod$coefficients %>%
                                           as.data.frame())['COS(specimen.collection.t_death_hour_mc)',]
                        poly4_co1 <- (mod$coefficients %>%
                                              as.data.frame())['poly(specimen.collection.t_exercise_hour_sqrt, df = 4)1',]
                        poly4_co2 <- (mod$coefficients %>%
                                              as.data.frame())['poly(specimen.collection.t_exercise_hour_sqrt, df = 4)2',]
                        poly4_co3 <- (mod$coefficients %>%
                                              as.data.frame())['poly(specimen.collection.t_exercise_hour_sqrt, df = 4)3',]
                        poly4_co4 <- (mod$coefficients %>%
                                              as.data.frame())['poly(specimen.collection.t_exercise_hour_sqrt, df = 4)4',]
                        mod_c <- (by_gene_df %>% 
                                          filter(metabolite == g) %>% ungroup() %>% 
                                          select(sin_model))[[1]][[1]]
                        mod_c$coefficients[[1]] <- res_co
                        mod_c$coefficients[[2]] <- sin_co
                        mod_c$coefficients[[3]] <- cos_co
                        mod_e <- (by_gene_df %>% 
                                          filter(metabolite == g) %>% ungroup() %>% 
                                          select(poly4_model))[[1]][[1]]
                        mod_e$coefficients[[1]] <- res_co
                        mod_e$coefficients[[2]] <- poly4_co1
                        mod_e$coefficients[[3]] <- poly4_co2
                        mod_e$coefficients[[4]] <- poly4_co3
                        mod_e$coefficients[[5]] <- poly4_co4
                        # Add the residual to the grid
                        grid$res_co <- res_co
                }
                # Add the predictions (YTotal) to the grid (from grid X values) using the model
                grid <- modelr::add_predictions(grid, mod, "pred") %>% as_tibble()
                # Rename the grid X values
                if( mdl %in% c('gam', 'poly4')){
                        names(grid)[1] <- "grid_t_exercise_hour_sqrt"
                        # Add the observed x values to the grid dataframe
                        grid$specimen.collection.t_exercise_hour_sqrt <- 
                                sub_data$specimen.collection.t_exercise_hour_sqrt
                }else if( mdl == 'sin'){
                        names(grid)[1] <- "grid_t_death_hour"
                        # Add the observed x values to the grid dataframe
                        grid$grid_t_death_hour <- round(grid$grid_t_death_hour, digits = 1)
                        grid$specimen.collection.t_death_hour_mc <- 
                                sub_data$specimen.collection.t_death_hour_mc
                }else if( mdl == 'ce2'){
                        # Add additional predictions
                        grid <- modelr::add_predictions(grid, mod_c, "pred_c") %>% as_tibble()
                        grid <- modelr::add_predictions(grid, mod_e, "pred_e") %>% as_tibble()
                        names(grid)[1] <- 'grid_t_exercise_hour_sqrt'
                        names(grid)[2] <- 'grid_t_death_hour'
                        # Add the observed x values to the grid dataframe
                        grid$specimen.collection.t_exercise_hour_sqrt <- 
                                sub_data$specimen.collection.t_exercise_hour_sqrt
                        grid$grid_t_death_hour <- round(grid$grid_t_death_hour, digits = 1)
                        grid$specimen.collection.t_death_hour_mc <- 
                                sub_data$specimen.collection.t_death_hour_mc
                }
                # Add the observed count values to the grid dataframe
                grid$count <- sub_data$count
                # Add additional annotation
                grid$animal.key.anirandgroup <- factor(sub_data$animal.key.anirandgroup, 
                                                       levels = ec_levels)
                grid$animal.registration.sex <- sub_data$animal.registration.sex
                # ggplot() +
                #         geom_point(data = grid,
                #                    aes(x = specimen.collection.t_exercise_hour_sqrt, y = count)) +
                #         geom_point(data = grid,
                #                    aes(x = grid_t_exercise_hour_sqrt, y = pred), 
                #                        color = 'blue')
                # geom_point(data = grid,
                #            aes(x = grid_t_exercise_hour_sqrt, y = pred_c, 
                #                color = 'red')) +
                # geom_point(data = grid,
                #         aes(x = grid_t_exercise_hour_sqrt, y = pred_e, 
                #                 color = 'green'))
                # ggplot() +
                #         geom_point(data = grid,
                #                    aes(x = specimen.collection.t_death_hour_mc, y = count)) +
                #         geom_point(data = grid,
                #                    aes(x = grid_t_death_hour, y = pred),
                #                        color = 'blue') +
                #         geom_point(data = grid,
                #                    aes(x = grid_t_death_hour, y = pred_c),
                #                        color = 'red') +
                #         geom_point(data = grid,
                #                    aes(x = grid_t_death_hour, y = pred_e),
                #                    color = 'green')
                pred_list[[model_pred]][[g]] <- grid
        }
        
        #Add the predictions from the model to an object in the top tibble
        by_gene_df <- by_gene_df %>%
                mutate(!!model_pred := pred_list[[model_pred]])
        # Per gene visualization
        for( g in genes){
                # Visualize the GAM model by gene
                model_pred_df <- by_gene_df %>%
                        select(metabolite, !!model_pred) %>%
                        ungroup() %>%
                        filter(metabolite == g) %>%
                        unnest(!!model_pred)
                # Collect the Control 7 hr data points
                hr7_df <- by_gene_df7 %>%
                        filter(metabolite == g) %>%
                        select(metabolite, data) %>%
                        ungroup() %>%
                        unnest(data) %>%
                        filter(animal.key.anirandgroup == "Control - 7 hr") %>%
                        select(metabolite, specimen.collection.t_exercise_hour_sqrt, 
                               specimen.collection.t_death_hour_mc, count)
                # Visualize the models
                if(mdl == 'poly4'){
                        # Raw counts with line
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_exercise_hour_sqrt, y = count,
                                               color = animal.key.anirandgroup)) +
                                scale_color_manual(values=ec_colors) +
                                geom_line(data = model_pred_df,
                                          aes(x = grid_t_exercise_hour_sqrt, pred),
                                          size = 1, alpha = 0.6, color = "firebrick") +
                                theme(legend.title = element_blank()) +
                                ggtitle(
                                        paste0(unique(model_pred_df$metabolite),
                                               " in ",TISSUE,"\n",poly4_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Square Root of Hours Post Exercise")
                        # Raw counts with Line with control 7
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_exercise_hour_sqrt, y = count,
                                               color = animal.key.anirandgroup)) +
                                geom_point(data = hr7_df, 
                                           mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
                                           color = ec_colors[9]) +
                                scale_color_manual(values=ec_colors, drop = FALSE) +
                                geom_line(data = model_pred_df, 
                                          aes(x = grid_t_exercise_hour_sqrt, pred), 
                                          size = 1, alpha = 0.6, color = "firebrick") +
                                theme(legend.title = element_blank()) +
                                ggtitle(
                                        paste0(unique(model_pred_df$SYMBOL_RAT),
                                               " in ",TISSUE,"\n",poly4_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Square Root of Hours Post Exercise")
                        
                }else if(mdl == 'sin'){
                        # Raw counts with line
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_death_hour_mc, y = count,
                                               color = animal.key.anirandgroup)) +
                                scale_color_manual(values=ec_colors) +
                                geom_line(data = model_pred_df, 
                                          aes(x = grid_t_death_hour, pred), 
                                          size = 1, alpha = 0.6, color = 'darkblue') +
                                theme(legend.title = element_blank()) +
                                ggtitle(
                                        paste0(unique(model_pred_df$metabolite),
                                               " in ",TISSUE,"\n",sin_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Hour of Death")
                        # Raw counts with Line with control 7
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_death_hour_mc, y = count,
                                               color = animal.key.anirandgroup)) +
                                geom_point(data = hr7_df, 
                                           mapping = aes(specimen.collection.t_death_hour_mc, count),
                                           color = ec_colors[9]) +
                                scale_color_manual(values=ec_colors, drop = F) +
                                geom_line(data = model_pred_df, 
                                          aes(x = grid_t_death_hour, pred), 
                                          size = 1, alpha = 0.6, color = 'darkblue') +
                                theme(legend.title = element_blank()) +
                                ggtitle(paste0(unique(model_pred_df$SYMBOL_RAT),
                                               " in ",TISSUE,"\n",sin_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Hour of Death")
                }else if(mdl == 'ce2'){
                        # Raw counts with line (TOD on x-axis)
                        # ggplot() +
                        #         geom_point(data = model_pred_df,
                        #                    aes(x = specimen.collection.t_death_hour_mc, y = count,
                        #                        color = animal.key.anirandgroup)) +
                        #         geom_line(data = model_pred_df, 
                        #                   aes(x = grid_t_death_hour, pred), 
                        #                   size = 1, alpha = 0.6, color = 'blue') +
                        #         theme(legend.title = element_blank()) +
                        #         ggtitle(
                        #                 paste0(unique(model_pred_df$SYMBOL_RAT),
                        #                        " in ",TISSUE,"\n",ce2_title)) +
                        #         ylab("Counts (Normalized)") + 
                        #         xlab("Hour of Death")
                        # # Raw counts with Line with control 7 (TOD on x-axis)
                        # ggplot() +
                        #         geom_point(data = model_pred_df,
                        #                    aes(x = specimen.collection.t_death_hour_mc, y = count,
                        #                        color = animal.key.anirandgroup)) +
                        #         geom_point(data = hr7_df, 
                        #         mapping = aes(specimen.collection.t_death_hour_mc, count),
                        #         color = ec_colors[9]) +
                        #         scale_color_manual(values=ec_colors, drop = F) +
                        #         geom_line(data = model_pred_df, 
                        #                   aes(x = grid_t_death_hour, pred), 
                        #                   size = 1, alpha = 0.6, color = 'blue') +
                        #         theme(legend.title = element_blank()) +
                        #         ggtitle(
                        #                 paste0(unique(model_pred_df$SYMBOL_RAT),
                        #                        " in ",TISSUE,"\n",ce2_title)) +
                        #         ylab("Counts (Normalized)") + 
                        #         xlab("Hour of Death")
                        # Ys from different model parts visualized (TOD on x-axis)
                        rb_df_e <- data.frame(x = 10.8, y=11.1)
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_death_hour_mc, y = count,
                                               color = animal.key.anirandgroup)) +
                                # geom_point(data = hr7_df, 
                                #            mapping = aes(specimen.collection.t_death_hour_mc, 
                                #                          count),color = ec_colors[9]) +
                                geom_line(data = model_pred_df, 
                                          aes(x = grid_t_death_hour, pred_c), 
                                          size = 1, alpha = 0.6, color = 'blue',
                                          stat = "smooth", method = "loess", se = F) +
                                geom_point(data = model_pred_df, 
                                           aes(x = grid_t_death_hour, pred_e), 
                                           size = 2, alpha = 1, color = "firebrick") +
                                # geom_point(data = rb_df_e, 
                                #            aes(x = x, y =y), 
                                #            size = 2, alpha = 1, color = "firebrick") +
                                geom_abline(intercept = res_co, slope = 0, 
                                            linetype = "dashed", alpha = 0.8) +
                                scale_color_manual(values=ec_colors, drop = F) +
                                theme(legend.title = element_blank()) +
                                ggtitle(
                                        paste0(unique(model_pred_df$metabolite),
                                               " in ",TISSUE,"\n",ce2_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Hour of Death")
                        
                        # Raw counts with line (HPE on x-axis)
                        # ggplot() +
                        #         geom_point(data = model_pred_df,
                        #                    aes(x = specimen.collection.t_exercise_hour_sqrt, 
                        #                        y = count, color = animal.key.anirandgroup)) +
                        #         geom_line(data = model_pred_df, 
                        #                   aes(x = grid_t_exercise_hour_sqrt, pred), 
                        #                   size = 1, alpha = 0.6, color = 'blue',
                        #                   stat = "smooth", method = "loess", se = F) +
                        #         theme(legend.title = element_blank()) +
                        #         ggtitle(
                        #                 paste0(unique(model_pred_df$SYMBOL_RAT),
                        #                        " in ",TISSUE,"\n",ce2_title)) +
                        #         ylab("Counts (Normalized)") + 
                        #         xlab("Square Root of Hours Post Exercise")
                        # # Raw counts with Line with control 7 (HPE on x-axis)
                        # ggplot() +
                        #         geom_point(data = model_pred_df,
                        #                    aes(x = specimen.collection.t_exercise_hour_sqrt, 
                        #                        y = count,color = animal.key.anirandgroup)) +
                        #         geom_point(data = hr7_df, 
                        #                    mapping = aes(specimen.collection.t_exercise_hour_sqrt,
                        #                                  y = count), color = ec_colors[9]) +
                        #         scale_color_manual(values=ec_colors, drop = F) +
                        #         geom_line(data = model_pred_df, 
                        #                   aes(x = grid_t_exercise_hour_sqrt, pred), 
                        #                   size = 1, alpha = 0.6, color = 'blue',
                        #                   stat = "smooth", method = "loess", se = F) +
                        #         theme(legend.title = element_blank()) +
                        #         ggtitle(
                        #                 paste0(unique(model_pred_df$SYMBOL_RAT),
                        #                        " in ",TISSUE,"\n",ce2_title)) +
                        #         ylab("Counts (Normalized)") + 
                        #         xlab("Square Root of Hours Post Exercise")
                        # Ys from different model parts visualized (HPE on x-axis)
                        #rb_df_c <- data.frame(x = 0.0817, y=11.8)
                        summary(mod)
                        ggplot() +
                                geom_point(data = model_pred_df,
                                           aes(x = specimen.collection.t_exercise_hour_sqrt, 
                                               y = count,color = animal.key.anirandgroup)) +
                                # geom_point(data = hr7_df, 
                                #            mapping = aes(specimen.collection.t_exercise_hour_sqrt,
                                #                          y = count), color = ec_colors[9]) +
                                scale_color_manual(values=ec_colors, drop = F) +
                                # geom_point(data = rb_df_c,
                                #            aes(x = x, y =y),
                                #            size = 2, alpha = 1, color = "blue") +
                                geom_point(data = model_pred_df, 
                                           aes(x = grid_t_exercise_hour_sqrt, pred_c), 
                                           size = 2, alpha = 0.6, color = 'blue') +
                                geom_line(data = model_pred_df, 
                                          aes(x = grid_t_exercise_hour_sqrt, pred_e), 
                                          size = 1, alpha = 0.6, color = 'firebrick',
                                          stat = "smooth", method = "loess", se = F) +
                                geom_abline(intercept = res_co, slope = 0, 
                                            linetype = "dashed", alpha = 0.8) +
                                theme(legend.title = element_blank()) +
                                ggtitle(
                                        paste0(unique(model_pred_df$metabolite),
                                               " in ",TISSUE,"\n",ce2_title)) +
                                ylab("Counts (Normalized)") + 
                                xlab("Square Root of Hours Post Exercise")
                }
        }
        
        
        
        
        # Raw counts with Line with control 7
        ec_pred_df %>%
                ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
                       color = metabolite) +
                geom_point() +
                geom_line(data = ec_pred_df, 
                          aes(grid_t_exercise_hour_sqrt_jit, pred), 
                          size = 1, alpha = 0.8, color = "blue") +
                geom_point(data = ec_hr7_df, 
                           mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
                           color = "red") +
                theme(legend.position = "none") +
                ggtitle(
                        paste0("Expression of ",
                               unique(ec_pred_df$SYMBOL_RAT),
                               ":\nExercise Groups & Control IPE (",TISSUE,")")) +
                ylab("Expression (Transformed/Normalized)") + 
                xlab("Hours Post Acute Exercise (Transformed)") +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
        ec_mod
        # Raw counts with Line with control 7
        ec_pred_df %>%
                ggplot(aes(specimen.collection.t_death_hour_mc, count), 
                       color = metabolite) +
                geom_point() +
                geom_line(data = ec_pred_df, 
                          aes(grid_t_death_hour, pred), 
                          size = 1, alpha = 0.8, color = "orange") +
                geom_point(data = ec_hr7_df, 
                           mapping = aes(specimen.collection.t_death_hour_mc, count),
                           color = "red") +
                theme(legend.position = "none") +
                ggtitle(
                        paste0("Expression of ",
                               unique(ec_pred_df$SYMBOL_RAT),
                               ":\nExercise Groups & Control IPE (",TISSUE,")")) +
                ylab("Expression (Transformed/Normalized)") + 
                xlab("Time of Death (Hour)")
        
        # Combined Model (CE)
        ########################
        # Visualize the raw counts, model predictions, and control 7 counts (x=Hours Post Exercise)
        # Collect the Control 7 hr data points
        ce_hr7_df <- by_gene_df7 %>%
                filter(metabolite == ce_gene) %>%
                select(metabolite, SYMBOL_RAT, data) %>%
                ungroup() %>%
                unnest(data) %>%
                select(metabolite, SYMBOL_RAT, 
                       specimen.collection.t_exercise_hour_sqrt,
                       specimen.collection.t_death_hour_mc, count)
        # Collect the Control 0 hr data points
        ce_hr1_df <- by_gene_df %>%
                filter(metabolite == ce_gene) %>%
                select(metabolite, SYMBOL_RAT, data) %>%
                ungroup() %>%
                unnest(data) %>% 
                filter(animal.key.anirandgroup == "Control - IPE") %>%
                select(metabolite, SYMBOL_RAT, 
                       specimen.collection.t_exercise_hour_sqrt,
                       specimen.collection.t_death_hour_mc, count)
        
        # Raw counts with Line with control 7
        ce_pred_df %>%
                ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
                       color = metabolite) +
                geom_point() +
                geom_line(data = ce_pred_df, 
                          aes(grid_t_exercise_hour_sqrt_jit, pred), 
                          size = 1, alpha = 0.8, color = "blue") +
                geom_point(data = ce_hr7_df, 
                           mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
                           color = "red") +
                theme(legend.position = "none") +
                ggtitle(
                        paste0("Expression of ",
                               unique(ce_pred_df$SYMBOL_RAT),
                               ":\nExercise Groups & Control IPE (",TISSUE,")")) +
                ylab("Expression (Transformed/Normalized)") + 
                xlab("Hours Post Acute Exercise (Transformed)") +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
        ce_mod
        # Raw counts with Line with control 7
        ce_pred_df %>%
                ggplot(aes(specimen.collection.t_death_hour_mc, count), 
                       color = metabolite) +
                geom_point() +
                geom_line(data = ce_pred_df, 
                          aes(grid_t_death_hour, pred), 
                          size = 1, alpha = 0.8, color = "orange") +
                geom_point(data = ce_hr7_df, 
                           mapping = aes(specimen.collection.t_death_hour_mc, count),
                           color = "red") +
                theme(legend.position = "none") +
                ggtitle(
                        paste0("Expression of ",
                               unique(ce_pred_df$SYMBOL_RAT),
                               ":\nExercise Groups & Control IPE (",TISSUE,")")) +
                ylab("Expression (Transformed/Normalized)") + 
                xlab("Time of Death (Hour)")
        
        
        
        ################################################################################
        
        