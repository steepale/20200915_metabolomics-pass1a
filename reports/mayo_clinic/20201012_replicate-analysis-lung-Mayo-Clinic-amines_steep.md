PASS1A (Rat) Metabolomics: Mayo Clinic Replicate Comparisons (2 reps)
================
Steep
11/26/2020

  - [Setup the Environment](#setup-the-environment)
  - [Per Run, Split the Data by duplicates and plot boxplots and matched
    density plot before
    normalization](#per-run-split-the-data-by-duplicates-and-plot-boxplots-and-matched-density-plot-before-normalization)
  - [Normalize and plot boxplots and matched density plot after
    normalization](#normalize-and-plot-boxplots-and-matched-density-plot-after-normalization)

## Setup the Environment

# Per Run, Split the Data by duplicates and plot boxplots and matched density plot before normalization

# Normalize and plot boxplots and matched density plot after normalization

``` r
# Assign Variables
################################################################################
tissue_study_institute_metab_families <- countdata_df %>%
  ungroup() %>%
  unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
  unite(TISSUE,STUDY_INSTITUTE_METAB_FAMILY, col = TISSUE_STUDY_INSTITUTE_METAB_FAMILY) %>%
  select(TISSUE_STUDY_INSTITUTE_METAB_FAMILY) %>% unlist() %>% as.character() %>% unique()

# DEV
# study_institute_metab_family <- "Mayo Clinic_ac"
# tissue <- 'liver'
# tissue_study_institute_metab_familiy <- tissue_study_institute_metab_families[1]
#for(tissue_study_institute_metab_family in tissue_study_institute_metab_families[1]){
  # Variables
tissue <- tissue_study_institute_metab_family %>% 
    str_split(pattern = '_') %>% unlist() %>% head(n=1)
study_institute_metab_family <- tissue_study_institute_metab_family %>%
    str_remove(paste0(tissue,'_'))
named <- 'named'
  
  # Collect the samples (remove QCs) & join sample order
sample_join <- countdata_df %>%
        ungroup() %>%
        unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
        filter(TISSUE == tissue) %>%
        filter(NAMED == named) %>%
        filter(STUDY_INSTITUTE_METAB_FAMILY == study_institute_metab_family) %>%
        unnest(SAMPLE_DATA) %>%
        filter(sample_type == 'Sample') %>%
        select(sample_id, sample_order)

  
# If the runs contain replicate samples, then proceed with the analysis
# Print the number of replicates
print(paste0('Number of replicates in ',tissue,' ',study_institute_metab_family,': ', length(reps)))
```

    ## [1] "Number of replicates in lung Mayo Clinic_amines: 2"

``` r
print(paste0('Replicates in ',tissue,' ',study_institute_metab_family,': ', paste(reps, collapse = ' ')))
```

    ## [1] "Replicates in lung Mayo Clinic_amines: 611 610"

``` r
# Collect the metabolites unique to run
metabolites <- countdata_df %>%
        ungroup() %>%
        unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
        filter(TISSUE == tissue) %>%
        filter(NAMED == named) %>%
        filter(STUDY_INSTITUTE_METAB_FAMILY == study_institute_metab_family) %>%
        unnest(METABOLITE_DATA) %>%
        select(refmet_name) %>% unlist() %>% as.character() %>% unique()

# Subset Data
################################################################################

# Collect abundances in long format and annotate replicate labelids
plot_df <- countdata_df %>%
        ungroup() %>%
        unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
        filter(TISSUE == tissue) %>%
        filter(NAMED == named) %>%
        filter(STUDY_INSTITUTE_METAB_FAMILY == study_institute_metab_family) %>%
        unnest(COUNT_DATA) %>%
        filter(viallabel %in% as.character(unlist(sample_join$sample_id))) %>%
        left_join(y = sample_join, by = c("viallabel" = "sample_id")) %>%
        unite(labelid, viallabel, col = "labelid_viallabel", remove = F) %>%
        mutate(REPLICATE = case_when(grepl(paste0(reps[1],"$"),labelid_viallabel) ~ reps[1],
                                     grepl(paste0(reps[2],"$"),labelid_viallabel) ~ reps[2],
                                     grepl(paste0(reps[3],"$"),labelid_viallabel) ~ reps[3]))

# Plot Boxplots faceted by shared metabolites
################################################################################
plot_df %>%
  ggplot(aes(y = METABOLITE_NAME, x = VALUE, color = REPLICATE)) +
  geom_boxplot(alpha = 0.8) +
  labs(title=paste0(tissue,' ',study_institute_metab_family,": Original Metabolite Abundances"),
             x = "Abundance", y = "") 
```

![](/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/mayo_clinic/20201012_replicate-analysis-lung-Mayo-Clinic-amines_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Normalization-1.png)<!-- -->

``` r
# Summary Statistics of Density plot distribution
################################################################################
# Iterate through the runs and collect vectors
df_all <- data.frame()
for(rep in c(reps)){
  # Collect numeric vector
  num_vec <- plot_df %>%
    filter(REPLICATE == rep) %>%
    select(VALUE) %>% unlist() %>% unname()
  df <- NumericSummaryStats(num_vec)
  row.names(df) <- rep
  df_all <- rbind(df_all, df)
}
df_all

# Plot the density plot for all the gene counts
################################################################################
plot_df %>%
  ggplot(aes(x = VALUE, color = REPLICATE)) +
  geom_density() +
  xlim(min(df_all$MIN), mean(df_all$Q3)) +
  labs(title=paste0(tissue,' ',study_institute_metab_family,": Original Metabolite Abundances"),
       x = "Abundance",
       y = "Density")
```

![](/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/mayo_clinic/20201012_replicate-analysis-lung-Mayo-Clinic-amines_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Normalization-2.png)<!-- -->

``` r
# Subset the data for sites to compare (scatterplots)
###################################################
plot_join_df <- countdata_df %>%
  ungroup() %>%
  unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
  filter(TISSUE == tissue) %>%
  filter(NAMED == named) %>%
  filter(STUDY_INSTITUTE_METAB_FAMILY == study_institute_metab_family) %>%
  unnest(COUNT_DATA) %>%
  filter(viallabel %in% as.character(unlist(sample_join$sample_id))) %>%
  unite(labelid, viallabel, col = "labelid_viallabel", remove = F) %>%
  mutate(REPLICATE = case_when(grepl(paste0(reps[1],"$"),labelid_viallabel) ~ reps[1],
                              grepl(paste0(reps[2],"$"),labelid_viallabel) ~ reps[2],
                              )) %>%
  select(REPLICATE,METABOLITE_NAME,labelid,VALUE) %>%
  mutate(REPLICATE = as.character(REPLICATE)) %>%
  split(.$REPLICATE) %>%
  map(~rename(., !!sym(unique(.$REPLICATE)) := "VALUE")) %>%
  map(~select(., -REPLICATE)) %>%
  purrr::reduce(left_join, by = c("labelid","METABOLITE_NAME"))

# Normalize
# Note in this step, we use labelid to account for mayo's duplicate samples
################################################################################
norm_df <- plot_df %>%
  unite(labelid, viallabel, col = "labelid_viallabel", remove = F) %>%
  split(.$REPLICATE) %>%
  map(select, labelid_viallabel, METABOLITE_NAME, VALUE) %>%
  map(pivot_wider, names_from = METABOLITE_NAME, values_from = VALUE) %>%
  map(column_to_rownames, var = "labelid_viallabel") %>%
  map(as.matrix) %>%
  map(AutoScaleMatrix) %>%
  map(as.data.frame) %>%
  map(rownames_to_column, var = "labelid_viallabel") %>%
  map(pivot_longer, names_to = 'METABOLITE_NAME', values_to = 'VALUE', cols = all_of(metabolites)) %>%
  map(mutate, REPLICATE = case_when(grepl(paste0(reps[1],"$"),labelid_viallabel) ~ reps[1],
                                     grepl(paste0(reps[2],"$"),labelid_viallabel) ~ reps[2],
                                    grepl(paste0(reps[3],"$"),labelid_viallabel) ~ reps[3])) %>%
  bind_rows()

# Plot Boxplots faceted by shared metabolites
################################################################################
norm_df %>%
  ggplot(aes(y = METABOLITE_NAME, x = VALUE, color = REPLICATE)) +
  geom_boxplot(alpha = 0.8) +
  labs(title=paste0(tissue,' ',study_institute_metab_family,": Normalized Metabolite Abundances"),
             x = "Abundance", y = "") 
```

![](/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/mayo_clinic/20201012_replicate-analysis-lung-Mayo-Clinic-amines_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Normalization-3.png)<!-- -->

``` r
# Normalize
# Subset the data for sites to compare (scatterplots)
###################################################
norm_join_df <- countdata_df %>%
  ungroup() %>%
  unite(STUDY_INSTITUTE,METAB_FAMILY, col = "STUDY_INSTITUTE_METAB_FAMILY") %>%
  filter(TISSUE == tissue) %>%
  filter(NAMED == named) %>%
  filter(STUDY_INSTITUTE_METAB_FAMILY == study_institute_metab_family) %>%
  unnest(COUNT_DATA) %>%
  filter(viallabel %in% as.character(unlist(sample_join$sample_id))) %>%
  unite(labelid, viallabel, col = "labelid_viallabel", remove = F) %>%
  mutate(REPLICATE = case_when(grepl(paste0(reps[1],"$"),labelid_viallabel) ~ reps[1],
                                     grepl(paste0(reps[2],"$"),labelid_viallabel) ~ reps[2],
                               grepl(paste0(reps[3],"$"),labelid_viallabel) ~ reps[3])) %>%
  split(.$REPLICATE) %>%
  map(select, labelid_viallabel, METABOLITE_NAME, VALUE) %>%
  map(pivot_wider, names_from = METABOLITE_NAME, values_from = VALUE) %>%
  map(column_to_rownames, var = "labelid_viallabel") %>%
  map(as.matrix) %>%
  map(AutoScaleMatrix) %>%
  map(as.data.frame) %>%
  map(rownames_to_column, var = "labelid_viallabel") %>%
  map(pivot_longer, names_to = 'METABOLITE_NAME', values_to = 'VALUE', cols = all_of(metabolites)) %>%
  map(mutate, REPLICATE = case_when(grepl(paste0(reps[1],"$"),labelid_viallabel) ~ reps[1],
                                     grepl(paste0(reps[2],"$"),labelid_viallabel) ~ reps[2],
                                    grepl(paste0(reps[3],"$"),labelid_viallabel) ~ reps[3])) %>%
  map(~rename(., !!sym(unique(.$REPLICATE)) := "VALUE")) %>%
  map(~select(., -REPLICATE)) %>%
  map(mutate, labelid = gsub(pattern = "_..*", replacement = '', labelid_viallabel)) %>%
  map(select, -labelid_viallabel) %>%
  purrr::reduce(left_join, by = c("labelid","METABOLITE_NAME"))

# Plot scatter plots
################################################################################
norm_join_df %>%
  ggplot(aes(x = !!sym(reps[1]), y = !!sym(reps[2])), color = ) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", size = 1) +
  geom_abline(linetype="dashed") +
  facet_wrap(~ METABOLITE_NAME) +
  expand_limits(x = 0, y = 0) +
  labs(title=paste0(tissue,' ',study_institute_metab_family,": Normalized Metabolite Abundances")) +
  coord_fixed(ratio=1)
```

![](/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/mayo_clinic/20201012_replicate-analysis-lung-Mayo-Clinic-amines_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Normalization-4.png)<!-- -->

``` r
# Plot scatter plots ordered by sample (include R2 values)
lm_df <- norm_join_df %>%
  rename(name = METABOLITE_NAME) %>%
  rename(x = reps[1]) %>%
  rename(y = reps[2])
# Iterate through each of the shared metabolites to collect summary info about their correlations
r2_df <- data.frame()
for(metab in metabolites){
  met_df <- lm_df %>%
    filter(name == metab)
  r2_df <- rbind(r2_df, lm_eqn(met_df))
}

# Display summaries
r2_df %>%
  arrange(METABOLITE)

r2_df %>%
  arrange(desc(R2)) %>%
  select(METABOLITE,R2)

# Plot the density plot for all the gene counts
################################################################################
norm_df %>%
  ggplot(aes(x = VALUE, color = REPLICATE)) +
  geom_density() +
  labs(title=paste0(tissue,' ',study_institute_metab_family,": Normalized Metabolite Abundances"),
       x = "Abundance",
       y = "Density")
```

![](/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/reports/mayo_clinic/20201012_replicate-analysis-lung-Mayo-Clinic-amines_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Normalization-5.png)<!-- -->

``` r
# Summary Statistics of Density plot distribution
################################################################################
# Iterate through the runs and collect vectors
df_all <- data.frame()
for(rep in reps){
  # Collect numeric vector
  num_vec <- norm_df %>%
    filter(REPLICATE == rep) %>%
    select(VALUE) %>% unlist() %>% unname()
  df <- NumericSummaryStats(num_vec)
  row.names(df) <- rep
  df_all <- rbind(df_all, df)
}
df_all
#}
```
