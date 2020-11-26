20201029\_pass1a-metabolomics-targeted-cross-tissue\_steep
================
Steep
10/29/2020

  - [Goals of Analysis](#goals-of-analysis)
  - [Visualize the Data Structure of Metabolomics
    Data](#visualize-the-data-structure-of-metabolomics-data)
      - [Setup the Environment](#setup-the-environment)
  - [Set Variables](#set-variables)
  - [Load the Metadata and Count
    Data](#load-the-metadata-and-count-data)
  - [Load & Incorporate the Phenotype
    Data](#load-incorporate-the-phenotype-data)
  - [Explore Metadata](#explore-metadata)
  - [Collect the Shared Metabolites and
    Samples](#collect-the-shared-metabolites-and-samples)
  - [Examine the Relationship between metadata variables: Try to infer
    batches](#examine-the-relationship-between-metadata-variables-try-to-infer-batches)
  - [Examine the Distribution of the Original
    Data](#examine-the-distribution-of-the-original-data)
  - [WARNING: Subsets of samples have been taken for development
    purposes](#warning-subsets-of-samples-have-been-taken-for-development-purposes)
  - [Ensure that QC samples are not used in
    scaling](#ensure-that-qc-samples-are-not-used-in-scaling)
  - [Impute Missing Values (Necessary for
    PCA)](#impute-missing-values-necessary-for-pca)
  - [Combine AutoScaled Datasets for PCA (Matrix; Columns:Metabolites,
    Rows: Unique\_ID)
    (labelid+INSTITUTE+TECH))](#combine-autoscaled-datasets-for-pca-matrix-columnsmetabolites-rows-unique_id-labelidinstitutetech)

## Goals of Analysis

# Visualize the Data Structure of Metabolomics Data

## Setup the Environment

``` r
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
#############################################################################
COMP <- "MacBook" 
#COMP <- "GreatLakes"       
if(COMP == "MacBook"){
        # MacBook
        WD <- '/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a'    
}else if(COMP == "GreatLakes"){
        # MichMed
        WD <- '/home/acsteep/motrpac/20200915_metabolomics-pass1a'
}

#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("org.Rn.eg.db")
#install.packages("kableExtra")

# Load dependencies
pacs...man <- c("tidyverse", "data.table","rlang","RColorBrewer","MASS","sfsmisc",
                "pheatmap","caret","ggforce","RANN","kableExtra")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

############################################################
##### Functions ############################################
############################################################

# Name functions
select <- dplyr::select
counts <- DESeq2::counts
map <- purrr::map
desc <- dplyr::desc
arrange <- dplyr::arrange
melt <- reshape2::melt

# Global options
options(dplyr.print_max = 100)
options(scipen=10000)

# Source the functions
source(paste0(WD,'/functions/not_in.R'))
source(paste0(WD,'/functions/rat_mouse_ortho.R'))
source(paste0(WD,'/functions/mouse2rat_ortho.R'))
source(paste0(WD,'/functions/lmp.R'))
source(paste0(WD,'/functions/cor_PC_1_6.R'))
source(paste0(WD,'/functions/elbow_finder.R'))
source(paste0(WD,'/functions/cor_outlier2.R'))
source(paste0(WD,'/functions/DissimilarityMatrixNominal.R'))
source(paste0(WD,'/functions/lm_eqn.R'))
source(paste0(WD,'/functions/reorder_cormat.R'))
source(paste0(WD,'/functions/save_pheatmap_pdf.R'))
source(paste0(WD,'/functions/modeav.R'))
source(paste0(WD,'/functions/count_na.R'))
source(paste0(WD,'/functions/se_mean.R'))
source(paste0(WD,'/functions/Mode.R'))
source(paste0(WD,'/functions/NumericSummaryStats.R'))
source(paste0(WD,'/functions/AutoScaleMatrix.R'))
source(paste0(WD,'/functions/ReverseAutoScaleMatrix.R'))
source(paste0(WD,'/functions/CenterMatrix.R'))
source(paste0(WD,'/functions/RangeScaleMatrix.R'))
source(paste0(WD,'/functions/ParetoScaleMatrix.R'))
source(paste0(WD,'/functions/VastScaleMatrix.R'))
source(paste0(WD,'/functions/LevelScaleMatrix.R'))
source(paste0(WD,'/functions/Log10Matrix.R'))
source(paste0(WD,'/functions/PowerMatrix.R'))
```

# Set Variables

``` r
################################################################################
####################### Set Variables ##########################################
################################################################################

# Dataset
###################################
dataset <- 'untargeted'

# Named
###################################
named <- 'named'

# Tissue
###################################
tissue <- 'liver'

# Metabolite family
###################################
metab_family <- 'hilicpos'

# Transformation Technique
###################################
# Original Data
Untransformed <- 'COUNT_DATA'

# Correlation
###################################
Corr <- 'Pearson'
Corr <- 'Spearman'
```

# Load the Metadata and Count Data

``` r
# Load the Data
################################################################################
# Files
metadata_rds <- paste0(WD, '/data/20201010_pass1a-metabolomics-metadata-table_steep.rds')
# Load Data
metadata_df <- readRDS(file = metadata_rds)

metadata_df %>%
  head(n =2) %>% kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

OME

</th>

<th style="text-align:left;">

DATASET

</th>

<th style="text-align:left;">

TISSUE

</th>

<th style="text-align:left;">

METAB\_FAMILY

</th>

<th style="text-align:left;">

NAMED

</th>

<th style="text-align:left;">

STUDY\_INSTITUTE

</th>

<th style="text-align:left;">

STUDY\_TITLE

</th>

<th style="text-align:left;">

STUDY\_TYPE

</th>

<th style="text-align:left;">

STUDY\_SUMMARY

</th>

<th style="text-align:left;">

STUDY\_DEPARTMENT

</th>

<th style="text-align:left;">

STUDY\_LABORATORY

</th>

<th style="text-align:left;">

STUDY\_LAST\_NAME

</th>

<th style="text-align:left;">

ST\_NUM\_GROUPS

</th>

<th style="text-align:left;">

SUBMIT\_DATE

</th>

<th style="text-align:left;">

STUDY\_COMMENTS

</th>

<th style="text-align:left;">

SUBJECT\_TYPE

</th>

<th style="text-align:left;">

SUBJECT\_SPECIES

</th>

<th style="text-align:left;">

SAMPLEPREP\_SUMMARY

</th>

<th style="text-align:left;">

SAMPLEPREP\_PROTOCOL\_FILENAME

</th>

<th style="text-align:left;">

CH\_CHROMATOGRAPHY\_TYPE

</th>

<th style="text-align:left;">

CH\_INSTRUMENT\_NAME

</th>

<th style="text-align:left;">

CH\_COLUMN\_NAME

</th>

<th style="text-align:left;">

CH\_METHODS\_FILENAME

</th>

<th style="text-align:left;">

CH\_SUMMARY

</th>

<th style="text-align:left;">

MS\_INSTRUMENT\_TYPE

</th>

<th style="text-align:left;">

MS\_INSTRUMENT\_NAME

</th>

<th style="text-align:left;">

MS\_TYPE

</th>

<th style="text-align:left;">

MS\_ION\_MODE

</th>

<th style="text-align:left;">

MS\_UNITS

</th>

<th style="text-align:left;">

MS\_COMMENTS

</th>

<th style="text-align:left;">

MS\_RESULTS\_FILE

</th>

<th style="text-align:left;">

ANALYSIS\_TYPE

</th>

<th style="text-align:left;">

ANALYSIS\_DETAILS

</th>

<th style="text-align:left;">

METABOLITE\_NAMES

</th>

<th style="text-align:right;">

METABOLITE\_N

</th>

<th style="text-align:left;">

SAMPLE\_NAMES

</th>

<th style="text-align:right;">

SAMPLE\_N

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

metabolomics

</td>

<td style="text-align:left;">

targeted

</td>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:left;">

3hib

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Duke University

</td>

<td style="text-align:left;">

PASS 1A

</td>

<td style="text-align:left;">

Animal

</td>

<td style="text-align:left;">

exercised vs. sedentary

</td>

<td style="text-align:left;">

Duke Molecular Physiology Institute

</td>

<td style="text-align:left;">

Metabolomics Core Laboratory

</td>

<td style="text-align:left;">

Ilkayeva

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

animal

</td>

<td style="text-align:left;">

rat

</td>

<td style="text-align:left;">

50 �l of plasma or 100 �l of tissue homogenate prepared at 50 mg of wet
tissue per 1 ml of homogenate in 50% acetonitrile/0.3% formic acid in
water containing an isotopically labeled internal standard
d6-2-hydroxyisobutyric acid (2-HIB) (CDN Isotopes) are precipitated with
400 �l of methanol. The methanol supernatants are dried under nitrogen
and econstituted in water.

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

RPC

</td>

<td style="text-align:left;">

Waters Acquity UPLC system

</td>

<td style="text-align:left;">

“Waters Acquity UPLC HSS T3 Column, 1.8 �m, 2.1 � 100 mm”

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

“The analytical column is used at 30�C, 10 �l of the sample are injected
onto the column, and eluted isocratically at 95% eluent A (0.1 % formic
acid in water) and 5% eluent B (acetonitrile) and a flow rate of 0.4
ml/min. The total run time is 6.5 min.”

</td>

<td style="text-align:left;">

triple quadrupole mass spectrometer

</td>

<td style="text-align:left;">

Waters Xevo TQ-S

</td>

<td style="text-align:left;">

MRM

</td>

<td style="text-align:left;">

negative

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

Mass transitions of m/z 103-\> 73 (3-hydroxyisobutyrate or 3-HIB) and
109 -\> 62 (d6-2-HIB) are monitored in a negative ion electrospray
ionization mode.

</td>

<td style="text-align:left;">

MS Results File(Required for untargeted expt)

</td>

<td style="text-align:left;">

targeted

</td>

<td style="text-align:left;">

“The ion ratio of 3-hydroxyisobutyrate to the internal standard 2-HIB-d6
is computed from extracted ion chromatograms using a software package
TargetLynx (Waters, Milford, MA). The ratios are converted to
concentrations using calibrators constructed from authentic 3-HIB
(Sigma, MO, USA) and Dialyzed Fetal Bovine Serum (Sigma, MO, USA). The
values are expressed in �M units for plasma samples or pmol/mg for
tissue samples.”

</td>

<td style="text-align:left;">

3-hydroxyisobutyrate

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

90008013107;90141013107;90047013107;90157013107;90134013107;90043013107;90028013107;90138013107;90031013107;90041013107;90110013107;90148013107;90044013107;90119013107;90027013107;90018013107;90046013107;90007013107;90133013107;90140013107;90013013107;90040013107;90050013107;90033013107;90017013107;90129013107;90042013107;90136013107;90023013107;90045013107;90029013107;90032013107;90124013107;90123013107;90126013107;90120013107;90053013107;90156013107;90117013107;90009013107;90128013107;90159013107;90037013107;90147013107;90160013107;90144013107;90005013107;90039013107;90115013107;90010013107;90025013107;90145013107;90109013107;90114013107;90020013107;90038013107;90150013107;90122013107;90152013107;90112013107;90121013107;90034013107;90026013107;90012013107;90015013107;90127013107;90001013107;90014013107;90143013107;90048013107;90155013107;90011013107;90139013107;90135013107;90052013107;90118013107;90142013107;90146013107

</td>

<td style="text-align:right;">

78

</td>

</tr>

<tr>

<td style="text-align:left;">

metabolomics

</td>

<td style="text-align:left;">

targeted

</td>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:left;">

aa

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Duke University

</td>

<td style="text-align:left;">

PASS 1A

</td>

<td style="text-align:left;">

Animal

</td>

<td style="text-align:left;">

exercised vs. sedentary

</td>

<td style="text-align:left;">

Duke Molecular Physiology Institute

</td>

<td style="text-align:left;">

Metabolomics Core Laboratory

</td>

<td style="text-align:left;">

Ilkayeva

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

animal

</td>

<td style="text-align:left;">

rat

</td>

<td style="text-align:left;">

“Amino Acids are analyzed using stable isotope dilution techniques by
flow injection tandem mass spectrometry. 100 microL of plasma or tissue
homogenate preprared at 50 mg/ml in 50% acetonitrile/0.3% formic acid
are spiked with a cocktail of heavy-isotope internal standards
(Cambridge Isotope Laboratories, MA, USA; CDN Isotopes, Canada) and
deproteinated with methanol. The methanol supernatants are dried and
esterified with acidified butanol.”

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

FIC

</td>

<td style="text-align:left;">

Waters Acquity UPLC system

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

“80% methanol is used as a mobile phase, the flow is 0.035 ml/min.”

</td>

<td style="text-align:left;">

triple quadrupole mass spectrometer

</td>

<td style="text-align:left;">

Waters Xevo TQD

</td>

<td style="text-align:left;">

neutral loss scan (MS/MS)

</td>

<td style="text-align:left;">

positive

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

“Mass spectra for amino acids are obtained using neutral loss scans of
m/z=102, 119, and 161. The spectra are acquired in a multi-channel
analyzer (MCA) mode to improve signal-to-noise. The data are generated
using a Waters Xevo TQD controlled by MassLynx 4.1 operating system
(Waters, Milford, MA).”

</td>

<td style="text-align:left;">

MS Results File(Required for untargeted expt)

</td>

<td style="text-align:left;">

targeted

</td>

<td style="text-align:left;">

“Ion ratios of amino acids to respective internal standards are computed
from centroided spectra using a software package NeoLynx (Waters,
Milford, MA). The ratios are converted to concentrations using
calibrators constructed from authentic amino acids (Sigma, MO, USA) and
Dialyzed Fetal Bovine Serum (Sigma, MO, USA). The values are expressed
in �M units for plasma and nmol/mg for tissue samples.”

</td>

<td style="text-align:left;">

alanine;arginine;aspartate;citrulline;glutamate;glycine;histidine;leucine+isoleucine;methionine;ornithine;phenylalanine;proline;serine;tyrosine;valine

</td>

<td style="text-align:right;">

15

</td>

<td style="text-align:left;">

90008013107;90141013107;90047013107;90157013107;90134013107;90043013107;90028013107;90138013107;90031013107;90041013107;90110013107;90148013107;90044013107;90119013107;90027013107;90018013107;90046013107;90007013107;90133013107;90140013107;90013013107;90040013107;90050013107;90033013107;90017013107;90129013107;90042013107;90136013107;90023013107;90045013107;90029013107;90032013107;90124013107;90123013107;90126013107;90120013107;90053013107;90156013107;90117013107;90009013107;90128013107;90159013107;90037013107;90147013107;90160013107;90144013107;90005013107;90039013107;90115013107;90010013107;90025013107;90145013107;90109013107;90114013107;90020013107;90038013107;90150013107;90122013107;90152013107;90112013107;90121013107;90034013107;90026013107;90012013107;90015013107;90127013107;90001013107;90014013107;90143013107;90048013107;90155013107;90011013107;90139013107;90135013107;90052013107;90118013107;90142013107;90146013107

</td>

<td style="text-align:right;">

78

</td>

</tr>

</tbody>

</table>

``` r
# Start the clock
ptm <- proc.time()
# Takes 21+ seconds (881.6 Mb) 1073+ seconds for (7.5Gb)
#countdata_rds <- paste0(WD, #'/data/20201010_pass1a-metabolomics-countdata-nested_steep.rds')

# Load a specific tissues count data
# Liver takes 230 seconds (4 min)
#countdata_rds <- paste0(WD, '/data/20201010_pass1a-metabolomics-countdata-nested-',tissue,'_steep.rds')

# Load a specific metabolite family count data
# ac takes 2 seconds
# hilicpos takes 453 seconds (on MB) or 260 seconds (on GL with 40 MB 2 cores)
# All data (original + autoscaled) takes 170 seconds on MB
countdata_rds <- paste0(WD, '/data/20201010_pass1a-metabolomics-countdata-nested-orig_steep.rds')

# Load the data
countdata_df <- readRDS(file = countdata_rds)

# Stop the clock
proc.time() - ptm
```

    ##    user  system elapsed 
    ##  87.493   1.209  90.426

# Load & Incorporate the Phenotype Data

``` r
# Load the phenotype data
pheno_obj <- paste0(WD,"/data/20201021_pass1a-06-pheno-viallabel_steep.rds")
pheno_df <- readRDS(pheno_obj)

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

pheno_df %>%
  head(n =2) %>% kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

pid

</th>

<th style="text-align:left;">

bid

</th>

<th style="text-align:left;">

labelid

</th>

<th style="text-align:left;">

viallabel

</th>

<th style="text-align:right;">

Key.protocol

</th>

<th style="text-align:right;">

Key.agegroup

</th>

<th style="text-align:left;">

Key.d\_arrive

</th>

<th style="text-align:left;">

Key.d\_sacrifice

</th>

<th style="text-align:right;">

Key.intervention

</th>

<th style="text-align:right;">

Key.sacrificetime

</th>

<th style="text-align:left;">

Key.anirandgroup

</th>

<th style="text-align:right;">

Key.sitename

</th>

<th style="text-align:left;">

Registration.d\_visit

</th>

<th style="text-align:right;">

Registration.days\_visit

</th>

<th style="text-align:right;">

Registration.staffid

</th>

<th style="text-align:right;">

Registration.siteid

</th>

<th style="text-align:left;">

Registration.formname

</th>

<th style="text-align:right;">

Registration.versionnbr

</th>

<th style="text-align:right;">

Registration.ratid

</th>

<th style="text-align:left;">

Registration.d\_arrive

</th>

<th style="text-align:right;">

Registration.days\_arrive

</th>

<th style="text-align:right;">

Registration.batchnumber

</th>

<th style="text-align:left;">

Registration.d\_reverselight

</th>

<th style="text-align:right;">

Registration.days\_reverselight

</th>

<th style="text-align:left;">

Registration.d\_birth

</th>

<th style="text-align:right;">

Registration.days\_birth

</th>

<th style="text-align:right;">

Registration.sex

</th>

<th style="text-align:right;">

Registration.weight

</th>

<th style="text-align:left;">

Registration.cagenumber

</th>

<th style="text-align:left;">

Registration.comments

</th>

<th style="text-align:left;">

Familiarization.d\_visit

</th>

<th style="text-align:right;">

Familiarization.days\_visit

</th>

<th style="text-align:right;">

Familiarization.staffid

</th>

<th style="text-align:right;">

Familiarization.siteid

</th>

<th style="text-align:left;">

Familiarization.formname

</th>

<th style="text-align:right;">

Familiarization.versionnbr

</th>

<th style="text-align:left;">

Familiarization.d\_treadmillbegin

</th>

<th style="text-align:right;">

Familiarization.days\_treadmillbegin

</th>

<th style="text-align:left;">

Familiarization.d\_treadmillcomplete

</th>

<th style="text-align:right;">

Familiarization.days\_treadmillcomplete

</th>

<th style="text-align:right;">

Familiarization.activity\_score

</th>

<th style="text-align:right;">

Familiarization.compliant

</th>

<th style="text-align:right;">

Familiarization.weight

</th>

<th style="text-align:left;">

Familiarization.fat

</th>

<th style="text-align:left;">

Familiarization.lean

</th>

<th style="text-align:left;">

Familiarization.comments

</th>

<th style="text-align:left;">

Acute.Test.d\_visit

</th>

<th style="text-align:right;">

Acute.Test.days\_visit

</th>

<th style="text-align:right;">

Acute.Test.staffid

</th>

<th style="text-align:right;">

Acute.Test.siteid

</th>

<th style="text-align:left;">

Acute.Test.formname

</th>

<th style="text-align:right;">

Acute.Test.versionnbr

</th>

<th style="text-align:left;">

Acute.Test.d\_start

</th>

<th style="text-align:right;">

Acute.Test.days\_start

</th>

<th style="text-align:left;">

Acute.Test.t\_start

</th>

<th style="text-align:right;">

Acute.Test.weight

</th>

<th style="text-align:right;">

Acute.Test.beginblood

</th>

<th style="text-align:right;">

Acute.Test.distance

</th>

<th style="text-align:right;">

Acute.Test.speed

</th>

<th style="text-align:right;">

Acute.Test.incline

</th>

<th style="text-align:left;">

Acute.Test.t\_complete

</th>

<th style="text-align:right;">

Acute.Test.endblood

</th>

<th style="text-align:right;">

Acute.Test.contactshock

</th>

<th style="text-align:left;">

Acute.Test.howlongshock

</th>

<th style="text-align:right;">

Acute.Test.timesshock

</th>

<th style="text-align:right;">

Acute.Test.intenseinitial

</th>

<th style="text-align:right;">

Acute.Test.intensefinal

</th>

<th style="text-align:left;">

Acute.Test.comments

</th>

<th style="text-align:left;">

Specimen.Collection.d\_visit

</th>

<th style="text-align:right;">

Specimen.Collection.days\_visit

</th>

<th style="text-align:right;">

Specimen.Collection.staffid

</th>

<th style="text-align:right;">

Specimen.Collection.siteid

</th>

<th style="text-align:left;">

Specimen.Collection.formname

</th>

<th style="text-align:right;">

Specimen.Collection.versionnbr

</th>

<th style="text-align:right;">

Specimen.Collection.anesthesiaid

</th>

<th style="text-align:left;">

Specimen.Collection.t\_anesthesia

</th>

<th style="text-align:left;">

Specimen.Collection.anesthesiacomments

</th>

<th style="text-align:left;">

Specimen.Collection.bloodtype

</th>

<th style="text-align:left;">

Specimen.Collection.bloodtube

</th>

<th style="text-align:right;">

Specimen.Collection.bloodcomplete

</th>

<th style="text-align:right;">

Specimen.Collection.bloodtechid

</th>

<th style="text-align:left;">

Specimen.Collection.t\_bloodstart

</th>

<th style="text-align:left;">

Specimen.Collection.t\_bloodstop

</th>

<th style="text-align:left;">

Specimen.Collection.t\_edtafill

</th>

<th style="text-align:left;">

Specimen.Collection.bloodcomments

</th>

<th style="text-align:left;">

Specimen.Collection.uterustype

</th>

<th style="text-align:left;">

Specimen.Collection.uterustechid

</th>

<th style="text-align:left;">

Specimen.Collection.uteruscomplete

</th>

<th style="text-align:left;">

Specimen.Collection.t\_uterusstart

</th>

<th style="text-align:left;">

Specimen.Collection.t\_uterusstop

</th>

<th style="text-align:left;">

Specimen.Collection.uterusweight

</th>

<th style="text-align:left;">

Specimen.Collection.uteruscomments

</th>

<th style="text-align:left;">

Specimen.Collection.t\_death

</th>

<th style="text-align:right;">

Specimen.Collection.deathtype

</th>

<th style="text-align:right;">

Specimen.Processing.versionnbr

</th>

<th style="text-align:left;">

Specimen.Processing.formname

</th>

<th style="text-align:right;">

Specimen.Processing.siteid

</th>

<th style="text-align:left;">

Specimen.Processing.sampletypedescription

</th>

<th style="text-align:right;">

Specimen.Processing.samplenumber

</th>

<th style="text-align:right;">

Specimen.Processing.timepoint

</th>

<th style="text-align:left;">

Specimen.Processing.aliquotdescription

</th>

<th style="text-align:left;">

Specimen.Processing.volume

</th>

<th style="text-align:left;">

Specimen.Processing.partialamt

</th>

<th style="text-align:right;">

Specimen.Processing.hemolyzed

</th>

<th style="text-align:left;">

Specimen.Processing.comments

</th>

<th style="text-align:left;">

Specimen.Processing.t\_collection

</th>

<th style="text-align:left;">

Specimen.Processing.t\_edtaspin

</th>

<th style="text-align:left;">

Specimen.Processing.t\_freeze

</th>

<th style="text-align:right;">

Specimen.Processing.techid

</th>

<th style="text-align:right;">

Calculated.Variables.wgt\_gain\_before\_acute

</th>

<th style="text-align:right;">

Calculated.Variables.lactate\_change\_dueto\_acute

</th>

<th style="text-align:right;">

Calculated.Variables.edta\_coll\_time

</th>

<th style="text-align:right;">

Calculated.Variables.deathtime\_after\_acute

</th>

<th style="text-align:right;">

Calculated.Variables.frozetime\_after\_acute

</th>

<th style="text-align:right;">

BICLabelData.shiptositeid

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

10029028

</td>

<td style="text-align:left;">

90001

</td>

<td style="text-align:left;">

90001013001

</td>

<td style="text-align:left;">

90001013001

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

Control - IPE

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

19JUN2018

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Registration

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

Jan 2018

</td>

<td style="text-align:right;">

\-120

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

345

</td>

<td style="text-align:left;">

3M

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

13JUL2018

</td>

<td style="text-align:right;">

73

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Familiarization

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

03JUL2018

</td>

<td style="text-align:right;">

63

</td>

<td style="text-align:left;">

13JUL2018

</td>

<td style="text-align:right;">

73

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

358

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Acute Test

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:left;">

09:20:00

</td>

<td style="text-align:right;">

356

</td>

<td style="text-align:right;">

1.3

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

09:50:00

</td>

<td style="text-align:right;">

1.6

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

1407

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Specimen Collection

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1295

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

Administration time missing. Animal woke up from under anesthesia, death
by cervical dislocation and decapitation.

</td>

<td style="text-align:left;">

EDTA

</td>

<td style="text-align:left;">

Draw

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1295

</td>

<td style="text-align:left;">

09:55:43

</td>

<td style="text-align:left;">

09:56:45

</td>

<td style="text-align:left;">

09:56:50

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

09:57:15

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

Animal Sample Processing

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

PaxGene RNA

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

1.1mL PaxGene RNA

</td>

<td style="text-align:left;">

Full

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

09:57:15

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

11:57:15

</td>

<td style="text-align:right;">

1407

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

0.3

</td>

<td style="text-align:right;">

410

</td>

<td style="text-align:right;">

435

</td>

<td style="text-align:right;">

7635

</td>

<td style="text-align:right;">

850

</td>

</tr>

<tr>

<td style="text-align:left;">

10029028

</td>

<td style="text-align:left;">

90001

</td>

<td style="text-align:left;">

90001013101

</td>

<td style="text-align:left;">

90001013104

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

Control - IPE

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

19JUN2018

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Registration

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

01MAY2018

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

Jan 2018

</td>

<td style="text-align:right;">

\-120

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

345

</td>

<td style="text-align:left;">

3M

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

13JUL2018

</td>

<td style="text-align:right;">

73

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Familiarization

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

03JUL2018

</td>

<td style="text-align:right;">

63

</td>

<td style="text-align:left;">

13JUL2018

</td>

<td style="text-align:right;">

73

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

358

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

1293

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Acute Test

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:left;">

09:20:00

</td>

<td style="text-align:right;">

356

</td>

<td style="text-align:right;">

1.3

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

09:50:00

</td>

<td style="text-align:right;">

1.6

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

18JUL2018

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

1407

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

Animal Specimen Collection

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1295

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

Administration time missing. Animal woke up from under anesthesia, death
by cervical dislocation and decapitation.

</td>

<td style="text-align:left;">

EDTA

</td>

<td style="text-align:left;">

Draw

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1295

</td>

<td style="text-align:left;">

09:55:43

</td>

<td style="text-align:left;">

09:56:45

</td>

<td style="text-align:left;">

09:56:50

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

09:57:15

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

Animal Sample Processing

</td>

<td style="text-align:right;">

910

</td>

<td style="text-align:left;">

EDTA Plasma

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

0.5 mL Plasma

</td>

<td style="text-align:left;">

Full

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

10:18:19

</td>

<td style="text-align:left;">

10:43:30

</td>

<td style="text-align:right;">

1406

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

0.3

</td>

<td style="text-align:right;">

410

</td>

<td style="text-align:right;">

435

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

870

</td>

</tr>

</tbody>

</table>

# Explore Metadata

``` r
# First cut is either by tissue or metabolite family/tech
metadata_df %>% 
  select(TISSUE,METAB_FAMILY) %>% table() %>% as.data.frame.matrix() %>%
  select("hilicpos","rppos","rpneg","lrppos","lrpneg","ionpneg","cer","ac",
         "oxylipneg","tca","amines","sphm","ka","aa","3hib","nuc","acoa","oa","baiba") %>%
  arrange(desc("hilicpos"),desc("rppos"),desc("rpneg"),desc("lrppos"),desc("lrpneg"),
          desc("ionpneg"),desc("cer"),desc("ac"),desc("oxylipneg"),desc("tca"),
          desc("amines"),desc("sphm"),desc("ka"),desc("aa"),desc("3hib"),desc("nuc"),
          desc("acoa"),desc("oa"),desc("baiba")) %>% kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

hilicpos

</th>

<th style="text-align:right;">

rppos

</th>

<th style="text-align:right;">

rpneg

</th>

<th style="text-align:right;">

lrppos

</th>

<th style="text-align:right;">

lrpneg

</th>

<th style="text-align:right;">

ionpneg

</th>

<th style="text-align:right;">

cer

</th>

<th style="text-align:right;">

ac

</th>

<th style="text-align:right;">

oxylipneg

</th>

<th style="text-align:right;">

tca

</th>

<th style="text-align:right;">

amines

</th>

<th style="text-align:right;">

sphm

</th>

<th style="text-align:right;">

ka

</th>

<th style="text-align:right;">

aa

</th>

<th style="text-align:right;">

3hib

</th>

<th style="text-align:right;">

nuc

</th>

<th style="text-align:right;">

acoa

</th>

<th style="text-align:right;">

oa

</th>

<th style="text-align:right;">

baiba

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

hippocampus

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

gastrocnemius

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

heart

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

kidney

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

lung

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

liver

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

brown-adipose

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

white-adipose

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

cortex

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

hypothalmus

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

vastus-lateralis

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

tibia

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

adrenal

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

colon

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

spleen

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

testes

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

ovaries

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

aorta

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

small-intestine

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

``` r
metadata_df %>% 
  filter(NAMED == 'named') %>%
  filter(METAB_FAMILY == metab_family) %>%
        select(DATASET,TISSUE,METAB_FAMILY,NAMED,STUDY_INSTITUTE,
               CH_CHROMATOGRAPHY_TYPE,MS_TYPE,MS_ION_MODE,DATASET) %>% 
        unite("CHROM_MS_ION", CH_CHROMATOGRAPHY_TYPE:MS_TYPE:MS_ION_MODE, remove = T) %>%
  kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

DATASET

</th>

<th style="text-align:left;">

TISSUE

</th>

<th style="text-align:left;">

METAB\_FAMILY

</th>

<th style="text-align:left;">

NAMED

</th>

<th style="text-align:left;">

STUDY\_INSTITUTE

</th>

<th style="text-align:left;">

CHROM\_MS\_ION

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

hippocampus

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

cortex

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

hypothalmus

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

gastrocnemius

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

vastus-lateralis

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

tibia

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

heart

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

kidney

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

adrenal

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

colon

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

spleen

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

testes

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

ovaries

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

aorta

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

lung

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

small-intestine

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

liver

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

brown-adipose

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

white-adipose

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

</tbody>

</table>

# Collect the Shared Metabolites and Samples

``` r
# There are 20 tissues
# Shared metabolites across all 20 tissues
all_tis_mets <- countdata_df %>%
  filter(METAB_FAMILY == metab_family) %>%
  filter(NAMED == named) %>%
  select(-COUNT_DATA, -'SAMPLE_DATA') %>%
  unnest(cols = 'METABOLITE_DATA') %>%
  filter(!is.na(refmet_name)) %>%
  group_by(refmet_name) %>%
  mutate(SHARED_METABOLITE_N = n()) %>%
  arrange(desc(SHARED_METABOLITE_N)) %>%
  ungroup() %>%
  filter(SHARED_METABOLITE_N == 20) %>%
  select(refmet_name) %>% unlist() %>% as.character() %>% unique()
  
all_tis_mets %>% kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

phenylalanine-d8 \[iSTD\]

</td>

</tr>

<tr>

<td style="text-align:left;">

valine-d8 \[iSTD\]

</td>

</tr>

<tr>

<td style="text-align:left;">

Deoxycytidine

</td>

</tr>

<tr>

<td style="text-align:left;">

Adenosine

</td>

</tr>

<tr>

<td style="text-align:left;">

Alanine

</td>

</tr>

<tr>

<td style="text-align:left;">

sn-Glycero-3-phosphocholine

</td>

</tr>

<tr>

<td style="text-align:left;">

Arginine

</td>

</tr>

<tr>

<td style="text-align:left;">

Asparagine

</td>

</tr>

<tr>

<td style="text-align:left;">

Betaine

</td>

</tr>

<tr>

<td style="text-align:left;">

Biliverdin

</td>

</tr>

<tr>

<td style="text-align:left;">

3-Dehydroxycarnitine

</td>

</tr>

<tr>

<td style="text-align:left;">

Carnitine

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(2:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(3:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(4:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(4:0(OH))

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(5:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(6:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(8:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(10:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(12:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(12:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(14:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(14:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(14:2)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(16:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(16:0(OH))

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(18:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(18:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(18:2)

</td>

</tr>

<tr>

<td style="text-align:left;">

CAR(20:4)

</td>

</tr>

<tr>

<td style="text-align:left;">

Choline

</td>

</tr>

<tr>

<td style="text-align:left;">

Citrulline

</td>

</tr>

<tr>

<td style="text-align:left;">

Creatine

</td>

</tr>

<tr>

<td style="text-align:left;">

Creatinine

</td>

</tr>

<tr>

<td style="text-align:left;">

Cyclohexylamine

</td>

</tr>

<tr>

<td style="text-align:left;">

Cytidine

</td>

</tr>

<tr>

<td style="text-align:left;">

Ectoine

</td>

</tr>

<tr>

<td style="text-align:left;">

Glutamic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Glutamine

</td>

</tr>

<tr>

<td style="text-align:left;">

Glycine

</td>

</tr>

<tr>

<td style="text-align:left;">

Guanidoacetic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Histamine

</td>

</tr>

<tr>

<td style="text-align:left;">

Histidine

</td>

</tr>

<tr>

<td style="text-align:left;">

Hydroxyproline

</td>

</tr>

<tr>

<td style="text-align:left;">

Hypotaurine

</td>

</tr>

<tr>

<td style="text-align:left;">

Imidazolepropionic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Inosine

</td>

</tr>

<tr>

<td style="text-align:left;">

Isoleucine

</td>

</tr>

<tr>

<td style="text-align:left;">

Leucine

</td>

</tr>

<tr>

<td style="text-align:left;">

Linoleoyl-EA

</td>

</tr>

<tr>

<td style="text-align:left;">

Lysine

</td>

</tr>

<tr>

<td style="text-align:left;">

Methionine

</td>

</tr>

<tr>

<td style="text-align:left;">

5’-Methylthioadenosine

</td>

</tr>

<tr>

<td style="text-align:left;">

N1-Acetylspermidine

</td>

</tr>

<tr>

<td style="text-align:left;">

N(6),N(6)-Dimethyl-lysine

</td>

</tr>

<tr>

<td style="text-align:left;">

N-6-Trimethyllysine

</td>

</tr>

<tr>

<td style="text-align:left;">

N6-Acetyllysine

</td>

</tr>

<tr>

<td style="text-align:left;">

N-Acetylornithine

</td>

</tr>

<tr>

<td style="text-align:left;">

Niacinamide

</td>

</tr>

<tr>

<td style="text-align:left;">

N-Lauroylglycine

</td>

</tr>

<tr>

<td style="text-align:left;">

Pantothenic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Phenylalanine

</td>

</tr>

<tr>

<td style="text-align:left;">

Choline phosphate

</td>

</tr>

<tr>

<td style="text-align:left;">

Pipecolic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Proline

</td>

</tr>

<tr>

<td style="text-align:left;">

Proline betaine

</td>

</tr>

<tr>

<td style="text-align:left;">

Putrescine

</td>

</tr>

<tr>

<td style="text-align:left;">

Serine

</td>

</tr>

<tr>

<td style="text-align:left;">

Serotonin

</td>

</tr>

<tr>

<td style="text-align:left;">

S-Methylcysteine S-oxide

</td>

</tr>

<tr>

<td style="text-align:left;">

Taurine

</td>

</tr>

<tr>

<td style="text-align:left;">

Thiamine

</td>

</tr>

<tr>

<td style="text-align:left;">

Threonine

</td>

</tr>

<tr>

<td style="text-align:left;">

Trigonelline

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimethylamine N-oxide

</td>

</tr>

<tr>

<td style="text-align:left;">

Tryptophan

</td>

</tr>

<tr>

<td style="text-align:left;">

Tyrosine

</td>

</tr>

<tr>

<td style="text-align:left;">

Urocanic acid

</td>

</tr>

<tr>

<td style="text-align:left;">

Valine

</td>

</tr>

<tr>

<td style="text-align:left;">

LPC(16:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPC(0:0/16:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPC(18:3)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPC(18:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPE(16:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPE(18:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

LPE(20:4)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(30:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(34:3)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(P-34:2)/PC(O-34:3)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(P-36:4)/PC(O-36:5)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(P-38:6)/PC(O-38:7)

</td>

</tr>

<tr>

<td style="text-align:left;">

PC(P-38:5)/PC(O-38:6)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(34:2)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(34:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(36:4)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(36:2)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(38:6)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(38:4)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(40:6)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-34:2)/PE(O-34:3)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-34:1)/PE(O-34:2)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-36:4)/PE(O-36:5)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-38:6)/PE(O-38:7)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-38:4)/PE(O-38:5)

</td>

</tr>

<tr>

<td style="text-align:left;">

PE(P-40:6)/PE(O-40:7)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/16:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/16:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/18:1)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/18:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/20:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

SM(d18:1/22:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

Sphingosine

</td>

</tr>

<tr>

<td style="text-align:left;">

Sphinganine

</td>

</tr>

<tr>

<td style="text-align:left;">

Cer(d18:1/16:0)

</td>

</tr>

<tr>

<td style="text-align:left;">

Cer(d18:0/24:1)

</td>

</tr>

</tbody>

</table>

``` r
# Create a metadata table elongated to labelid (just the major metadata variables) from the count data df
meta_labelid_df <- countdata_df %>%
  filter(METAB_FAMILY == metab_family) %>%
  filter(NAMED == named) %>%
  select(-COUNT_DATA, -'METABOLITE_DATA') %>%
  unnest(cols = 'SAMPLE_DATA') %>%
  mutate(viallabel = sample_id) %>%
  select(-raw_file, -sample_id) %>%
  select(TISSUE, sample_type, viallabel) %>%
  unique()
```

# Examine the Relationship between metadata variables: Try to infer batches

``` r
metadata_df %>%
filter(METAB_FAMILY == metab_family) %>%
  filter(NAMED == named) %>%
  select(DATASET,NAMED,TISSUE,METAB_FAMILY,STUDY_INSTITUTE,CH_CHROMATOGRAPHY_TYPE,MS_TYPE,MS_ION_MODE) %>%
  unite(CH_CHROMATOGRAPHY_TYPE,MS_TYPE,MS_ION_MODE, col = "TECHNOLOGY") %>%
  arrange(TISSUE,STUDY_INSTITUTE,TECHNOLOGY) %>% kbl() %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

DATASET

</th>

<th style="text-align:left;">

NAMED

</th>

<th style="text-align:left;">

TISSUE

</th>

<th style="text-align:left;">

METAB\_FAMILY

</th>

<th style="text-align:left;">

STUDY\_INSTITUTE

</th>

<th style="text-align:left;">

TECHNOLOGY

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

hippocampus

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

gastrocnemius

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

heart

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

kidney

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

lung

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

liver

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

brown-adipose

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

white-adipose

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

cortex

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

hypothalmus

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

vastus-lateralis

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

tibia

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

adrenal

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

colon

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

spleen

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

testes

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

ovaries

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

aorta

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

<tr>

<td style="text-align:left;">

untargeted

</td>

<td style="text-align:left;">

named

</td>

<td style="text-align:left;">

small-intestine

</td>

<td style="text-align:left;">

hilicpos

</td>

<td style="text-align:left;">

Broad Institute

</td>

<td style="text-align:left;">

HILIC\_ESI\_positive

</td>

</tr>

</tbody>

</table>

# Examine the Distribution of the Original Data

# WARNING: Subsets of samples have been taken for development purposes

# Ensure that QC samples are not used in scaling

``` r
# Examine the Distributions for original data and transformed data
for(TFORM in c('Untransformed','AutoScaled')){
  #Starting with matrix: convert to long format dataframe for plotting purposes
  # When transforming, ensure that QC samples are not present
  if(TFORM == 'Untransformed'){
    tissue_long_df <- countdata_df %>%
    ungroup() %>%
    filter(METAB_FAMILY == metab_family) %>%
    filter(NAMED == named) %>%
    select(TISSUE, COUNT_DATA) %>%
    unnest(COUNT_DATA) %>%
    filter(METABOLITE_NAME %in% all_tis_mets) %>%
    filter(!grepl('QC', viallabel) ) %>%
    select(TISSUE,viallabel, METABOLITE_NAME,VALUE)
    
  }else if(TFORM == 'AutoScaled'){
    tissue_long_int <- countdata_df %>%
    ungroup() %>%
    filter(METAB_FAMILY == metab_family) %>%
    filter(NAMED == named) %>%
    select(TISSUE, COUNT_DATA) %>%
    unnest(COUNT_DATA) %>%
    filter(METABOLITE_NAME %in% all_tis_mets) %>%
    filter(!grepl('QC', viallabel) ) %>%
    select(viallabel, METABOLITE_NAME,VALUE) %>%
    pivot_wider(names_from = METABOLITE_NAME, values_from = VALUE) %>%
    tibble::column_to_rownames(var = "viallabel") %>%
    as.matrix() %>% 
    AutoScaleMatrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "viallabel") %>%
    pivot_longer(names_to = "METABOLITE_NAME", values_to = "VALUE", cols = all_of(all_tis_mets))
    
    # Join the tissues
    tissue_join <- countdata_df %>%
      ungroup() %>%
      filter(METAB_FAMILY == metab_family) %>%
      filter(NAMED == named) %>%
      select(TISSUE, COUNT_DATA) %>%
      unnest(COUNT_DATA) %>%
      filter(METABOLITE_NAME %in% all_tis_mets) %>%
      filter(!grepl('QC', viallabel) ) %>%
      select(viallabel,TISSUE) %>%
      unique()
    tissue_long_df <- left_join(tissue_long_int, tissue_join, by = c('viallabel'))
  }
  # Plot Boxplots for all shared metabolites faceted by tissue
  ################################################################################
  tissue_long_df$TISSUE <- factor(tissue_long_df$TISSUE)
  p <- tissue_long_df %>%
    slice_sample(prop = 0.1) %>%
    ggplot(aes(y = METABOLITE_NAME, x = VALUE)) +
    geom_boxplot(aes(color = TISSUE), alpha = 0.5) +
    labs(title="Box Plots of Metabolite Abundances Across Tissues",
               x = "Abundance", y = "") +
    theme(axis.text.y=element_text(size=4)) +
    theme(axis.text.x=element_text(size=8)) +
    theme(legend.text = element_text(size=8)) +
    facet_wrap(~ TISSUE)
  print(p)
# Plot the density plot for all the gene counts
################################################################################
  if(TFORM == 'Untransformed'){
    p <- tissue_long_df %>%
      slice_sample(prop = 0.1) %>%
      ggplot(aes(x = VALUE, color = TISSUE)) +
      geom_density() +
      labs(title="Density Plot of Metabolite Abundances Across Tissues",
           x = "Abundance",y = "Density") +
      facet_wrap(~ TISSUE) +
      xlim(0,1000000)
    }else if(TFORM == 'AutoScaled'){
      p <- tissue_long_df %>%
        slice_sample(prop = 0.1) %>%
        ggplot(aes(x = VALUE, color = TISSUE)) +
        geom_density() +
        labs(title="Density Plot of Metabolite Abundances Across Tissues",
             x = "Abundance",y = "Density") +
        facet_wrap(~ TISSUE)
      }
    print(p)

    # Summary Statistics of Density Plot Distribution Across Tissues
    ################################################################################
    tissues <- tissue_long_df$TISSUE %>% as.character() %>% unique()
    # Iterate through the cas sites and collect vectors
    #cas_site <- 'duke'
    df_all <- data.frame()
    #tt <- tissues[1]
    for(tt in tissues){
      # Collect numeric vector
      num_vec <- tissue_long_df %>%
        filter(TISSUE == tt) %>%
        select(VALUE) %>% unlist() %>% unname()
      df <- NumericSummaryStats(num_vec)
      row.names(df) <- tt
      df_all <- rbind(df_all, df)
    }
    df_all %>%
      arrange(MEDIAN) %>% kbl() %>% kable_styling() %>% print()
}
```

![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Distribution-1.png)<!-- -->![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Distribution-2.png)<!-- -->

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

NA\_COUNT

</th>

<th style="text-align:right;">

NA\_FREQ

</th>

<th style="text-align:right;">

MEAN

</th>

<th style="text-align:right;">

MEDIAN

</th>

<th style="text-align:right;">

MAX

</th>

<th style="text-align:right;">

MIN

</th>

<th style="text-align:right;">

MID\_RANGE

</th>

<th style="text-align:right;">

VARIANCE

</th>

<th style="text-align:right;">

STD\_DEV

</th>

<th style="text-align:right;">

SE\_MEAN

</th>

<th style="text-align:right;">

Q1

</th>

<th style="text-align:right;">

Q3

</th>

<th style="text-align:right;">

KURTOSIS

</th>

<th style="text-align:right;">

SKEW

</th>

<th style="text-align:left;">

BC\_LAMBDA

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

tibia

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0005388

</td>

<td style="text-align:right;">

2187807

</td>

<td style="text-align:right;">

256363.0

</td>

<td style="text-align:right;">

36714679

</td>

<td style="text-align:right;">

1404

</td>

<td style="text-align:right;">

36713275

</td>

<td style="text-align:right;">

20721207375584

</td>

<td style="text-align:right;">

4552055

</td>

<td style="text-align:right;">

105690.4

</td>

<td style="text-align:right;">

42740.5

</td>

<td style="text-align:right;">

2011147

</td>

<td style="text-align:right;">

11.85

</td>

<td style="text-align:right;">

3.25

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

cortex

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

0.0003316

</td>

<td style="text-align:right;">

14485024

</td>

<td style="text-align:right;">

670277.0

</td>

<td style="text-align:right;">

1218821839

</td>

<td style="text-align:right;">

543

</td>

<td style="text-align:right;">

1218821296

</td>

<td style="text-align:right;">

4620459002146344

</td>

<td style="text-align:right;">

67973958

</td>

<td style="text-align:right;">

714723.8

</td>

<td style="text-align:right;">

83622.0

</td>

<td style="text-align:right;">

5314712

</td>

<td style="text-align:right;">

117.89

</td>

<td style="text-align:right;">

10.13

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

hypothalmus

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0.0003079

</td>

<td style="text-align:right;">

17488744

</td>

<td style="text-align:right;">

751069.5

</td>

<td style="text-align:right;">

1969956397

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

1969956336

</td>

<td style="text-align:right;">

11447891579210928

</td>

<td style="text-align:right;">

106994820

</td>

<td style="text-align:right;">

1327720.4

</td>

<td style="text-align:right;">

115901.0

</td>

<td style="text-align:right;">

5938063

</td>

<td style="text-align:right;">

145.57

</td>

<td style="text-align:right;">

11.55

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

0.0005526

</td>

<td style="text-align:right;">

41057857

</td>

<td style="text-align:right;">

754593.0

</td>

<td style="text-align:right;">

2213409660

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:right;">

2213409574

</td>

<td style="text-align:right;">

26269405914514300

</td>

<td style="text-align:right;">

162078394

</td>

<td style="text-align:right;">

1704389.5

</td>

<td style="text-align:right;">

131164.0

</td>

<td style="text-align:right;">

4112904

</td>

<td style="text-align:right;">

33.15

</td>

<td style="text-align:right;">

5.33

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

white-adipose

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

8760822

</td>

<td style="text-align:right;">

813384.0

</td>

<td style="text-align:right;">

387625426

</td>

<td style="text-align:right;">

120

</td>

<td style="text-align:right;">

387625306

</td>

<td style="text-align:right;">

577091954242424

</td>

<td style="text-align:right;">

24022738

</td>

<td style="text-align:right;">

252563.3

</td>

<td style="text-align:right;">

130632.5

</td>

<td style="text-align:right;">

5539146

</td>

<td style="text-align:right;">

36.51

</td>

<td style="text-align:right;">

5.26

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

aorta

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

12043869

</td>

<td style="text-align:right;">

838226.0

</td>

<td style="text-align:right;">

641767540

</td>

<td style="text-align:right;">

1846

</td>

<td style="text-align:right;">

641765694

</td>

<td style="text-align:right;">

1783047030759710

</td>

<td style="text-align:right;">

42226142

</td>

<td style="text-align:right;">

446793.7

</td>

<td style="text-align:right;">

135029.2

</td>

<td style="text-align:right;">

5682594

</td>

<td style="text-align:right;">

66.77

</td>

<td style="text-align:right;">

7.25

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

testes

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

88942732

</td>

<td style="text-align:right;">

896425.5

</td>

<td style="text-align:right;">

13784709527

</td>

<td style="text-align:right;">

1782

</td>

<td style="text-align:right;">

13784707745

</td>

<td style="text-align:right;">

748701445207312128

</td>

<td style="text-align:right;">

865275358

</td>

<td style="text-align:right;">

12864503.7

</td>

<td style="text-align:right;">

264302.8

</td>

<td style="text-align:right;">

4096090

</td>

<td style="text-align:right;">

129.27

</td>

<td style="text-align:right;">

11.29

</td>

<td style="text-align:left;">

\-0.0999999999999999

</td>

</tr>

<tr>

<td style="text-align:left;">

hippocampus

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

24773032

</td>

<td style="text-align:right;">

958025.0

</td>

<td style="text-align:right;">

2107703008

</td>

<td style="text-align:right;">

1102

</td>

<td style="text-align:right;">

2107701906

</td>

<td style="text-align:right;">

13418532781175000

</td>

<td style="text-align:right;">

115838391

</td>

<td style="text-align:right;">

1217868.0

</td>

<td style="text-align:right;">

106608.0

</td>

<td style="text-align:right;">

8078399

</td>

<td style="text-align:right;">

105.06

</td>

<td style="text-align:right;">

9.67

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

vastus-lateralis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001120

</td>

<td style="text-align:right;">

16059015

</td>

<td style="text-align:right;">

1036148.0

</td>

<td style="text-align:right;">

391480236

</td>

<td style="text-align:right;">

1078

</td>

<td style="text-align:right;">

391479158

</td>

<td style="text-align:right;">

1834072929160206

</td>

<td style="text-align:right;">

42826078

</td>

<td style="text-align:right;">

453167.0

</td>

<td style="text-align:right;">

164709.0

</td>

<td style="text-align:right;">

5595453

</td>

<td style="text-align:right;">

15.84

</td>

<td style="text-align:right;">

3.77

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

lung

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

15707896

</td>

<td style="text-align:right;">

1242464.0

</td>

<td style="text-align:right;">

754252253

</td>

<td style="text-align:right;">

1783

</td>

<td style="text-align:right;">

754250470

</td>

<td style="text-align:right;">

2441250615189642

</td>

<td style="text-align:right;">

49409014

</td>

<td style="text-align:right;">

519462.1

</td>

<td style="text-align:right;">

231774.0

</td>

<td style="text-align:right;">

5938368

</td>

<td style="text-align:right;">

45.88

</td>

<td style="text-align:right;">

6.05

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

gastrocnemius

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

18770911

</td>

<td style="text-align:right;">

1262529.0

</td>

<td style="text-align:right;">

609880603

</td>

<td style="text-align:right;">

789

</td>

<td style="text-align:right;">

609879814

</td>

<td style="text-align:right;">

2922465702951464

</td>

<td style="text-align:right;">

54059834

</td>

<td style="text-align:right;">

568358.6

</td>

<td style="text-align:right;">

184358.0

</td>

<td style="text-align:right;">

7471157

</td>

<td style="text-align:right;">

26.19

</td>

<td style="text-align:right;">

4.68

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

brown-adipose

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

21646709

</td>

<td style="text-align:right;">

1550356.5

</td>

<td style="text-align:right;">

1610797842

</td>

<td style="text-align:right;">

220

</td>

<td style="text-align:right;">

1610797622

</td>

<td style="text-align:right;">

5464572948536833

</td>

<td style="text-align:right;">

73922750

</td>

<td style="text-align:right;">

777144.6

</td>

<td style="text-align:right;">

276187.2

</td>

<td style="text-align:right;">

13331682

</td>

<td style="text-align:right;">

166.86

</td>

<td style="text-align:right;">

10.79

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

ovaries

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

0.0044209

</td>

<td style="text-align:right;">

25367871

</td>

<td style="text-align:right;">

1569144.5

</td>

<td style="text-align:right;">

2066806092

</td>

<td style="text-align:right;">

1789

</td>

<td style="text-align:right;">

2066804303

</td>

<td style="text-align:right;">

13999569677997834

</td>

<td style="text-align:right;">

118319777

</td>

<td style="text-align:right;">

1763023.7

</td>

<td style="text-align:right;">

275202.2

</td>

<td style="text-align:right;">

6875162

</td>

<td style="text-align:right;">

97.08

</td>

<td style="text-align:right;">

8.93

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

adrenal

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

54835444

</td>

<td style="text-align:right;">

1609857.5

</td>

<td style="text-align:right;">

7941380778

</td>

<td style="text-align:right;">

2175

</td>

<td style="text-align:right;">

7941378603

</td>

<td style="text-align:right;">

115196580367340992

</td>

<td style="text-align:right;">

339406217

</td>

<td style="text-align:right;">

3568153.2

</td>

<td style="text-align:right;">

251052.8

</td>

<td style="text-align:right;">

6548580

</td>

<td style="text-align:right;">

131.52

</td>

<td style="text-align:right;">

10.10

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

kidney

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

28490863

</td>

<td style="text-align:right;">

1649135.0

</td>

<td style="text-align:right;">

2589629023

</td>

<td style="text-align:right;">

342

</td>

<td style="text-align:right;">

2589628681

</td>

<td style="text-align:right;">

18856026229012708

</td>

<td style="text-align:right;">

137317247

</td>

<td style="text-align:right;">

1443686.2

</td>

<td style="text-align:right;">

293075.5

</td>

<td style="text-align:right;">

8532598

</td>

<td style="text-align:right;">

134.99

</td>

<td style="text-align:right;">

10.64

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

liver

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

36352419

</td>

<td style="text-align:right;">

1722578.5

</td>

<td style="text-align:right;">

1527061499

</td>

<td style="text-align:right;">

1605

</td>

<td style="text-align:right;">

1527059894

</td>

<td style="text-align:right;">

15411382322404590

</td>

<td style="text-align:right;">

124142589

</td>

<td style="text-align:right;">

1305102.1

</td>

<td style="text-align:right;">

266326.8

</td>

<td style="text-align:right;">

18371312

</td>

<td style="text-align:right;">

40.50

</td>

<td style="text-align:right;">

5.93

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

spleen

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

68780464

</td>

<td style="text-align:right;">

1843465.5

</td>

<td style="text-align:right;">

8809244240

</td>

<td style="text-align:right;">

2741

</td>

<td style="text-align:right;">

8809241499

</td>

<td style="text-align:right;">

270266119836188416

</td>

<td style="text-align:right;">

519871253

</td>

<td style="text-align:right;">

5465369.2

</td>

<td style="text-align:right;">

382272.8

</td>

<td style="text-align:right;">

10618617

</td>

<td style="text-align:right;">

121.38

</td>

<td style="text-align:right;">

10.84

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

small-intestine

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

65961991

</td>

<td style="text-align:right;">

1928974.0

</td>

<td style="text-align:right;">

8411139147

</td>

<td style="text-align:right;">

3374

</td>

<td style="text-align:right;">

8411135773

</td>

<td style="text-align:right;">

191380120543171776

</td>

<td style="text-align:right;">

437470137

</td>

<td style="text-align:right;">

4599092.2

</td>

<td style="text-align:right;">

353652.0

</td>

<td style="text-align:right;">

9258160

</td>

<td style="text-align:right;">

132.19

</td>

<td style="text-align:right;">

10.96

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

heart

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

31872599

</td>

<td style="text-align:right;">

2045525.0

</td>

<td style="text-align:right;">

2244720624

</td>

<td style="text-align:right;">

313

</td>

<td style="text-align:right;">

2244720311

</td>

<td style="text-align:right;">

18879520184957976

</td>

<td style="text-align:right;">

137402766

</td>

<td style="text-align:right;">

1444585.3

</td>

<td style="text-align:right;">

415950.0

</td>

<td style="text-align:right;">

15376785

</td>

<td style="text-align:right;">

100.80

</td>

<td style="text-align:right;">

9.39

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

colon

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

29644259

</td>

<td style="text-align:right;">

2078473.0

</td>

<td style="text-align:right;">

2765325256

</td>

<td style="text-align:right;">

4678

</td>

<td style="text-align:right;">

2765320578

</td>

<td style="text-align:right;">

23809529141715324

</td>

<td style="text-align:right;">

154303367

</td>

<td style="text-align:right;">

1622180.2

</td>

<td style="text-align:right;">

373254.8

</td>

<td style="text-align:right;">

10447980

</td>

<td style="text-align:right;">

114.85

</td>

<td style="text-align:right;">

10.18

</td>

<td style="text-align:left;">

0

</td>

</tr>

</tbody>

</table>

![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Distribution-3.png)<!-- -->![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Boxplots%20and%20Density%20Plots%20to%20Measure%20Distribution-4.png)<!-- -->

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

NA\_COUNT

</th>

<th style="text-align:right;">

NA\_FREQ

</th>

<th style="text-align:right;">

MEAN

</th>

<th style="text-align:right;">

MEDIAN

</th>

<th style="text-align:right;">

MAX

</th>

<th style="text-align:right;">

MIN

</th>

<th style="text-align:right;">

MID\_RANGE

</th>

<th style="text-align:right;">

VARIANCE

</th>

<th style="text-align:right;">

STD\_DEV

</th>

<th style="text-align:right;">

SE\_MEAN

</th>

<th style="text-align:right;">

Q1

</th>

<th style="text-align:right;">

Q3

</th>

<th style="text-align:right;">

KURTOSIS

</th>

<th style="text-align:right;">

SKEW

</th>

<th style="text-align:left;">

BC\_LAMBDA

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

tibia

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0005388

</td>

<td style="text-align:right;">

\-0.74

</td>

<td style="text-align:right;">

\-0.73

</td>

<td style="text-align:right;">

1.85

</td>

<td style="text-align:right;">

\-2.16

</td>

<td style="text-align:right;">

4.01

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

0.42

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-1.00

</td>

<td style="text-align:right;">

\-0.48

</td>

<td style="text-align:right;">

2.93

</td>

<td style="text-align:right;">

0.54

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

white-adipose

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

\-0.32

</td>

<td style="text-align:right;">

\-0.46

</td>

<td style="text-align:right;">

14.80

</td>

<td style="text-align:right;">

\-1.67

</td>

<td style="text-align:right;">

16.47

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.81

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.73

</td>

<td style="text-align:right;">

\-0.25

</td>

<td style="text-align:right;">

28.18

</td>

<td style="text-align:right;">

3.81

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

cortex

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

0.0003316

</td>

<td style="text-align:right;">

\-0.38

</td>

<td style="text-align:right;">

\-0.45

</td>

<td style="text-align:right;">

3.21

</td>

<td style="text-align:right;">

\-1.97

</td>

<td style="text-align:right;">

5.18

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:right;">

0.60

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.70

</td>

<td style="text-align:right;">

\-0.25

</td>

<td style="text-align:right;">

5.18

</td>

<td style="text-align:right;">

1.95

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

hypothalmus

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

0.0003079

</td>

<td style="text-align:right;">

\-0.35

</td>

<td style="text-align:right;">

\-0.44

</td>

<td style="text-align:right;">

8.77

</td>

<td style="text-align:right;">

\-1.82

</td>

<td style="text-align:right;">

10.59

</td>

<td style="text-align:right;">

0.49

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.71

</td>

<td style="text-align:right;">

\-0.24

</td>

<td style="text-align:right;">

27.39

</td>

<td style="text-align:right;">

3.81

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

aorta

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

\-0.40

</td>

<td style="text-align:right;">

\-0.44

</td>

<td style="text-align:right;">

3.01

</td>

<td style="text-align:right;">

\-1.39

</td>

<td style="text-align:right;">

4.40

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

0.40

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

\-0.28

</td>

<td style="text-align:right;">

8.10

</td>

<td style="text-align:right;">

2.08

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

plasma

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

0.0005526

</td>

<td style="text-align:right;">

0.08

</td>

<td style="text-align:right;">

\-0.39

</td>

<td style="text-align:right;">

24.34

</td>

<td style="text-align:right;">

\-1.96

</td>

<td style="text-align:right;">

26.30

</td>

<td style="text-align:right;">

1.98

</td>

<td style="text-align:right;">

1.41

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.75

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:right;">

12.03

</td>

<td style="text-align:right;">

2.17

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

testes

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

\-0.18

</td>

<td style="text-align:right;">

\-0.39

</td>

<td style="text-align:right;">

13.90

</td>

<td style="text-align:right;">

\-2.06

</td>

<td style="text-align:right;">

15.96

</td>

<td style="text-align:right;">

0.81

</td>

<td style="text-align:right;">

0.90

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.62

</td>

<td style="text-align:right;">

\-0.10

</td>

<td style="text-align:right;">

24.52

</td>

<td style="text-align:right;">

3.48

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

hippocampus

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

\-0.07

</td>

<td style="text-align:right;">

\-0.36

</td>

<td style="text-align:right;">

6.34

</td>

<td style="text-align:right;">

\-1.64

</td>

<td style="text-align:right;">

7.98

</td>

<td style="text-align:right;">

0.92

</td>

<td style="text-align:right;">

0.96

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.57

</td>

<td style="text-align:right;">

\-0.01

</td>

<td style="text-align:right;">

6.80

</td>

<td style="text-align:right;">

2.46

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

vastus-lateralis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001120

</td>

<td style="text-align:right;">

\-0.12

</td>

<td style="text-align:right;">

\-0.33

</td>

<td style="text-align:right;">

6.09

</td>

<td style="text-align:right;">

\-1.45

</td>

<td style="text-align:right;">

7.54

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

0.87

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.05

</td>

<td style="text-align:right;">

5.47

</td>

<td style="text-align:right;">

2.05

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

lung

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

\-0.08

</td>

<td style="text-align:right;">

\-0.31

</td>

<td style="text-align:right;">

11.63

</td>

<td style="text-align:right;">

\-1.95

</td>

<td style="text-align:right;">

13.58

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

0.82

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.50

</td>

<td style="text-align:right;">

0.03

</td>

<td style="text-align:right;">

12.39

</td>

<td style="text-align:right;">

2.63

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

gastrocnemius

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

\-0.04

</td>

<td style="text-align:right;">

\-0.27

</td>

<td style="text-align:right;">

5.33

</td>

<td style="text-align:right;">

\-1.30

</td>

<td style="text-align:right;">

6.63

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.84

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.59

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:right;">

3.39

</td>

<td style="text-align:right;">

1.66

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

adrenal

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.05

</td>

<td style="text-align:right;">

\-0.24

</td>

<td style="text-align:right;">

7.54

</td>

<td style="text-align:right;">

\-1.24

</td>

<td style="text-align:right;">

8.78

</td>

<td style="text-align:right;">

0.89

</td>

<td style="text-align:right;">

0.94

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.43

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:right;">

9.08

</td>

<td style="text-align:right;">

2.65

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

ovaries

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

0.0044209

</td>

<td style="text-align:right;">

\-0.07

</td>

<td style="text-align:right;">

\-0.24

</td>

<td style="text-align:right;">

6.35

</td>

<td style="text-align:right;">

\-1.31

</td>

<td style="text-align:right;">

7.66

</td>

<td style="text-align:right;">

0.54

</td>

<td style="text-align:right;">

0.73

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.43

</td>

<td style="text-align:right;">

0.07

</td>

<td style="text-align:right;">

16.82

</td>

<td style="text-align:right;">

3.35

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

brown-adipose

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

\-0.21

</td>

<td style="text-align:right;">

11.43

</td>

<td style="text-align:right;">

\-3.37

</td>

<td style="text-align:right;">

14.80

</td>

<td style="text-align:right;">

0.95

</td>

<td style="text-align:right;">

0.97

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.42

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:right;">

17.47

</td>

<td style="text-align:right;">

3.17

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

spleen

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

\-0.18

</td>

<td style="text-align:right;">

9.30

</td>

<td style="text-align:right;">

\-1.42

</td>

<td style="text-align:right;">

10.72

</td>

<td style="text-align:right;">

1.16

</td>

<td style="text-align:right;">

1.08

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.38

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:right;">

11.53

</td>

<td style="text-align:right;">

2.98

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

kidney

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

0.09

</td>

<td style="text-align:right;">

\-0.14

</td>

<td style="text-align:right;">

7.20

</td>

<td style="text-align:right;">

\-1.75

</td>

<td style="text-align:right;">

8.95

</td>

<td style="text-align:right;">

0.64

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.36

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:right;">

4.04

</td>

<td style="text-align:right;">

1.66

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

colon

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:right;">

\-0.13

</td>

<td style="text-align:right;">

11.39

</td>

<td style="text-align:right;">

\-1.59

</td>

<td style="text-align:right;">

12.98

</td>

<td style="text-align:right;">

0.77

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.34

</td>

<td style="text-align:right;">

0.44

</td>

<td style="text-align:right;">

10.99

</td>

<td style="text-align:right;">

2.46

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

liver

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.52

</td>

<td style="text-align:right;">

\-0.10

</td>

<td style="text-align:right;">

10.12

</td>

<td style="text-align:right;">

\-1.43

</td>

<td style="text-align:right;">

11.55

</td>

<td style="text-align:right;">

2.12

</td>

<td style="text-align:right;">

1.46

</td>

<td style="text-align:right;">

0.02

</td>

<td style="text-align:right;">

\-0.44

</td>

<td style="text-align:right;">

1.17

</td>

<td style="text-align:right;">

2.70

</td>

<td style="text-align:right;">

1.56

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

small-intestine

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:right;">

\-0.07

</td>

<td style="text-align:right;">

11.15

</td>

<td style="text-align:right;">

\-1.43

</td>

<td style="text-align:right;">

12.58

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.87

</td>

<td style="text-align:right;">

0.01

</td>

<td style="text-align:right;">

\-0.32

</td>

<td style="text-align:right;">

0.66

</td>

<td style="text-align:right;">

9.65

</td>

<td style="text-align:right;">

1.93

</td>

<td style="text-align:left;">

None

</td>

</tr>

<tr>

<td style="text-align:left;">

heart

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.0001105

</td>

<td style="text-align:right;">

0.51

</td>

<td style="text-align:right;">

\-0.05

</td>

<td style="text-align:right;">

12.30

</td>

<td style="text-align:right;">

\-1.56

</td>

<td style="text-align:right;">

13.86

</td>

<td style="text-align:right;">

2.04

</td>

<td style="text-align:right;">

1.43

</td>

<td style="text-align:right;">

0.02

</td>

<td style="text-align:right;">

\-0.33

</td>

<td style="text-align:right;">

1.03

</td>

<td style="text-align:right;">

11.39

</td>

<td style="text-align:right;">

2.68

</td>

<td style="text-align:left;">

None

</td>

</tr>

</tbody>

</table>

# Impute Missing Values (Necessary for PCA)

# Combine AutoScaled Datasets for PCA (Matrix; Columns:Metabolites, Rows: Unique\_ID) (labelid+INSTITUTE+TECH))

\-PREFA and PREFB are pooled samples for QC purposes \# PCA Analysis
(Shared Samples and Shared Metabolites)

``` r
for(TFORM in c('Untransformed','AutoScaled')){
  if(TFORM == 'Untransformed'){
    tissue_mat <- countdata_df %>%
    ungroup() %>%
    filter(METAB_FAMILY == metab_family) %>%
    filter(NAMED == named) %>%
    select(TISSUE, COUNT_DATA) %>%
    unnest(COUNT_DATA) %>%
    filter(METABOLITE_NAME %in% all_tis_mets) %>%
    mutate(viallabel = ifelse(grepl('QC_',viallabel), paste0(viallabel,'__',TISSUE), viallabel)) %>%
    select(viallabel, METABOLITE_NAME,VALUE) %>%
    pivot_wider(names_from = METABOLITE_NAME, values_from = VALUE) %>%
    tibble::column_to_rownames(var = "viallabel") %>%
    as.matrix()
  }else if(TFORM == 'AutoScaled'){
    tissue_mat <- countdata_df %>%
    ungroup() %>%
    filter(METAB_FAMILY == metab_family) %>%
    filter(NAMED == named) %>%
    select(TISSUE, COUNT_DATA) %>%
    unnest(COUNT_DATA) %>%
    filter(METABOLITE_NAME %in% all_tis_mets) %>%
    mutate(viallabel = ifelse(grepl('QC_',viallabel), paste0(viallabel,'__',TISSUE), viallabel)) %>%
    select(viallabel, METABOLITE_NAME,VALUE) %>%
    pivot_wider(names_from = METABOLITE_NAME, values_from = VALUE) %>%
    tibble::column_to_rownames(var = "viallabel") %>%
    as.matrix() %>% AutoScaleMatrix()
  }
  
  # Determine the proportion of missing data
  na_df <- data.frame()
  for(i in 1:ncol(tissue_mat)){
    # Collect numeric vector
    num_vec <- tissue_mat[,i] %>% unname()
    df <- NumericSummaryStats(num_vec)
    row.names(df) <- colnames(tissue_mat)[i]
    na_df <- rbind(na_df, df)
  }
  na_df %>% select(NA_COUNT, NA_FREQ) %>%
    filter(NA_FREQ > 0) %>%
    arrange(desc(NA_FREQ)) %>% kbl() %>% kable_styling()

  # Collect the major columns to be adjsuted
  imputed_cols10 <- na_df %>% select(NA_COUNT, NA_FREQ) %>%
    filter(NA_FREQ > 0) %>%
    arrange(desc(NA_FREQ)) %>%
    head(n = 10) %>% row.names()

  # Use a knn method to impute missing values (need to come back and adjust)
  preProcValues <- preProcess(tissue_mat,
                              method = c("knnImpute"),
                            k = floor(sqrt(ncol(tissue_mat))),
                            knnSummary = mean)
  # To get the normalizaed data
  impute_tissue_mat <- predict(preProcValues, tissue_mat, na.action = na.pass)
  # To denormalize ("deautoscale")
  tissue_mat_knn <- ReverseAutoScaleMatrix(impute_tissue_mat)

  # Examine the post impute data
  #summary(tissue_mat)
  #summary(tissue_mat_knn)

  # Plot a PCA with outliers labeled
  #dim(tissue_mat_knn)
  #tissue_mat_sam <- tissue_mat_knn[sample(1:nrow(tissue_mat_knn), 200),]

  pca <- prcomp(tissue_mat_knn, scale. = F)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  PC <- pca$x %>% as.data.frame()
  PC$viallabel_tissue <- as.character(row.names(pca$x))
  PC <- PC %>%
    tidyr::separate(col = viallabel_tissue, into = c('viallabel', 'QCtissue'),
    sep = "__", remove = T)
  # Join the phenotype data and additional metadata
  PC <- left_join(PC, pheno_df, by = 'viallabel')
  meta_labelid_df1 <- meta_labelid_df %>%
    filter(sample_type == 'Sample') %>%
    mutate(sample_type = as.character(sample_type))
  meta_labelid_df2 <- meta_labelid_df %>%
    filter(sample_type != 'Sample')
  PC <- left_join(PC, meta_labelid_df1, by = 'viallabel') %>%
    mutate(sample_type = ifelse(grepl('QC_PREFA',viallabel), 'QC-DriftCorrection', sample_type)) %>%
    mutate(sample_type = ifelse(grepl('QC_PREFB',viallabel), 'QC-Pooled', sample_type)) %>%
    mutate(sample_type = ifelse(grepl('QC_Sed',viallabel), 'QC-Reference', sample_type)) %>%
    mutate(TISSUE = ifelse(grepl('QC_PREFA',viallabel), QCtissue, TISSUE)) %>%
    mutate(TISSUE = ifelse(grepl('QC_PREFB',viallabel), QCtissue, TISSUE)) %>%
    mutate(TISSUE = ifelse(grepl('QC_Sed',viallabel), QCtissue, TISSUE))
  table(PC$sample_type) %>% kbl() %>% kable_styling()

  PC %>%
    filter(sample_type %in% c('Sample','QC-Pooled')) %>%
    group_by(sample_type,TISSUE) %>%
      mutate(SAMPLE_N = n()) %>%
      select(sample_type,TISSUE,SAMPLE_N) %>%
      unique() %>% pivot_wider(names_from = sample_type, values_from = SAMPLE_N)

  # Clean up annotations
  PC$sample_type <- as.character(PC$sample_type)
  PC$TISSUE <- as.character(PC$TISSUE)
  PC_DRIFT <- PC %>%
    filter(sample_type %in% c('QC-DriftCorrection','Sample')) %>%
    mutate(QC_DRIFT = ifelse(sample_type == 'QC-DriftCorrection', 'QC-DriftCorrection' ,sample_type))
  PC_Pooled <- PC %>%
    filter(sample_type %in% c('QC-Pooled','Sample')) %>%
    mutate(QC_DRIFT = ifelse(sample_type == 'QC-Pooled', 'QC-Pooled' ,sample_type))
  PC_REF <- PC %>%
    filter(sample_type %in% c('QC-Reference','Sample')) %>%
    mutate(QC_DRIFT = ifelse(sample_type == 'QC-Reference', 'QC-Reference' ,sample_type))

  # Plot the pca (one outlier removed)
  p<- PC_Pooled %>%
    filter(viallabel != '90038016903') %>%
    mutate(dot_size = ifelse(grepl('QC',sample_type), 1.5, 1)) %>%
        ggplot(aes(x = PC1, y = PC2, color=TISSUE, shape = sample_type, size = dot_size)) +
        geom_point(alpha = 0.5) +
        ggtitle(paste0(metab_family,',',named,',',dataset,',','QC_Pooled',',',TFORM)) + 
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance'))
        #theme(aspect.ratio=1)
  print(p)
  p<- PC_REF %>%
    filter(viallabel != '90038016903') %>%
    mutate(dot_size = ifelse(grepl('QC',sample_type), 1.5, 1)) %>%
        ggplot(aes(x = PC1, y = PC2, color=TISSUE, shape = sample_type, size = dot_size)) +
        geom_point(alpha = 0.5) +
        ggtitle(paste0(metab_family,',',named,',',dataset,',','QC_REF',',',TFORM)) + 
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance'))
        #theme(aspect.ratio=1)
  print(p)
 
}
```

![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Impute%20Missing%20Values%20+%20PCA-1.png)<!-- -->![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Impute%20Missing%20Values%20+%20PCA-2.png)<!-- -->![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Impute%20Missing%20Values%20+%20PCA-3.png)<!-- -->![](20201029_pass1a-metabolomics-targeted-cross-tissue_steep_files/figure-gfm/Impute%20Missing%20Values%20+%20PCA-4.png)<!-- -->
