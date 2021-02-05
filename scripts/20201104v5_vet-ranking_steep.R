# Vetting script
# TISSUE: brown-adipose
# NORM_ORDER: feat_sam
# N_LOG: none
# S_CENTER: none
# S_SCALE: none
# F_CENTER: mean
# F_SCALE: sd

# Load dependencies
pacs...man <- c("tidyverse", "data.table","rlang","RColorBrewer","MASS","sfsmisc",
                "pheatmap","caret","ggforce","RANN","kableExtra","gtools","combinat",
                "knitr","markdown","rmarkdown","gridExtra","cowplot","devtools")
for(pac in pacs...man){
  suppressWarnings(suppressPackageStartupMessages(library(pac, character.only = TRUE)))
}

# Load data
d_dir <- '/Volumes/Frishman_4TB/motrpac/20200915_metabolomics-pass1a/files_to_share'
um1 <- readRDS(paste0(d_dir, '/20201103_named-site-comparison-0001_brown-adipose_m70-s77_university-of-michigan-rppos_steep.RDS'))
br1 <- readRDS(paste0(d_dir, '/20201103_named-site-comparison-0001_brown-adipose_m70-s77_broad-institute-hilicpos_steep.RDS'))

# Ensure that orders of columns match
if(!all(colnames(br1) == colnames(um1))){
  stop('Columns out of order')
}
# Ensure that orders of rownames match
if(!all(rownames(br1) == rownames(um1))){
  stop('Rownames out of order')
}

# Dimensions
dim(um1)
dim(br1)

# Visualize cross correlation matrixes before any adjustments
################################
paletteLength <- 100
redblue100 <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# Visualize the feature feature matrix
image(cor(cbind(um1,br1),method="spearman"),col=redblue100)
image(cor(cbind(t(um1),t(br1)),method="spearman"),col=redblue100)

# Visualize Abundances with boxplots
################################
boxplot(data.frame(log(br1)))
boxplot(data.frame(log(um1)))
boxplot(data.frame(log(t(br1))))
boxplot(data.frame(log(t(um1))))

# Count Zeros across features and samples
##########################
CountZeros <- function(x){
  sum(x==0, na.rm = T)
}
# Features
f0.br1 <- apply(br1, 2, CountZeros) %>% sort() %>% rev()
f0.br1
f0.um1 <- apply(um1, 2, CountZeros) %>% sort() %>% rev()
f0.um1
# Samples
s0.br1 <- apply(br1, 1, CountZeros) %>% sort() %>% rev()
s0.br1
s0.um1 <- apply(um1, 1, CountZeros) %>% sort() %>% rev()
s0.um1

# Count missing values across features and samples
##########################
CountNAs <- function(x){
  sum(is.na(x))
}

# Features
is.na(um1) %>% table()
fna.br1 <- apply(br1, 2, CountNAs) %>% sort() %>% rev()
fna.br1
fna.um1 <- apply(um1, 2, CountNAs) %>% sort() %>% rev()
fna.um1
# Samples
sna.br1 <- apply(br1, 1, CountNAs) %>% sort() %>% rev()
sna.br1
sna.um1 <- apply(um1, 1, CountNAs) %>% sort() %>% rev()
sna.um1

# Remove Features
##############################
# Remove Glutathione, CAR(5:0)
rm.br1 <- match(c('Glutathione', 'CAR(5:0)'), colnames(br1))
rm.um1 <- match(c('Glutathione', 'CAR(5:0)'), colnames(um1))
if(!all(rm.br1 == rm.um1)){
  stop('Indexes out of order')
}
br2<- br1[,-rm.br1]
um2<- um1[,-rm.um1]

# Impute Features
#############################
# No need yet

# Examine the feature/sample mean/median
#########################################
# Features
f.mean.um<-apply(um2,2,mean)
f.mean.br<-apply(br2,2,mean)
f.median.um<-apply(um2,2,median)
f.median.br<-apply(br2,2,median)
# Samples
s.median.um<-apply(um2,1,median)
s.median.br<-apply(br2,1,median)
s.mean.um<-apply(um2,1,mean)
s.mean.br<-apply(br2,1,mean)
# Plots
plot(s.median.br,s.median.um)
plot(f.median.br,f.median.um)
plot(s.mean.br,s.mean.um)
plot(f.mean.br,f.mean.um)

# Perform transformations
#############################
# # Log transform
# um3 <- log(um2)
# br3 <- log(br2)

# linear
um3 <- um2
br3 <- br2

# By Feature:
# Center/Scale
#f.median.um3 <- apply(um3, 2, median)
f.mean.um3 <- apply(um3, 2, mean)
f.sd.um3 <- apply(um3, 2, sd)
f.mean.br3 <- apply(br3, 2, mean)
f.sd.br3 <- apply(br3, 2, sd)
for(i in 1:dim(um2)[2]){
  um3[,i] <- (um3[,i]-f.mean.um3[i])/f.sd.um3[i]
  br3[,i] <- (br3[,i]-f.mean.br3[i])/f.sd.br3[i]
}

# # By Sample:
# # Center/Scale
# s.median.um3 <- apply(um3, 1, median)
# s.sd.um3 <- apply(um3, 1, sd)
# s.median.br3 <- apply(br3, 1, median)
# s.sd.br3 <- apply(br3, 1, sd)
# for(i in 1:dim(um2)[1]){
#   um3[i,] <- (um3[i,]-s.median.um3[i])/s.sd.um3[i]
#   br3[i,] <- (br3[i,]-s.median.br3[i])/s.sd.br3[i]
# }

# Visualize cross correlation matrixes after adjustments
################################
# Visualize the feature feature matrix
image(cor(cbind(um3,br3),method="spearman"),col=redblue100)
image(cor(cbind(t(um3),t(br3)),method="spearman"),col=redblue100)

# Visualize Abundances with boxplots
################################
boxplot(data.frame(br3))
boxplot(data.frame(um3))
boxplot(data.frame(t(br3)))
boxplot(data.frame(t(um3)))

# Examine the feature/sample mean/median
#########################################
# Features
f.mean.um3<-apply(um3,2,mean)
f.mean.br3<-apply(br3,2,mean)
f.median.um3<-apply(um3,2,median)
f.median.br3<-apply(br3,2,median)
# Samples
s.median.um3<-apply(um3,1,median)
s.median.br3<-apply(br3,1,median)
s.mean.um3<-apply(um3,1,mean)
s.mean.br3<-apply(br3,1,mean)
# Plots
plot(s.median.br3,s.median.um3)
plot(f.median.br3,f.median.um3)
plot(s.mean.br3,s.mean.um3)
plot(f.mean.br3,f.mean.um3)

# Perform Rank calculation
###############################
# Features
######################
fxf136 <- cor(cbind(um3,br3), method = 'spearman')
dim(fxf136)
rank68.a<-apply(fxf136[1:68,69:136],1,rank)
rank68.b<-apply(fxf136[1:68,69:136],2,rank)

# verified as linear:
plot(fxf136[1:68,69],rank68.b[,1])
rank68.a2<-rank68.b2<-rep(1,68)
for(i in 1:68){
  rank68.a2[i]<-rank68.a[i,i]
  rank68.b2[i]<-rank68.b[i,i]
}
table(rank68.a2,rank68.b2)
rank68.b2
rank68.a2
rank68.m2<-(rank68.a2+rank68.b2)/136
rank68.m3<-((rank68.a2/68)+(rank68.b2/68))/2
plot(rep(1,68),pch=16,cex=(rank68.m2-0.4)*2.5)

mean(rank68.m2)
# Samples
######################
sxs154 <- cor(cbind(t(um3),t(br3)), method = 'spearman')
dim(sxs154)
rank77.a<-apply(sxs154[1:77,78:154],1,rank)
rank77.b<-apply(sxs154[1:77,78:154],2,rank)

# verified as linear:
plot(sxs154[1:77,78],rank77.b[,1])
rank77.a2<-rank77.b2<-rep(1,77)
for(i in 1:77){
  rank77.a2[i]<-rank77.a[i,i]
  rank77.b2[i]<-rank77.b[i,i]
}
table(rank77.a2,rank77.b2)
rank77.b2
rank77.a2
rank77.s2<-(rank77.a2+rank77.b2)/154
plot(rep(1,length(rank77.s2)),pch=16,cex=(rank77.s2-0.4)*2.5)

mean(rank68.m2)
mean(rank77.s2)

