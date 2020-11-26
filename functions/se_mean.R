# Calculate the standard error of the mean
se_mean <- function(x, na.rm = FALSE) {
        if(na.rm==FALSE){
               sd(x)/sqrt(length(x)) 
        }else{
                sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)])) 
        }}

