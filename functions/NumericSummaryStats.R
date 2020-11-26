# A summary function that takes a vector of numeric values and provides summary statistics
NumericSummaryStats <- function(num_vec){
  # Load libraries
  ##############################################################################
  library(caret)
  library(tidyverse)
  library(e1071)
  
  # Source sub functions
  ##############################################################################
  # Counts the number of NAs in a vector/column
  count_na <- function(y) sum(length(which(is.na(y))))
  
  # method = "mode" [default]: calculates the mode for a unimodal vector, else returns an NA
  # method = "nmodes": calculates the number of modes in the vector
  # method = "modes": lists all the modes for a unimodal or polymodal vector
  # https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
  modeav <- function (x, method = "mode", na.rm = FALSE)
  {
    x <- unlist(x)
    if (na.rm)
      x <- x[!is.na(x)]
    u <- unique(x)
    n <- length(u)
    #get frequencies of each of the unique values in the vector
    frequencies <- rep(0, n)
    for (i in seq_len(n)) {
      if (is.na(u[i])) {
        frequencies[i] <- sum(is.na(x))
      }
      else {
        frequencies[i] <- sum(x == u[i], na.rm = TRUE)
      }
    }
    #mode if a unimodal vector, else NA
    if (method == "mode" | is.na(method) | method == "")
    {return(ifelse(length(frequencies[frequencies==max(frequencies)])>1,NA,u[which.max(frequencies)]))}
    #number of modes
    if(method == "nmode" | method == "nmodes")
    {return(length(frequencies[frequencies==max(frequencies)]))}
    #list of all modes
    if (method == "modes" | method == "modevalues")
    {return(u[which(frequencies==max(frequencies), arr.ind = FALSE, useNames = FALSE)])}  
    #error trap the method
    warning("Warning: method not recognised.  Valid methods are 'mode' [default], 'nmodes' and 'modes'")
    return()
  }
  # Calculate the standard error of the mean
  se_mean <- function(x, na.rm = FALSE) {
        if(na.rm==FALSE){
               sd(x)/sqrt(length(x)) 
        }else{
                sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)])) 
        }}
  
  # Generate a summary table to summazrize the density plots resulting from each transfrmation   technique
  ################################################################################
  # Summary Datafrmae to spit out
  summary_out <- data.frame()

  # NA Values
  ################################################################
  # The number of NA values
  NA_COUNT <- count_na(num_vec)
  # The percent NA values
  NA_FREQ <- count_na(num_vec)/length(num_vec) %>% round(digits = 2)
                
  # Central Tendency
  ################################################################
  # Calculate the mean for all numeric data
  MEAN <- mean(num_vec, na.rm = TRUE)  %>% round(digits = 2)
  # Calculate the trimmed mean for all numeric data
  TRIMMED_MEAN <- mean(num_vec, na.rm = TRUE, trim = 0.02) %>% round(digits = 2)
  # Calculate the median
  MEDIAN <- median(num_vec, na.rm = TRUE) %>% round(digits = 2)
  # Calculate the mode(s)
  # MODE[[i]] <- df %>%
  #         map(.f = Mode)
  # Calculate the min & max
  MAX <- max(num_vec, na.rm = TRUE) %>% round(digits = 2)
  MIN <- min(num_vec, na.rm = TRUE) %>% round(digits = 2)
  # Calculate the midrange
  MID_RANGE <- MAX - MIN %>% round(digits = 2)

  # Variability (Spread)
  ################################################################
  # Calculate the variance
  VARIANCE <- var(num_vec, na.rm = TRUE) %>% round(digits = 2)
  # Calculate the standard deviation
  STD_DEV <- sqrt(VARIANCE) %>% round(digits = 2)
  # Calculate the standard error of the mean
  SE_MEAN <- se_mean(num_vec, na.rm = TRUE) %>% round(digits = 2)
                
  # Relative Standing (Distribution)
  ################################################################
  # Calculate the quantiles
  Q1 <- quantile(num_vec, na.rm = TRUE, probs = 0.25) %>% round(digits = 2)
  Q3 <- quantile(num_vec, na.rm = TRUE, probs = 0.75) %>% round(digits = 2)
  IQR <- paste0(Q1, ' - ', Q3)
  # Calculate the skew
  if(MAX == MIN){
    SKEW <- 0 %>% round(digits = 2)
  }else{
    SKEW <- e1071::skewness(num_vec, na.rm = TRUE) %>% round(digits = 2)
  }
  KURTOSIS <- e1071::kurtosis(num_vec, na.rm = TRUE) %>% round(digits = 2)
  
  # Transformations (consider for uncentered and scaled counts)
  ##############################################################################
  if(!is.null((preProcess(as.data.frame(num_vec), method = c("BoxCox", "center", "scale"),
                  na.remove = TRUE))$bc[['num_vec']]$lambda) ){
    BC_LAMBDA <- (preProcess(as.data.frame(num_vec), 
                  method = c("BoxCox", "center", "scale"),
                  na.remove = TRUE))$bc[['num_vec']]$lambda %>% as.character()
  }else{
    BC_LAMBDA <- "None"
  }
  # Create an output dataframe of features
  ##############################################################################
  summary_out <- data.frame(NA_COUNT,NA_FREQ, MEAN, MEDIAN,
                     MAX, MIN, MID_RANGE, VARIANCE, STD_DEV, SE_MEAN, Q1, Q3,
                     KURTOSIS, SKEW,BC_LAMBDA)
  summary_out
}
