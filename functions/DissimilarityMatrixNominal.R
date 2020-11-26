# Function to calculate the dissimilarity matrix between nominal variables
DissimilarityMatrixNominal <- function(df){
  # Variables
  f_set <- names(df)
  d <- list()
  # Iterate through each feature
  for(f in f_set){
    d[[f]] <- matrix(nrow = length(df[[f]]), ncol = length(df[[f]]))
    for(i in 1:length(df[[f]])){
      for(j in 1:length(df[[f]])){
        p = 1
        if(df[[f]][i] == df[[f]][j]){
          m = 1
        }else{
          m = 0
        }
        d[[f]][i,j] <- (p - m)/p
      }
    }
  }
  # After iteration
  Reduce('+', d)/length(f_set)
}
