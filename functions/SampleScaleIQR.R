################################################
# Scale a matrix (genes as columns, samples as rows)
SampleScaleIQR <- function(mat){
        # Ensure data is matrix
        if(!is.matrix(mat)){
                stop("Data must be matrix")
        }
        # Ensure data is numeric
        if(!is.numeric(mat)){
                stop("Data must be numeric")
        }
        # Iterate over rows
        mat_out <- mat
        for(n in 1:nrow(mat)){
          # Scale
          Q1 <- quantile(mat[n,], na.rm = TRUE, probs = 0.25)
          Q3 <- quantile(mat[n,], na.rm = TRUE, probs = 0.75)
          IQR <- Q3-Q1
          mat_out[n,] <- mat[n,]/IQR
        }
        # Output data
        mat_out
}
