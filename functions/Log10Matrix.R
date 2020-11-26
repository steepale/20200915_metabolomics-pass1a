# Log 10 transform function
################################################
# Log transform a matrix (genes as columns, samples as rows)
# Note: adds 1 to values prior to log10
Log10Matrix <- function(mat){
        # Ensure data is matrix
        if(!is.matrix(mat)){
                stop("Data must be matrix")
        }
        # Ensure data is numeric
        if(!is.numeric(mat)){
                stop("Data must be numeric")
        }
        # Iterate over columns
        mat_out <- mat
        for(n in 1:ncol(mat)){
                # Apply transformation
                mat_out[,n] <- log10(mat[,n] + 1)
        }
        # Output data
        mat_out
}

