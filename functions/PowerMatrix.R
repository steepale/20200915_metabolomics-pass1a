# Power transform function
################################################
# Power transform a matrix (genes as columns, samples as rows)
PowerMatrix <- function(mat){
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
                mat_out[,n] <- sqrt(mat[,n])
        }
        # Output data
        mat_out
}
