# Range scaling function
################################################
# Range scale a matrix (genes as columns, samples as rows)
RangeScaleMatrix <- function(mat){
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
                mat_out[,n] <- 
                        (mat[,n] - mean(mat[,n],na.rm=T))/(max(mat[,n],na.rm=T) - min(mat[,n],na.rm=T))
        }
        # Output data
        mat_out
}
