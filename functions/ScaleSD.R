################################################
# Scale a matrix (genes as columns, samples as rows)
ScaleSD <- function(mat){
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
          # Scale
          mat_out[,n] <- mat[,n]/sd(mat[,n],na.rm=T)
        }
        # Output data
        mat_out
}
