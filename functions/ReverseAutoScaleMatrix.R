# Reverse Auto scaling function
################################################
# Reverse Auto scale a matrix (genes as columns, samples as rows)
ReverseAutoScaleMatrix <- function(mat_org, mat_norm){
        # Ensure data is matrix
        #if(!is.matrix(mat)){
        #        stop("Data must be matrix")
        #}
        # Ensure data is numeric
        #if(!is.numeric(mat)){
        #        stop("Data must be numeric")
        #}
        # Iterate over columns
        mat_out <- mat_norm
        for(n in 1:ncol(mat_norm)){
                # Transform
                mat_out[,n] <- 
                        (mat_norm[,n] * sd(mat_org[,n],na.rm=T))  + mean(mat_org[,n],na.rm=T)
        }
        # Output data
        mat_out
}
