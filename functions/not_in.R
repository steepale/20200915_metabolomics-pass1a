# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
  !('%in%'(x,y))
}
################################################################################

