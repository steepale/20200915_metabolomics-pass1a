# Collects the range of a vector and outputs the range as a pretty string
range_chr <- function(x){
  rng <- range(x, na.rm = T)
  paste0(rng[1],'-',rng[2])
}
