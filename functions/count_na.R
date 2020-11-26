# Counts the number of NAs in a vector/column
count_na <- function(y) sum(length(which(is.na(y))))
