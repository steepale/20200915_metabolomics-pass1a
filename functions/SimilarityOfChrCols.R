# A function that takes 2 nested dataframes of a dataframe as input and calculates their degree of similarity
SimilarityOfChrCols <- function(x,y){
  total <- (unname(unlist(x)) %in% unname(unlist(y))) %>% table() %>% sum()
  tru_n <- (unname(unlist(x)) %in% unname(unlist(y))) %>% sum()
  perc_sim <- round(tru_n/total, digits = 3)
  perc_sim
}
