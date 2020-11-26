# Grab a linear equation from simple linear regression and robust
######################################################
lm_eqn <- function(met_df){
library(MASS)
  library(sfsmisc)
  met <- unique(met_df$name)
  mod <- NULL
  try(mod <- lm(y ~ x, met_df), silent = T)
  b <- NULL
  try(b <- format(unname(coef(mod)[1]), digits = 2), silent = T)
  m <- NULL
  try(m <- format(unname(coef(mod)[2]), digits = 2), silent = T)
  r2 <- NULL
  try(r2 <- format(summary(mod)$r.squared, digits = 3), silent = T)
  eqn <- NULL
  try(eqn <- paste0('y = ',m,'x + ',b), silent = T)
  rlm_mod <- NULL
  try(summary(rlm_mod <- rlm(y ~ x, met_df)), silent = T)
  try(pval <- f.robftest(rlm_mod, var = "x")$p.value, silent = T)
  try(fstat <- f.robftest(rlm_mod, var = "x")$statistic %>% unname(), silent = T)
  if(is.null(rlm_mod) & is.null(mod)){
    out <- data.frame(METABOLITE = met, EQUATION = NA, R2 = NA, F_TEST_RLM = NA, PVAL_RLM = NA)
  }else if(is.null(rlm_mod)){
    out <- data.frame(METABOLITE = met, EQUATION = eqn, R2 = r2, F_TEST_RLM = NA, PVAL_RLM = NA)
  }else{
    out <- data.frame(METABOLITE = met, EQUATION = eqn, R2 = r2, F_TEST_RLM = fstat, PVAL_RLM = pval)
  }
  out
}
