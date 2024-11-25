## UTILITIES

require(tidyverse)
require(magrittr)
require(terra)
require(sf)



##Function to omit uneccesary covariates in the following function


fn_omit <- function(id, mcol, data) {
  f1 <- paste("(count/total) ~", as.character(mcol$Var1[id]))
  f2 <- paste("(count/total) ~", as.character(mcol$Var2[id]))
  
  f1_mod <- glm(
    formula = f1,
    weights = total,
    family = binomial,
    data = data
  )
  
  f2_mod <- glm(
    formula = f2,
    weights = total,
    family = binomial,
    data = data
  )
  
  f1_BIC <- BIC(f1_mod)
  f2_BIC <- BIC(f2_mod)
  
  ifelse(f1_BIC < f2_BIC,
         return(as.character(mcol$Var2[id])),
         return(as.character(mcol$Var1[id])))
  
}

covariate_selection <- function(input) {
  cormat <- cor(input[, 14:ncol(input)])
  cormat[upper.tri(cormat, diag = T)] <- 0
  
  # Check for multicollinearity from the Pearson's correlation
  # coefficient matrix by observing pairs that return > 0.8
  
  mcol <- reshape2::melt(cormat) %>%
    filter(abs(value) > 0.8)
  
  # Compare the BIC of the GLM fitted with the individual covariates
  # from the problematic pairs. Omit the covariate from the model
  # that gives a larger BIC value
  
  omit <- unique(sapply(1:nrow(mcol), function(x)
    fn_omit(x, mcol, input)))
  
  fdat <- input %>% select(4,5,14:ncol(input)) %>% dplyr::select(-all_of(omit))
  
  # Use the step function in the backward direction to select the
  # optimal subset of covariates. Object f_sel is the formula we
  # will use for the INLA part
  
  mod <- glm(
    formula = count / total ~ .,
    weights = total,
    family = "binomial",
    data =fdat
  ) %>%
    step(.,
         direct = "backward",
         k = log(nrow(fdat)),
         trace = 0)
  
  return(mod)
  
}







variable_scale <- function(untransformed_data,round) {
  # Input:
  # Data frame consisting with covariate names appended
  # with the corresponding power transformation number
  
  scaled_cluster_data <- untransformed_data
  scaled_raster_data<- list()
  
  
  covariates <- untransformed_data %>% dplyr::select(-c(1:13))
  scaled_covs <- covariates
  
  for (ii in 1:ncol(covariates)) {
    print(ii)
    covariate <- covariates %>% pull(names(covariates)[[ii]])
    covariate_name <- names(covariates)[[ii]]
    
    grid <- list.files(path = paste0("raster/",round,"/"),
                       pattern = covariate_name,
                       full.names = T) %>%
      rast() %>%
      values()
    
    
    
    scaled_grid <-
      scale(grid)
    
    scaling_attributes <- attributes(scaled_grid)
    
    scaled_covariate <-
      (covariate - scaling_attributes$`scaled:center`) / scaling_attributes$`scaled:scale`
    
    
    
    scaled_cluster_data[[covariate_name]] <-
      scaled_covariate
    scaled_raster_data[[covariate_name]] <-
      as.vector(scaled_grid)
    
    
  }
  
  # Output:
  # A list of scaled raster and cluster level covariates
  
  scaled_cluster_data
  
  return(list(
    scaled_cluster_data,
    as_tibble(scaled_raster_data)
  ))
  
}



