#optimization function

optimize <- function(RPPs){#vector of RPPs to be optimized
  
  #source functions needed
  library(dplyr)
  source("scripts/REVEALS/REVEALSinR_modified.R", echo = FALSE)
  source("scripts/REVEALS/reveals_and_PSA_functions.R", echo = FALSE)
  source("scripts/supportive/reconstruct_correct_validate_functions.R", echo=FALSE)
  
  #print(RPPs)
  #replace first ten RPPs with input RPPs in parameters df
  REVEALS_parameters <- parameters
  REVEALS_parameters$PPE[1:m] <- RPPs
  REVEALS_parameters <- REVEALS_parameters %>%
    select(-Continent) %>%
    rename(species = Taxon,
           fallspeed = Fallspeed, 
           PPE = PPE, 
           PPE.error = PPE_Error)
  
  #run REVEALS on all Dataset_IDs
  #source("scripts/supportive/run_REVEALS.R", echo=FALSE)
  #actual REVEALS analysis for all sites
  IDs <- unique(pollen_df$Dataset_ID)
  IDs <- IDs[IDs %in% unique(meta_df$Dataset_ID)]
  source("scripts/REVEALS/reveals_per_site.R", echo=FALSE)
  # cl <- makeCluster(n)
  # clusterExport(cl,
  #               list("parameters",
  #                    "pollen_df",
  #                    "RS",
  #                    "urban",
  #                    "optimize",
  #                    "meta_df",
  #                    "taxa",
  #                    "class",
  #                    "REVEALS_parameters",
  #                    "reveals_analysis_per_site",
  #                    "rmultinom_reveals",
  #                    "rnorm_reveals",
  #                    "K_factors"),
  #               envir = environment())
  # 
  # REVEALS_df <- clusterApplyLB(cl,
  #                              IDs,
  #                             reveals_analysis_per_site,
  #                             params = REVEALS_parameters,
  #                             parallel = FALSE)%>%
  #   bind_rows()
  
  REVEALS_df <- lapply(IDs,
                       reveals_analysis_per_site,
                       params = REVEALS_parameters) %>% 
    bind_rows()
  
  REVEALS_df[is.na(REVEALS_df)] <- 0
  
  #reconstruct forest
  forest_df <- reconstruct_forest(REVEALS_df)
  #correct for urban areas
  correct_forest_df <- correct_urban(forest_df)
  #validate
  validation_result <- validate_forest(correct_forest_df)
  
  return(validation_result$RSS)
}
