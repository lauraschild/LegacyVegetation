#use final parameters in validation and plot
source("scripts/REVEALS/REVEALSinR_modified.R", echo=FALSE)
source("scripts/REVEALS/reveals_and_PSA_functions.R", echo = FALSE)
source("scripts/supportive/reconstruct_correct_validate_functions.R", echo=FALSE)
source("scripts/REVEALS/reveals_per_site.R", echo=TRUE)


IDs <- unique(pollen_df$Dataset_ID)
IDs <- IDs[IDs %in% unique(meta_df$Dataset_ID)]

parameters$PPE[1:m] <- results$par

parameters <- parameters %>%
  select(-Continent) %>%
  rename(species = Taxon,
         fallspeed = Fallspeed, 
         PPE = PPE, 
         PPE.error = PPE_Error)
#run REVEALS in parallel for last validation
library(parallel)
cl <- makeCluster(n)
clusterExport(cl,
              list("parameters",
                   "pollen_df",
                   "reveals_analysis_per_site",
                   "meta_df",
                   "taxa",
                   "rmultinom_reveals",
                   "rnorm_reveals",
                   "K_factors"))
REVEALS_df <- clusterApplyLB(cl = cl,
                             IDs,
                             reveals_analysis_per_site,
                             params = parameters)%>%
  bind_rows()

stopCluster(cl)

REVEALS_df[is.na(REVEALS_df)] <- 0


#reconstruct forest
forest_df <- reconstruct_forest(REVEALS_df)
#correct for urban areas
correct_forest_df <- correct_urban(forest_df)
#validate
validation_result <- validate_forest(correct_forest_df)
#plot
plot_valid(validation_result, convergence = results$convergence)
