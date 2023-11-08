#script to run REVEALs with "optimized" parameters
#and reconstruct forest cover for original and optimized REVEALS run
rm(list = ls())

continent <- "Europe"

#load meta data
source("scripts/supportive/edit_meta.R")

#load classification file
class <- read.csv2("input/classification.csv")

#load look-up table for K factors
K_factors <- read.csv("input/K_factors.csv")

source("scripts/reconstruct/load_pollen_complete.R")

source("scripts/supportive/get_new_parameters.R")
parameters$Continent <- continent
head(parameters)

write.csv(parameters,
          paste0("output/PANGAEA/RPP_",continent,".csv"),
          row.names = FALSE)
source("scripts/REVEALS/run_all_sites.R")

write.csv(REVEALS_df,
          paste0("output/PANGAEA/optimized_REVEALS_",continent,".csv"),
          row.names = FALSE)

