#main
rm(list = ls())
library(dplyr)
#set number of clusters
n <- 20

#numer of taxa to optimize
m <- 10

#loop through all continents
continents <- c("Indopacific","Asia","Africa","South_America","Europe","North_America")

#load meta data
source("scripts/supportive/edit_meta.R", echo=FALSE)

#load remote sensing tree cover
RS <- read.csv("data/landsat_treecover.csv")%>%
  select(Dataset_ID, mean)

#load urban correction
urban <- read.csv("data/openness.csv")

#load classification file
class <- read.csv2("data/classification.csv")

#load look-up table for K factors
K_factors <- read.csv("data/K_factors.csv")

for(continent in continents){
  #load pollen files
  source("scripts/supportive/load_pollen.R", echo=FALSE)
  
  #prep parameter df and bounds
  source("scripts/supportive/prep_pollen_parameters.R", echo=FALSE)
  
  #run optimization
  source("scripts/supportive/run_optimization.R", echo=TRUE)
  
  #save optimized parameters
  write.csv(data.frame(Taxa = parameters$Taxon[1:m],
                       PPE = results$par),
            paste0("output/optimized_RPP_",continent,".csv"),
            row.names= FALSE)
  
  #run REVEALS with final parameters and reconstruct tree cover (SCRIPT)
  #validate and save final plot with MAE and RSS (SCRIPT)
  source("scripts/supportive/reconstruct_correct_validate_functions.R", echo=FALSE)
  
  #run reveals for all sites with new parameters
  source("scripts/supportive/final_valid.R", echo=FALSE)
  
  
}

