#run all sites with REVEALS
library(parallel)
#run REVEALS on all Dataset_IDs
#source("scripts/supportive/run_REVEALS.R", echo=FALSE)
#actual REVEALS analysis for all sites
IDs <- unique(pollen_df$Dataset_ID)
IDs <- IDs[IDs %in% unique(meta_df$Dataset_ID)]

parameters <- parameters %>%
  select(-Continent) %>%
  rename(species = Taxon,
         fallspeed = Fallspeed, 
         PPE = PPE, 
         PPE.error = PPE_Error)

source("scripts/REVEALS/reveals_per_site_complete.R", echo=TRUE)
source("scripts/REVEALS/reveals_and_PSA_functions.R", echo = FALSE)

cl <- makeCluster(n)
clusterExport(cl,
              list("parameters",
                   "pollen_df",
                   "optimize",
                   "meta_df",
                   "taxa",
                   "class",
                   "parameters",
                   "reveals_analysis_per_site_complete",
                   "rmultinom_reveals",
                   "rnorm_reveals",
                   "K_factors"),
              envir = environment())

REVEALS_df <- clusterApplyLB(cl,
                             IDs,
                             reveals_analysis_per_site_complete,
                             params = parameters)%>%
  bind_rows()

stopCluster(cl)
REVEALS_df[is.na(REVEALS_df)] <- 0
head(REVEALS_df)
