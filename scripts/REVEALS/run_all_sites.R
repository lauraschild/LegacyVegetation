#run all sites with REVEALS
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

REVEALS_df <- lapply(IDs,
                     reveals_analysis_per_site_complete,
                     params = parameters) %>% 
  bind_rows()

REVEALS_df[is.na(REVEALS_df)] <- 0
head(REVEALS_df)
