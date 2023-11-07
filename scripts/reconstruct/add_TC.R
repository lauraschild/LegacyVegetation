#add reconstructed tree cover to existing REVEALS files and new opti files
rm(list = ls())
library(tidyverse)

meta <- read.csv("data/pollen_meta/metadata.txt", sep = "\t")
new_events <- unique(meta$Event[meta$Project !='LegacyPollen_1.0'])

continents <- c("Africa","Asia","Indopacific","Europe","South_America","North_America")

for(continent in continents){
  #load REVEALS
  dir <- ifelse(continent %in% c("Europe", "North_America","Asia"),
                "NH_sites_orig",
                "SH_sites_orig")
  path <- paste0("reveals_and_psas/output/PANGAEA/",continent,"_original_RPP.csv")
  
  REVEALS_df <- data.table::fread(path)%>%
    as.data.frame()
  
  tree_cover <- data.table::fread(paste0("output/reconstructions/forest_",continent,"_REVEALS.csv"))%>%
    as.data.frame()
  
  merge_cols <- c("Dataset_ID", "Depth [cm]")
  if(continent != "Asia"){
    TC_added <- merge(REVEALS_df,
                      tree_cover[,c(merge_cols,"forest")],
                      by = merge_cols) %>%
      mutate(Event = ifelse(Event %in% new_events,
                            paste0(Event,"_Pollen"),
                            Event)) %>% 
      rename('forest [cover in %]' = forest)
  }else{
    TC_added <- cbind(REVEALS_df,forest = tree_cover$forest)%>%
      mutate(Event = ifelse(Event %in% new_events,
                            paste0(Event,"_Pollen"),
                            Event)) %>% 
      rename('forest [cover in %]' = forest)
  }

  
  
  data.table::fwrite(TC_added,
                     paste0("output/PANGAEA/REVEALS_and_forest_",continent,".csv"),
                     row.names = FALSE)
  # tree_cover %>%
  #   select(forest, Dataset_ID, 'Depth [cm]')%>%
  #   rename('forest [cover in %]' = forest)%>%
  #   data.table::fwrite(paste0("output/PANGAEA/just_forest_",continent,".csv"),
  #                      row.names = FALSE)

  #load optimized REVEALS
  opti <- data.table::fread(paste0("output/reconstructions/optimized_REVEALS_",continent,".csv"))%>%
    as.data.frame()

  rm(tree_cover)

  tree_cover <- data.table::fread(paste0("output/reconstructions/forest_",continent,"_opti.csv"))%>%
    as.data.frame()

  rm(TC_added)

  if(continent != "Asia"){
    TC_added <- merge(opti,
                      tree_cover[,c(merge_cols,"forest")],
                      all.x = TRUE)%>%
      mutate(Event = ifelse(Event %in% new_events,
                            paste0(Event,"_Pollen"),
                            Event)) %>% 
      rename('forest [cover in %]' = forest)
  }else{
    TC_added <- cbind(opti, forest = tree_cover$forest)%>%
      mutate(Event = ifelse(Event %in% new_events,
                            paste0(Event,"_Pollen"),
                            Event)) %>% 
      rename('forest [cover in %]' = forest)
  }


  data.table::fwrite(TC_added,
                     paste0("output/PANGAEA/opti_and_forest_",continent,".csv"),
                     row.names = FALSE)

}
