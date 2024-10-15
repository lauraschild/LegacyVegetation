#script to reconstruct forest for all three data types
#(Pollen, original REVEALS, optimized REVEALS)

rm(list = ls())
library(dplyr)
source("scripts/supportive/reconstruct_correct_validate_functions.R", echo=TRUE)
continents <- c("Europe")
class <- read.csv2("input/classification_revised.csv")

path_check <- dir.exists("output/reconstructions")

if(!path_check) dir.create("output/reconstructions")

#load meta data with basin type
meta <- data.table::fread("reveals_and_psas/input/metadata2.csv") %>% 
  as.data.frame() %>% 
  dplyr::select(Dataset_ID, Continent, Basin_Type,Latitude,Longitude, Basin_Area) %>% 
  filter(Continent %in% c("Eastern North America","Western North America","Asia","Europe")) %>% 
  mutate(record_type = ifelse(Basin_Type == "Lake" & Basin_Area >= 500000,
                              "large lake",
                              ifelse(Basin_Type == "Lake",
                                     "small lake",
                                     ifelse(Basin_Type == "Peat",
                                            "peatland",
                                            "no use"))),
         record_type = ifelse(is.na(record_type),
                              "small lake",
                              record_type),
         Continent = ifelse(Continent %in% c("Europe","Asia"),
                            Continent,
                            "North America")) %>% 
  filter(record_type != "no use", !is.na(record_type)) 


#### POLLEN ####
get_pollen <- function(continent){
  source("scripts/reconstruct/load_pollen_complete.R",
         local = TRUE)
  
  pollen_df[,-(1:13)] <- t(apply(pollen_df[,-(1:13)],
                                 1,
                                 function(x) 100*x/sum(x,na.rm = TRUE)))
  
  continental_df <- pollen_df %>% 
    rename(Age_BP = `Age_mean [yrs BP]`) %>% 
    dplyr::select(-Event,-Pollen_Data_Source, -Site_ID,
                  -Age_Model_Source,-ends_with("]")) %>% 
    merge(meta[,c("Dataset_ID","record_type")],
          by = "Dataset_ID",
          all.x = TRUE) %>% 
    reconstruct_forest()
  
  #save forest cover and composition in cont df in reconstruction directory
  data.table::fwrite(continental_df,
                     paste0("output/reconstructions/pollen_forest_",
                            continent,".csv"))
  
  continental_df %>% 
    dplyr::select(1:5,record_type, forest) %>% 
    return()
}

pollen_forest_site <-lapply(continents,
                            get_pollen) %>% 
  bind_rows() 


#### REVEALS####
#retrival function
get_reveals <- function(continent){
  source("scripts/reconstruct/load_REVEALS_with_SD.R",
         local = TRUE)
  
  #make df with composition and forest cover
  continental_df <- REVEALS_df %>% 
    rename(Age_BP = `Age_mean [yrs BP]`) %>% 
    dplyr::select(-Event,-Pollen_Data_Source, -Site_ID,
                  -Age_Model_Source,-ends_with("]")) %>% 
    merge(meta[,c("Dataset_ID","record_type")]) %>% 
    reconstruct_forest_with_SD()
  
  #save in reconstruction directory
  data.table::fwrite(continental_df,
                     paste0("output/reconstructions/REVEALS_with_SD_forest_",
                            continent,".csv"))
  
  #only return with forest cover
  continental_df %>% 
    dplyr::select(1:5,record_type, forest_mean) %>% 
    return()
}

reveals_forest_site <- lapply(continents,
                              get_reveals) %>% 
  bind_rows()

#save just reconstructed forest cover
reveals_forest_site %>%
  data.table::fwrite("output/reconstructions/REVEALS_with_SD_forest_complete.csv")

#filter pollen df to only include same records as reveals
pollen_forest_site <- pollen_forest_site %>% 
  filter(Dataset_ID %in% unique(reveals_forest_site$Dataset_ID))

#save just reconstructed forest cover
pollen_forest_site %>% 
  data.table::fwrite("output/reconstructions/pollen_forest_complete.csv")
