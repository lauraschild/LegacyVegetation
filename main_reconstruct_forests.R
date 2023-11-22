#script to reconstruct forest for all three data types
#(Pollen, original REVEALS, optimized REVEALS)

rm(list = ls())
library(dplyr)
source("scripts/supportive/reconstruct_correct_validate_functions.R", echo=TRUE)
continents <- c("Europe")
class <- read.csv2("input/classification.csv")

for(continent in continents){
  #load pollen complete
  source("scripts/reconstruct/load_pollen_complete.R")
  pollen_df <- cbind(pollen_df[,1:13],
                     t(apply(pollen_df[,-(1:13)],
                             1,
                             function(x) 100*x/sum(x, na.rm = TRUE))))
  
  #load REVEALS
  source("scripts/reconstruct/load_REVEALS_complete.R")
  
  #load opti
  source("scripts/reconstruct/load_opti_complete.R")
  
  composition_list <- list(pollen_df,
                           REVEALS_df,
                           opti_df)
  
  #lapply on all dfs to reconstruct tree cover
  forests <- lapply(composition_list,
                    reconstruct_forest_past)
  
  #save as continental
  for(i in 1:length(forests)){
    df <- forests[[i]]
    types <- c("Pollen","REVEALS","opti")
    if(i == 3){
      file <- paste0("output/PANGAEA/composition_forest_",continent,"_",types[i],".csv")
      data.table::fwrite(df,
                         file,
                         row.names = FALSE)
    }
    if(i == 2){
      data.table::fread("output/PANGAEA/Europe_original_RPP.csv") %>% 
        as.data.frame() %>% 
        mutate(`forest [cover in %]` = df$forest) %>% 
        data.table::fwrite("output/PANGAEA/PANGAEA_original_REVEALS_and_forest_Europe.csv")
    }
    if(i==3){
      data.table::fread("output/PANGAEA/optimized_REVEALS_Europe.csv") %>% 
        as.data.frame() %>% 
        mutate(`forest [cover in %]` = df$forest) %>% 
        data.table::fwrite("output/PANGAEA/PANGAEA_optimized_REVEALS_and_forest_Europe.csv")
    }
  }
}


