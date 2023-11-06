#check some things out

rm(list = ls())
library(tidyverse)

#how much meta data is missing
continents <- c("Asia","North_America","Europe","South_America","Indopacific","Africa")

taxa_list <- c()
ID_list <- c()

for(continent in continents){
  df <- read.csv(paste0("data/pollen_data/LegacyPollen2_counts_",tolower(continent),".csv"),
                   sep = "\t")
  
  taxa <- names(df)[-(1:16)]
  taxa_list <- c(taxa_list, taxa)
  
  IDs <- unique(df$Dataset_ID)
  ID_list <- c(ID_list, IDs)
}

#load classification file
class <- read.csv2("data/classification2.csv")
missing <- unique(taxa_list)[!(unique(taxa_list) %in% class$taxa)]

source("scripts/supportive/edit_meta.R", echo=FALSE)
sum(!(ID_list %in% meta_df$Dataset_ID))
