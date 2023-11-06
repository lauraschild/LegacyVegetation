#load meta data
#check which Events were present in LegaPollen1.0
#rename other Events in new data set
#append '_Pollen'

rm(list = ls())
library(dplyr)
meta <- read.csv("data/pollen_meta/metadata.txt", sep = "\t")
names(meta)
new_events <- unique(meta$Event[meta$Project !='LegacyPollen_1.0'])
new_events

continents <- c("Asia","Europe","Indopacific","South_America","North_America")

for(continent in continents){
  
  file_reveals <- paste0("output/PANGAEA/REVEALS_and_forest_",continent,".csv")
  data.table::fread(file_reveals)%>%
    as.data.frame()%>%
    mutate(Event = ifelse(Event %in% new_events,
                          paste0(Event,"_Pollen"),
                          Event))%>%
    data.table::fwrite(paste0("output/PANGAEA/REVEALS_and_forest_newEvents_",continent,".csv"),
                       row.names = FALSE)
  
  file_opti <- paste0("output/PANGAEA/opti_and_forest_",continent,".csv")
  data.table::fread(file_opti)%>%
    as.data.frame()%>%
    mutate(Event = ifelse(Event %in% new_events,
                          paste0(Event,"_Pollen"),
                          Event))%>%
    data.table::fwrite(paste0("output/PANGAEA/opti_and_forest_newEvents_",continent,".csv"),
                       row.names = FALSE)
  
}