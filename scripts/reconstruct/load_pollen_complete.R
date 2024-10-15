#script to load continental pollen data

pollen_file <- paste0("reveals_and_psas/input/LegacyPollen2_counts_",tolower(continent),".csv")
unnecessary_cols <- c("Site_Name","Site_Name_Clean","Sample_ID")

pollen_df <- 
  data.table::fread(pollen_file)%>%
  as.data.frame()%>%
  rename(Age_Model_Source = Age_Source)%>%
  dplyr::select(-any_of(unnecessary_cols))

#rescale to percent to sort by abundance
rescale <- t(apply(pollen_df[,-(1:13)],
                   1,
                   function(x) 100*x/sum(x, na.rm = TRUE)))
#remove any columns with no observations
counts <- colSums(rescale)
keep_taxa <- names(sort(counts[counts>0], decreasing = TRUE))
#and put Poaceae at the end to keep from optimizing
keep_taxa <- c(keep_taxa[!(keep_taxa %in% c("Poaceae","Indeterminable"))],"Poaceae","Indeterminable")

#and keep the metadata cols
columns <- c(names(pollen_df[,1:13]),keep_taxa)
pollen_df <- pollen_df[,columns] %>% 
  filter(!Dataset_ID %in% c(17324,
                            17326))

names(pollen_df)[1:13] <- c("Event","Pollen_Data_Source","Site_ID","Dataset_ID","Longitude","Latitude",
                            "Continent","Age_Model_Source","Age_min [yrs BP]", "Age_max [yrs BP]", "Age_median [yrs BP]",
                            "Age_mean [yrs BP]", "Depth [cm]")
rm(rescale)
rm(keep_taxa)