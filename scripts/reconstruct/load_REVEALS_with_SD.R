#load REVEALS complete
#script to read REVEALS data


path <- paste0("output/PANGAEA/",continent,"_revision.csv")

REVEALS_df <- data.table::fread(path)%>%
  as.data.frame()%>%
  dplyr::select(1:13,
         contains("[mean of cover"),
         contains("[sd of cover")) %>% 
  filter(!Dataset_ID %in% c(17324,
                            17326) )

names(REVEALS_df) <- gsub(" [mean of cover in %]",
                          "_mean",
                          names(REVEALS_df),
                          fixed = TRUE)
names(REVEALS_df) <- gsub(" [sd of cover in %]",
                          "_sd",
                          names(REVEALS_df),
                          fixed = TRUE)