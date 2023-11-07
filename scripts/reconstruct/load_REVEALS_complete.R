#load REVEALS complete
#script to read REVEALS data


path <- paste0("reveals_and_psas/output/PANGAEA/",continent,"_original_RPP.csv")

REVEALS_df <- data.table::fread(path)%>%
  as.data.frame()%>%
  dplyr::select(any_of(names(pollen_df[1:13])),
         contains("[mean of cover"))

names(REVEALS_df) <- gsub(" [mean of cover in %]",
                          "",
                          names(REVEALS_df),
                          fixed = TRUE)
