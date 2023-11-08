#load opti complete
filepath <- paste0("output/PANGAEA/optimized_REVEALS_",continent,".csv")

opti_df <- data.table::fread(filepath)%>%
  as.data.frame()%>%
  dplyr::select(any_of(names(pollen_df[1:13])),
         contains("[mean"))

names(opti_df) <- gsub(" [mean cover in %]",
                          "",
                          names(opti_df),
                          fixed = TRUE)
