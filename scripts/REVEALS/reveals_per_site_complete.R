#reveals_analysis_per_site complete

reveals_analysis_per_site_complete = function(ID,
                                     params){
  
  library(dplyr)
  
  #get parameters
  source("scripts/REVEALS/REVEALSinR_modified.R", local = TRUE)
  source("scripts/REVEALS/reveals_and_PSA_functions.R", local = TRUE)
  source("scripts/REVEALS/reveals_settings.R",  local = TRUE)

  # prepare tables for site
  pollen_site <- pollen_df %>%
    filter(Dataset_ID == ID)
  
  #replace spaces with underscores
  names(pollen_site)[-(1:13)] <- gsub(" ","_",names(pollen_site)[-(1:13)])
  
  #replace NA values with 0
  pollen_site[is.na(pollen_site)] <- 0
  #remove empty columns 
  #except for age!!
  pollen_site <- cbind(pollen_site[,1:13],
                       pollen_site[,-(1:13)][,colSums(pollen_site[,-(1:13)])!=0])
  
  # adapt to REVEALS conventions
  tBasin <- filter(meta_df, Dataset_ID == ID)%>%
    pull(Basin_Type)
  if(tBasin == "other") tBasin <- "peatland"
  ABasin <- filter(meta_df,Dataset_ID == ID)%>%
    pull(Basin_Area)
  dBasin <- round(sqrt((4*ABasin)/pi))
  
  #drop age from pollen_site and add sample_ID
  pollen_reveals <- pollen_site %>%
    select(any_of(taxa))%>%
    mutate(sample_ID = row_number())%>%
    relocate(sample_ID)
  
  #run REVEALS
  results <- run_reveals_with_PSA(pollen = pollen_reveals,
                                  params = params,
                                  tBasin = tBasin,
                                  dBasin = dBasin,
                                  dwm = dwm,
                                  n = n_ppe_variations,
                                  regionCutoff = regionCutoff,
                                  ppefun = ppefun,
                                  pollenfun = pollenfun,
                                  verbose = FALSE,
                                  ID = ID,
                                  influx_thresh = NULL,
                                  compute_psr = FALSE)
  
  #organize results
  REVEALS_site <- results$reveals_results
  REVEALS_site <- cbind(pollen_site[,1:13],
                        REVEALS_site)%>% 
    select(-sample_ID,-Basin.type, -Distance.weighting)
   
  names(REVEALS_site)[14:ncol(REVEALS_site)] <-stringi::stri_replace_all_regex(names(REVEALS_site)[14:ncol(REVEALS_site)],
                                                                 pattern = c(".mean",".median",".q90",".q10",".sd"),
                                                                 replacement = c(" [mean cover in %]"," [median cover in %]",
                                                                                 " [10_percentile of cover in %]", " [90_percentile of cover in %]",
                                                                                 " [sd of cover in %]"),
                                                                 vectorize_all = FALSE)
  
  
  return(REVEALS_site)
}


