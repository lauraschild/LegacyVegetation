#reveals_analysis_per_site

reveals_analysis_per_site = function(ID,
                                     params){
  
  library(dplyr)
  
  #get parameters
  source("scripts/REVEALS/reveals_settings.R", echo=FALSE)
  source("scripts/REVEALS/REVEALSinR_modified.R", echo=FALSE)
  source("scripts/REVEALS/reveals_and_PSA_functions.R", echo = FALSE)

  # prepare tables for site
  pollen_site <- pollen_df %>%
    filter(Dataset_ID == ID)
  
  #replace spaces with underscores
  names(pollen_site) <- gsub(" ","_",names(pollen_site))
      
  #replace NA values with 0
  pollen_site[is.na(pollen_site)] <- 0
  #remove empty columns 
  #except for age!!
  pollen_site <- cbind(pollen_site[,1:4],
                       pollen_site[,-(1:4)][,colSums(pollen_site[,-(1:4)])!=0])
 
  # adapt to REVEALS conventions
  tBasin <- filter(meta_df, Dataset_ID == ID)%>%
    pull(Basin_Type)
  if(tBasin == "other") tBasin <- "peatland"
  ABasin <- filter(meta_df, Dataset_ID == ID)%>%
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
  REVEALS_site <- cbind(pollen_site[,1:4],
                        REVEALS_site[,grepl(".mean",names(REVEALS_site))])
  names(REVEALS_site) <- gsub(".mean","",names(REVEALS_site))
  
  return(REVEALS_site)
}
        
      
     