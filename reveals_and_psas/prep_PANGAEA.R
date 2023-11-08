#script to make REVEALS runs ready for PANGAEA upload

rm(list = ls())

library(tidyverse)

main_dir <- "reveals_and_psas"

meta <- read.csv(file.path(main_dir,"input/meta_df_from_googlesheet.txt"),
                 sep = "\t")

#function to run on one REVEALS run (run dir and final name as input)
prep_PANGAEA <- function(run_dir,
                         name_to = run_dir){
  
  run_path <- paste0(main_dir,"/output/", run_dir,"/Tables")
  continents <- list.dirs(run_path,
                          full.names = FALSE,
                          recursive = FALSE)
  
  continents <- continents[continents %in% c("Africa","Asia","Europe","North_America",
                                             "South_America","Indopacific")]
  
  for (continent in continents){
    sites <- list.files(paste0(main_dir,"/output/",run_dir,"/Tables/",continent))
    
    if(length(sites) == 0) next
    
    #load pollen data 
    pollen_file <- "reveals_and_psas/input/LegacyPollen2_counts_"
    

    
    if(continent == "North_America") {
      pollen_df <- lapply(c("north_america_east","north_america_west"),
                          function(x) read_tsv(paste0(pollen_file,x,".csv"))%>%
                            mutate(Sample_ID = as.character(Sample_ID)))%>%
        bind_rows()
    }else{
      pollen_df <- read.csv(paste0(pollen_file,tolower(continent),".csv"),sep = "\t")%>%
        mutate(Sample_ID = as.character(Sample_ID))
    }
    meta_cols <- c("Sample_ID", "Event", "Data_Source","Site_ID", "Dataset_ID", "Longitude",
                   "Latitude", "Age_Source", "minAgeBP","maxAgeBP",
                   "medianAgeBP", "meanAgeBP", "Depth")

    # combine sites into continental dfs + add necessary metadata

    cont_df <- lapply(sites,
                      function(x) data.frame(data.table::fread(file.path(run_path,
                                                                         continent,
                                                                         x)))%>%
                        mutate(Sample_ID = as.character(Sample_ID),
                               Basin_Area = as.integer(Basin_Area)))%>%
      bind_rows()%>%
      merge(pollen_df[,meta_cols],
            by = c("Dataset_ID","Sample_ID"),
            all.x = TRUE)%>%
      relocate(meta_cols[-1])%>%
      select(-Sample_ID, -Distance_weighting, -RPP_Set, -Basin_Source)

    # rename Taxa columns
    new_names <- stringi::stri_replace_all_regex(names(cont_df),
                                                  pattern = c("[.]","q90", "q10"),
                                                  replacement = c(" [", "90_percentile", "10_percentile"),
                                                  vectorize = FALSE)
    
    new_names[grepl("[[]", new_names)] <- paste0(new_names[grepl("[[]",new_names)]," of cover in %]")
    
    names(cont_df) <- new_names
    
    #rename other columns
    cont_df <- cont_df %>%
      select(-any_of("Age_BP")) %>%
      rename(Pollen_Data_Source = Data_Source,
             Age_Model_Source = Age_Source,
             `Depth [cm]` = Depth,
             `Age_min [yrs BP]` = minAgeBP,
             `Age_max [yrs BP]` = maxAgeBP,
             `Age_median [yrs BP]` = medianAgeBP,
             `Age_mean [yrs BP]` = meanAgeBP,
             `Basin_Area [m^2]` = Basin_Area,
             `Basin_Diameter [m]` = Basin_Diameter,
             `Pollen_Source_Radius [m]` = Pollen_Source_Radius,
             `Pollen_Source_Area [m^2]` = Pollen_Source_Area,
             `relative_influx_at_PollenSourceRadius [median in %]` = Median_rel_influx_at_PSR,
             `relative_influx_at_PollenSourceRadius [standard deviation in %]` = Std_rel_influx_at_PSR)
    
    last_meta_col <- which(names(cont_df)== "relative_influx_at_PollenSourceRadius [standard deviation in %]")
    
    #alphabetically order taxa cols
    taxa_df <- cont_df[,(order(names(cont_df[,(last_meta_col+1):ncol(cont_df)]))+last_meta_col)]
    #change NAs in taxa cols to 0
    taxa_df[is.na(taxa_df)] <- 0
    
    cont_df <- cbind(cont_df[,1:last_meta_col],
                     taxa_df)
    PANGAEA_path <- "output/PANGAEA"
    if(!dir.exists(PANGAEA_path)) dir.create(PANGAEA_path)
    
    data.table::fwrite(cont_df,
                      paste0(PANGAEA_path,"/",
                             continent,"_",name_to,".csv"),
                      row.names = FALSE)

  }
  
}

prep_PANGAEA(run_dir = "original_RPP",
             name_to = "original_RPP")


