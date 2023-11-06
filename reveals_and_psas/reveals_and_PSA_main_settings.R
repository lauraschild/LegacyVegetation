
#  in RStudio press ctrl+Shift+o to see the list of contents 

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# RUN SETTINGS =====================================================================================
{
    ## preprocessing settings ====
    {
        # update google sheet meta data before preprocessing data, works only from pc not cluster so far
        # meta data is never copied back to google table.
        # last updated: 12.04.23 (Peter)
        update_meta_df = F
    }
    
    # subset definitions ====
    output_folder = 'original_RPP'
    {
        # define a subset by creating a logical filter over rows in the meta data table as to be found in
        # https://docs.google.com/spreadsheets/d/1-_0PGA0Cqpil2Z6a2Pt3c_2yFfNM28XlA-j8Dj5gJN8/edit#gid=527521502
        # filter continents in run_continents below as this is also used for other purposes than filtering meta data
        # subset definitions
        {
            if (output_folder == 'original_RPP'){
                subset_def = str2expression(text="!is.na(meta_df$Basin_Area)")
            }
            
            if (output_folder == 'Comp_Time_Test'){
                subset_def = str2expression(text="!is.na(meta_df$Basin_Area)")
            }
            
            if (output_folder == "params_test"){
              subset_def = str2expression(text = "meta_df$Dataset_ID %in% c(41602, 42673, 45097)")
            }
            
          
          if (is.null(subset_def)){ stop(paste0("Subset, ", output_folder, " is not defined.")) }
          
        }
        
        skip_sites = c() # e.g. to exclude problematic dataset_IDs
        run_continents = c(
                           #"Asia"
                           #,"North America"
                           #,
                           "Europe"
                           #,
                           # "South America"
                           # ,"Indopacific"
                           # ,"Africa"
                           )  
    }
    
    # reveals settings ====
    reveals_analysis = T
    
    {
      #picking parameters
      parameter_set = "orig"  #"opti" or "orig"
      hemisphere = "northern" # "northern" or "southern"
      
      rpp_set = paste0("NH_", parameter_set)
      if(hemisphere == "northern") rpp_set = paste0("CF_",parameter_set)
      
      rpp_shortname = rpp_set
      for_southern_sites_use_params_from = "Northern Hemisphere"

        
        parallel = T # NEEDS TO BE UPDATED. run parallel on cluster, remember to update google sheet locally before, if changes has been made to it
            nr_parallel_max = 20
        large_basins_last = T # for psa computation a few very large lakes may take many hours, put those to the end of ID list to get most results quick
        overwrite_results = F # set F for resuming computations, set T to make a fresh version of all sites to this subset folder
        n_ppe_variations = 2000 # 2000 for final run, not relevant for pollen source radius
        
        compare_k_factors = F # see folder subset/Plots/compare_step_sizes_lake_model
        use_fixed_steps_lake_model = T # modified version by Peter as Theuerkaufs parametrization produced unnecessary many steps for large lakes
            n_steps_lake_model = 500 # 500 for final runs # tested that this is enough steps between basin center and basin radius such that influx integration over basin area does not deviate more than 1% from Theuerkaufs solution
    }
        
    # pollen source area settings ====
    compute_PSA = T
    {
        region_steps = 40 # 40 for final runs. go from regionCutoff to basin radius in region_steps steps, break at 80% influx from regioncutoff (~= infinity)
        regionCutoff = 1e+06 # m
        influx_thresh = 80 # %
    }
}

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# GENERAL DEFINITONS ===============================================================================
{
    ## reveals parameter ====
    {
        run_other_as_peatland = T 
        
        dwm = "gpm unstable"
        dwm_shortname = "US"
        
        peat_area = 100 
        repl_val_NA_PPE = 1
        repl_val_NA_PPE_error  = 0.25
        min_fallspeed = 0.01 # m/s # required by Theuerkaufs REVEALS version. Lower values in our dataset are set to this minimum value.
        max_fallspeed = 0.15 # m/s # required by Theuerkaufs REVEALS version. Higher values in our dataset are set to this minimum value.
        
        plot_errors = "sd"
    }
    
    # basin types ====
    {
        basin_types  = 
            list("Peatland" = c("Peat","Peatland", "Bog","Fen",
                                "Marsh","Mire","Swamp",
                                "Spring","Creek","Wetland",
                                "Stream","p","P"),
                 "Lake" = c("Lake","Pond","Paleo-lake","l","L"),
                 "Other" = c("River","r", "R",
                             "Regularly flooded",
                             "Fluvial","Terrace",
                             "Alluvial", "Alluvial fans", "Alluvial/Fluvial",
                             "Ocean","o","O","Marine","m","M",
                             "Estuary","Fjord","Lagoon",
                             "Others", "Other", "other", "Marl", "Farmland",
                             "Soil", "Cliff", "Travertine",
                             "Coast", "Loess", "Loess-Paleosol",  "Loess/Palaeosol",
                             "Valley", "Cave", "Watershed", "Animal Midden",
                             "Meadow", "Forest","Dune","Paleosol", "loess/palaeosol", 
                             "Archaeological", "Archaeological profile",
                             "Delta", "?", "", "Ice cap", "Ice Cap", "Ice"))
    }
    
    # other ====
    {
        continents = c("Asia",
                       "North America",
                       "Europe",
                       "South America",
                       "Indopacific",
                       "Africa") # the order must stay like this
        
        age_cols_old = c("meanAgeBP", "mean_Age"); age_col_new = "Age_BP"
        model_setting = paste0(rpp_shortname,"_",dwm_shortname)
        
        meta_cols_pollen = c("Event", "Data_Source", "Site_ID", "Dataset_ID", "Site_Name", 
                             "Site_Name_Clean", "Longitude", "Latitude", "Continent",
                             "Age_Source", "minAgeBP", "maxAgeBP",         
                             "medianAgeBP", "meanAgeBP", "Depth",            
                             "Sample_ID", "Age_AD_BC", "Age_BP")
    }
}


# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# PATH definitions ===========================================================================================
{
    input_dir = file.path(main_dir, "input")
    output_dir = file.path(main_dir, "output")
    subset_dir = file.path(output_dir, output_folder)
    
    ## pollen data ====
    {
        # according to Thomas Böhmers mail from 28.10.2022:
        # bioing\data\LegacyPollen\Datasets_assigned_ages\Counts
        # also see 'pollen data information.txt' in the input folder
        # date of last copy to REVEALS folder: 08.04.23 (Peter)
        pollen_asia_file = file.path(input_dir, "LegacyPollen2_counts_asia.csv")
        pollen_europe_file = file.path(input_dir, "LegacyPollen2_counts_europe.csv")
        pollen_north_america_east_file = file.path(input_dir, "LegacyPollen2_counts_north_america_east.csv")
        pollen_north_america_west_file = file.path(input_dir, "LegacyPollen2_counts_north_america_west.csv")
        pollen_south_america_file = file.path(input_dir, "LegacyPollen2_counts_south_america.csv")
        pollen_indopacific_file = file.path(input_dir, "LegacyPollen2_counts_indopacific.csv")
        pollen_africa_file = file.path(input_dir, "LegacyPollen2_counts_africa.csv")
    }
    
    # google sheet with meta data ====
    {
        meta_df_link = 'https://docs.google.com/spreadsheets/d/1-_0PGA0Cqpil2Z6a2Pt3c_2yFfNM28XlA-j8Dj5gJN8/edit#gid=527521502'
        meta_df_file = file.path(input_dir, "meta_df_from_googlesheet.txt")
        meta_df_columns = c("Project", "Event", "Data_Source", 
                            "Site_Name", "Site_Name_Clean",
                            "Site_ID", "Dataset_ID", "Latitude", "Longitude", 
                            "Elevation", "Continent",
                            "Archive_Type", "Basin_Source", "Basin_Comment", 
                            "Basin_Area", "Data_Type", "Number_Pollen_Samples", "AgeModel_Source", 
                            "Number_Datings")
    }
    
    # basin areas ====
    {
        basin_areas_file = file.path(input_dir, "220923_basin_areas_all.csv")
    }
    
    # rpp and fall speed data ====
    {
        parameter_file_name = file.path(input_dir, "parameters.xlsx")
    }
    
    # log files
    {
        log_file = file.path(main_dir, "output", output_folder, "Logs",
                                    paste0(output_folder,"_logging.txt"))
        log_file_IDs_to_start = file.path(main_dir, "output", output_folder, "Logs",
                                          paste0(output_folder, "_IDs_to_start.log"))
        log_file_IDs_started = file.path(main_dir, "output", output_folder, "Logs",
                                         paste0(output_folder, "_IDs_started.log"))
        log_file_IDs_finished = file.path(main_dir, "output", output_folder, "Logs",
                                          paste0(output_folder, "_IDs_finished.log"))
    }
    
}  

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# PREPARE FOLDERS ==================================================================================
{
    ## definition of (sub) folders ====
    
    subfolders = list("Plots" = c("site_overviews",
                                  "compare_step_sizes_lake_model",
                                  "REVEALS_results"=list(continents),
                                  "REVEALS_results_and_pollen"=list(continents)),
                      "Tables" = c(continents, "Comp_times_per_site"),
                      "Program" = c(),
                      "Logs" = c())
    
    # create folders ====
    {
        if (file.exists(subset_dir)==F){ dir.create(subset_dir, recursive = TRUE) }
        cat(paste0("\n\nCreate output folders for subset \'", output_folder, 
                   "\' if necessary:"))
        
        create_subfolders(subfolders, subset_dir)
    }
    
}
