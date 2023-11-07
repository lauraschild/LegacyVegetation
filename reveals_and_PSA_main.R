
#  in RStudio press ctrl+Shift+o to see the list of contents 

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# INITIALISATION ===================================================================================
{

    rm(list=ls()); gc()
    
    # set main path depending on where started the script from
    {
      main_dir <- "reveals_and_psas"
    }
    
    # source scripts ====
    {
        # functions and other scripts to source
        scripts_dir = file.path(main_dir, "source_files")
        source(file.path(scripts_dir,"reveals_and_PSA_functions.R"))
        
        # settings file
        source_file(file.path(main_dir,
                              "reveals_and_PSA_main_settings.R"))
    }
    
    # # back up program settings ====
    # {
    #     # create back up of the scripts that created the output data
    #     file.copy(file.path(main_dir, "reveals_and_PSA_main_settings.R"), 
    #               file.path(subset_dir, "Program", "settings.R"),
    #               copy.date = T)
    #     file.copy(file.path(main_dir, "reveals_and_PSA_main.R"), 
    #               file.path(subset_dir, "Program", "main.R"),
    #               copy.date = T)
    #     file.copy(file.path(main_dir, "source_files"), 
    #               file.path(subset_dir, "Program"),
    #               copy.date = T, recursive = T, overwrite = T)
    # }
}

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# PREPROCESSING ====================================================================================
{
    
    # load RPP and fallspeed parameters
    {
        cat("\nLoad PPE and fall speed parameters for set", rpp_set)
        
        parameters = readxl::read_xlsx(path = parameter_file_name,
                                       sheet = rpp_set)
        
        # REVEALS requires fallspeeds in the range (min_fallspeed, max_fallspeed)
        # set under/over to min_fallspeed/max_fallspeed
        parameters$Fallspeed[parameters$Fallspeed < min_fallspeed] = min_fallspeed
        parameters$Fallspeed[parameters$Fallspeed > max_fallspeed] = max_fallspeed
    }
    
    # load meta data ====
    {
        cat("\nLoad meta data for records.")
        if (update_meta_df){
            meta_df = update_meta_data(meta_df_columns)
        } else {
            # when running from cluster load the new df from disk (thus update before on local pc and changes have been made!)
            meta_df = suppressWarnings(read.csv(file=meta_df_file, stringsAsFactors = F, sep="\t"))
        }
        meta_df = modify_meta_df_for_reveals(meta_df)
    }
    
    # load continental pollen data frames ====
    {
        cat("\nLoad continental pollen data for", paste0(run_continents, "."))
        pollen_dfs = load_contintental_pollen_frames("Europe")
        meta_cols_pollen = pollen_dfs[[2]]
        pollen_dfs = pollen_dfs[[1]] 
    }
    
    # filter subset to run from meta and pollen data ====
    {
        meta_df = filter_subset(meta_df, subset_def)
        dataset_IDs = meta_df[[2]]
        meta_df = meta_df[[1]]
        
        meta_df$Basin_Radius_km = sqrt(meta_df$Basin_Area/pi)/1000
    }
    
#     # estimate computation times - so far only based on basin radius and median number of ages
#     {
#         if (output_folder == 'Comp_Time_Test'){
#             
#             # show number of samples per site
#             hist(meta_df$Number_Pollen_Samples, 
#                  breaks = seq(min(meta_df$Number_Pollen_Samples, na.rm = T),
#                               max(meta_df$Number_Pollen_Samples, na.rm = T)+9, 10))
#             median(meta_df$Number_Pollen_Samples, na.rm = T)
#             mean(meta_df$Number_Pollen_Samples, na.rm = T)
#             
#             # select sites to measure computation times, based on the distribution of basin radii
#             percentiles = quantile(meta_df$Basin_Radius_km, c(seq(0, .009, .001),
#                                                               seq(.01, .5, .01),
#                                                               seq(.55, .95, .05),
#                                                               seq(.955, 1, 0.001)))
#             percentiles[length(percentiles)] = percentiles[length(percentiles)] + 1
#             meta_all = meta_df %>% select(Dataset_ID, Number_Pollen_Samples, Basin_Radius_km)
#             for (jj in 1:(length(percentiles)-1)){
#                 measure_site = meta_all  %>% 
#                     filter((Basin_Radius_km >= percentiles[jj]) &
#                            (Basin_Radius_km <= percentiles[jj + 1]))
#                 meta_all = meta_all %>% filter(!(Dataset_ID %in% measure_site$Dataset_ID))
#                                                
#                 measure_site = measure_site %>%
#                     mutate(diff = abs(Number_Pollen_Samples - median(Number_Pollen_Samples, na.rm = T))) %>%
#                     arrange(diff) 
#                 measure_site = measure_site %>% slice(1) %>% select(Dataset_ID)
#                 
#                 if (jj == 1){
#                     measure_sites = c(as.integer(measure_site))
#                 } else{
#                     measure_sites = c(measure_sites, as.integer(measure_site))
#                 }
#             }
#             measure_sites = measure_sites[!is.na(measure_sites)]
#             
#             head(meta_df %>% filter(Dataset_ID %in% measure_sites)  %>% arrange(Basin_Radius_km))
#             tail(meta_df %>% filter(Dataset_ID %in% measure_sites)  %>% arrange(Basin_Radius_km))
#             
#             meta_df = meta_df %>% filter(Dataset_ID %in% measure_sites)  %>% arrange(Basin_Radius_km)
#             dataset_IDs = meta_df$Dataset_ID
#             
#             hist(meta_df$Basin_Radius_km)
#             
#         } else{
#             
#         
#             lookup_path = file.path(main_dir, "output", "Comp_Time_Test", "Tables", 
#                                     "Comp_Time_Test_comp_time_lookup.rda")
#             
#             if (file.exists(lookup_path)){
#                 # load lookup table for computation times
#                 comp_time_df = load_rda(lookup_path)
#                 
#                 
#                 # Basin radius and number of sites for current job
#                 {
#                     breaks = c(0)
#                     for (ii in 1:(nrow(comp_time_df)-1)){
#                         breaks = c(breaks, 
#                                    (comp_time_df$Basin_Radius_km[ii] + comp_time_df$Basin_Radius_km[ii+1])/2)
#                     }
#                     breaks = c(breaks, ceiling(max(comp_time_df$Basin_Radius_km)))
#                     
#                     hist_radii = hist(meta_df$Basin_Radius_km, plot = F, 
#                                       breaks = breaks) 
#                     # c(seq(0, 100, 10), 
#                     #        seq(150, round(max(meta_df$Basin_Radius_km)+49, -2), 50)))
#                     hist_radii = data.frame(cbind(Basin_Radius_km = hist_radii$breaks[2:length(hist_radii$breaks)], #hist_radii$mids, 
#                                                   Nr_Sites = hist_radii$counts)) %>% 
#                         filter(Nr_Sites != 0)   
#                 }
#                 
#                 # estimate computation times from lookup table previously computed
#                 {
#                     # so far for measured for sites with approx 50 ages 
#                     # and neglecting the impact of the number of taxa present at site with its individual combination of fallspeeds
#                     
#                     # optimistic estimation
#                     opt_inds = .bincode(hist_radii$Basin_Radius_km, comp_time_df$Basin_Radius_km)
#                     if (is.na(opt_inds[length(opt_inds)])){ opt_inds[length(opt_inds)] = nrow(comp_time_df)}
#                     opt_comp_times = round(sum(hist_radii$Nr_Sites * comp_time_df$Comp_Time_Min_With_PSA[opt_inds],
#                                                na.rm = T) / 60, 0)
#                     cat("\nOptimistic estimation of comp. times: ", opt_comp_times, " hours (", 
#                         round(opt_comp_times/24, 1), " days) for single core run.")
#                     
#                     # pessimistic estimation
#                     cons_inds = .bincode(hist_radii$Basin_Radius_km, comp_time_df$Basin_Radius_km) + 1
#                     if (is.na(cons_inds[length(cons_inds)])){ cons_inds[length(cons_inds)] = nrow(comp_time_df)}
#                     cons_comp_times = round(sum(hist_radii$Nr_Sites * comp_time_df$Comp_Time_Min_With_PSA[cons_inds],
#                                                 na.rm = T) / 60, 0)
#                     cat("\nVery pessimistic estimation of comp. times: ", cons_comp_times, " hours (", 
#                         round(cons_comp_times/24, 1), " days) for single core run.")
#                 }
#             }
#         }
#     }
}

# ///////////////////////////////////////////////////////////////////////////////////////////// ====
# REVEALS and PSA ==================================================================================
{
    if (length(dataset_IDs)>0){
        
        # # definition of loggers
        # {
        #     parallelLogger_available = require("ParallelLogger")
        #     if (parallelLogger_available){
        #         lib("ParallelLogger")
        #         logger_definitions = list("parallel" = list("overwrite_log_file" = T,
        #                                                     "path"=log_file),
        #                                   "IDs_to_start" = list("overwrite_log_file" = T,
        #                                                         "path"=log_file_IDs_to_start),
        #                                   "IDs_started" = list("overwrite_log_file" = T,
        #                                                        "path"=log_file_IDs_started),
        #                                   "IDs_finished" = list("overwrite_log_file" = T,
        #                                                         "path"=log_file_IDs_finished))
        #         list_of_loggers = create_loggers(logger_definitions)
        #         log_msg(logger_name = "IDs_to_start",
        #                 msg = paste0(dataset_IDs, collapse = " "))
        #     }
        # }
        
        if (parallel & (length(dataset_IDs) > 1)){
            
            ## parallel computations of sites ====
            cat("\nStart computations per site:\n")
            
            # definitions of export variables
            {
                
                export_variables = c(
                    
                    # static settings:
                    "overwrite_results", "compute_PSA", "age_col_new",
                    "plot_errors", "regionCutoff",
                    "rpp_shortname", "rpp_set", "dwm_shortname",
                    "dwm", "continents", "model_setting",
                    "n_ppe_variations",
                    "region_steps", "influx_thresh",
                    "main_dir", "output_folder", "subset_dir",
                    "repl_val_NA_PPE", "repl_val_NA_PPE_error",
                    "min_fallspeed", "max_fallspeed", "peat_area",
                    "plot_errors", "for_southern_sites_use_params_from",
                    "use_fixed_steps_lake_model", "n_steps_lake_model",
                    "scripts_dir", "compare_k_factors", "large_basins_last",
                    
                    # variables:
                    "dataset_IDs", "pollen_dfs", "parameters", "meta_df", "meta_cols_pollen",
                    c("dataset_IDs", "list_of_loggers")[parallelLogger_available+1],
                    
                    # functions (alternatively, delete single functions here and
                    #            import where necessary in the functions run parallely)
                    
                        # helper 
                        "lib", "source_file", "return_taxa", "replace_spaces", 
                        "replace_underscores", "clear", "is.scalar", 
                        
                        # parallel computing
                        "log_msg",
                        
                        # preprocessing
                        "load_contintental_pollen_frames", "update_meta_data", 
                        "modify_meta_df_for_reveals", "filter_subset",
                        
                        # reveals and psa
                        "checkREVEALS", "gpm", "LakeModel", 
                        "DispersalFactorK", "DispersalFactorK_nrSteps_fixed",
                        "rnorm_reveals", "REVEALS", "rmultinom_reveals", 
                        "run_reveals_with_PSA", "compute_PSR", "plotREVEALS_mod",
                        "reveals_analysis_per_site" 
                )
            }
            
            # run jobs
            {
                parallel_results = run_parallel(parallel_IDs = dataset_IDs,
                                                job_function = reveals_analysis_per_site,
                                                export_variables = export_variables,
                                                list_of_loggers = list_of_loggers,
                                                nr_parallel_max = nr_parallel_max)
            }
        } else{
            
            # sequential computation of sites ====
            for (ID in dataset_IDs){
                cat("\n\n\n", which(dataset_IDs == ID), "/", 
                    length(dataset_IDs))
                reveals_analysis_per_site(ID)
            }
        }

        # compare csv files found with IDs that should have been started ====
        {
            lib("readr")
            all_csvs = list.files(file.path(subset_dir, "Tables"), recursive = T)
            all_csvs = all_csvs[grep(".csv", all_csvs)]
            all_csvs = unique(sort(parse_number(all_csvs)))
            
            to_run = sort(meta_df$Dataset_ID)
            missing_sites = setdiff(to_run, all_csvs)
            cat("csv files missing in output files:\n\t", paste0(missing_sites, ", "))
            cat("Meta data of missing sites:")
            meta_df %>% filter(Dataset_ID %in% missing_sites)
        }
        
        # in test runs, save comp. times for future estimations
        {
            if (output_folder == 'Comp_Time_Test'){
                for (ID in dataset_IDs){
                    
                    comp_time_site = load_rda(file.path(subset_dir, "Tables", "Comp_times_per_site", 
                                                        paste0("site_", ID, 
                                                               "_computation_time_min_with_PSA_",
                                                               compute_PSA, ".rda")))
                    
                    if (ID == dataset_IDs[1]){
                        comp_times = comp_time_site
                    } else{
                        comp_times = rbind(comp_times, comp_time_site)
                    }
                }
                
                if (compute_PSA){ 
                    comp_times = comp_times %>% select(-Including_PSA) %>% rename(Comp_Time_Min_With_PSA = Comp_Time_Min, 
                                                                                  Dataset_ID = ID)
                } else{
                    comp_times = comp_times %>% select(-Including_PSA) %>% rename(Comp_Time_Min_only_REVEALS = Comp_Time_Min, 
                                                                                  Dataset_ID = ID)
                }
                meta_df = left_join(meta_df, comp_times)
                
                comp_time_df = meta_df %>% select(Dataset_ID, Number_Pollen_Samples, Basin_Radius_km, Comp_Time_Min_With_PSA) 
                comp_time_df
                
                # save such that after this test run with a few sites (output_folder = 'Comp_Time_Test' in settings) 
                # in future runs with many sites, the computation time can be roughly estimated beforehand 
                save(comp_time_df, file = file.path(subset_dir, "Tables", 
                                                    paste0(output_folder, "_comp_time_lookup.rda")))
            }
        }
    }
}
