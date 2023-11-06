
# HELPER FUNCTIONS
{
    load_rda = function(file_name){
        # loads an rda file, and returns it
        # consider saving and loading rds otherwise, where loading into a variable works 
        # without this workaround
        
        load(file_name)
        get(ls()[ls() != "file_name"]) 
    }
    
    create_subfolders = function(x, folder_path, indent=1) {
        for (ii in 1:(length(x))){
            folders = x[[ii]]
            if (is.list(x)){
                folder_name = replace_spaces(names(x)[ii])
                if(folder_name == ""){folder_name = as.character(x[ii])}
                
            } else if ((is.vector(folders) == T) & (length(folders) == 1)) {
                folders = x[ii]
                folder_name = replace_spaces(as.character(x[ii]))
            } 
            cat("\n",paste0(rep("\t",indent,collapse=T)), folder_name)
            final_path = file.path(folder_path, replace_spaces(folder_name))
            if (file.exists(final_path)==F){
                dir.create(final_path)
            }
            folders
            if ((typeof(folders) == "list") | ((is.vector(folders) == T) & (length(folders) > 1))){
                create_subfolders(folders, file.path(folder_path, folder_name), indent = indent+1)
            } else if (is.list(folders) & (length(folders[[1]])==1)){
                create_subfolders(folders, file.path(folder_path, folder_name), indent = indent+1)
            }
        }
        cat("\n")
    }
    
    lib = function(lib_name){
        for (name in lib_name){
            suppressPackageStartupMessages(library(name, character.only = T))
        }
    }
    
    source_file = function(file_path){
        suppressMessages(suppressWarnings(suppressPackageStartupMessages(
            source(file_path))))
    }
    
    return_taxa = function(df,meta){
        # return all cols of df but those specified in meta
        taxa_colnames = list()
        for (jj in colnames(df)){
            if (jj %in% meta==F){
                taxa_colnames = append(taxa_colnames,jj)
            }
        }
        return(unlist(taxa_colnames))
    }
    
    replace_spaces = function(x) { gsub(" ","_",x) }
    
    replace_underscores = function(x) { gsub("_"," ",x) }
    
    clear = function(){
        # clear plots and opened but not closed plot files prior to to new plots
        if (as.numeric(dev.cur()!=1)){
            for (ii in 1:(as.numeric(dev.cur())-1)){
                dev.off()
            }
        }
    }
    
    is.scalar = function(x) is.atomic(x) && length(x) == 1L
}

# PREPROCESSING
{
    load_contintental_pollen_frames = function(cont_names = NULL){
        
        # load all continents if non were passed
        if (is.null(cont_names)==T){cont_names = continents}
        
        # load only specified continents if passed
        {
            # asia
            if ("Asia" %in% cont_names){
                pollen_asia = read.csv(pollen_asia_file,
                                       sep = "\t", stringsAsFactors = F)
                if (any(age_cols_old %in% colnames(pollen_asia))){
                    old_ind = which(age_cols_old %in% colnames(pollen_asia))
                    pollen_asia[,age_col_new] = as.numeric(pollen_asia[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_asia) == T){
                    pollen_asia[,age_col_new] = as.numeric(pollen_asia[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_asia) == T){
                    pollen_asia[,"Sample_ID"] = as.character(pollen_asia[,"Sample_ID"])
                } 
            } else{
                pollen_asia = NULL
            }
            
            # north america (combine east and west first)
            if ("North America" %in% cont_names){
                pollen_north_america_east = read.csv(pollen_north_america_east_file,
                                                     sep = "\t", stringsAsFactors = F)
                pollen_north_america_west = read.csv(pollen_north_america_west_file,
                                                     sep = "\t", stringsAsFactors = F)
                pollen_north_america = suppressMessages(dplyr::full_join(pollen_north_america_east,
                                                                         pollen_north_america_west))
                pollen_north_america$Continent = "North America"
                if (any(age_cols_old %in% colnames(pollen_north_america))){
                    old_ind = which(age_cols_old %in% colnames(pollen_north_america))
                    pollen_north_america[,age_col_new] = as.numeric(pollen_north_america[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_north_america) == T){
                    pollen_north_america[,age_col_new] = as.numeric(pollen_north_america[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_north_america) == T){
                    pollen_north_america[,"Sample_ID"] = as.character(pollen_north_america[,"Sample_ID"])
                }   
            } else{
                pollen_north_america = NULL
            }
            
            # europe
            if ("Europe" %in% cont_names){
                pollen_europe = read.csv(pollen_europe_file,
                                         sep = "\t", stringsAsFactors = F)
                if (any(age_cols_old %in% colnames(pollen_europe))){
                    old_ind = which(age_cols_old %in% colnames(pollen_europe))
                    pollen_europe[,age_col_new] = as.numeric(pollen_europe[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_europe) == T){
                    pollen_europe[,age_col_new] = as.numeric(pollen_europe[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_europe) == T){
                    pollen_europe[,"Sample_ID"] = as.character(pollen_europe[,"Sample_ID"])
                }   
            } else{
                pollen_europe = NULL
            }
            
            # south america
            if ("South America" %in% cont_names){
                pollen_south_america = read.csv(pollen_south_america_file,
                                                sep = "\t", stringsAsFactors = F)
                if (any(age_cols_old %in% colnames(pollen_south_america))){
                    old_ind = which(age_cols_old %in% colnames(pollen_south_america))
                    pollen_south_america[,age_col_new] = as.numeric(pollen_south_america[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_south_america) == T){
                    pollen_south_america[,age_col_new] = as.numeric(pollen_south_america[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_south_america) == T){
                    pollen_south_america[,"Sample_ID"] = as.character(pollen_south_america[,"Sample_ID"])
                }   
            } else{
                pollen_south_america = NULL
            }
            
            # indopacific
            if ("Indopacific" %in% cont_names){
                pollen_indopacific = read.csv(pollen_indopacific_file,
                                              sep = "\t", stringsAsFactors = F)
                if (any(age_cols_old %in% colnames(pollen_indopacific))){
                    old_ind = which(age_cols_old %in% colnames(pollen_indopacific))
                    pollen_indopacific[,age_col_new] = as.numeric(pollen_indopacific[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_indopacific) == T){
                    pollen_indopacific[,age_col_new] = as.numeric(pollen_indopacific[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_indopacific) == T){
                    pollen_indopacific[,"Sample_ID"] = as.character(pollen_indopacific[,"Sample_ID"])
                }   
            } else{
                pollen_indopacific = NULL
            }
            
            # africa
            if ("Africa" %in% cont_names){
                pollen_africa = read.csv(pollen_africa_file,
                                         sep = "\t", stringsAsFactors = F)
                if (any(age_cols_old %in% colnames(pollen_africa))){
                    old_ind = which(age_cols_old %in% colnames(pollen_africa))
                    pollen_africa[,age_col_new] = as.numeric(pollen_africa[,age_cols_old[old_ind]])
                }
                if (age_col_new %in% colnames(pollen_africa) == T){
                    pollen_africa[,age_col_new] = as.numeric(pollen_africa[,age_col_new])
                } 
                if ("Sample_ID" %in% colnames(pollen_africa) == T){
                    pollen_africa[,"Sample_ID"] = as.character(pollen_africa[,"Sample_ID"])
                }   
            } else{
                pollen_africa = NULL
            }
        }
        
        #     
        pollen_dfs = list("Asia"=pollen_asia, "North_America"=pollen_north_america,
                          "Europe"=pollen_europe,"South_America"=pollen_south_america,
                          "Indopacific"=pollen_indopacific, "Africa"=pollen_africa)
        
        meta_cols_pollen = meta_cols_pollen[meta_cols_pollen %in% 
                                            as.character(unlist(lapply(pollen_dfs, colnames)))]
        return(list("pollen_dfs"=pollen_dfs, "meta_cols_pollen"=meta_cols_pollen))
    }
    
    update_meta_data = function(cols){
        
        # writing back to google sheet is not implemented yet
        
        lib(c("dplyr","googlesheets4"))
        
        # start from pc with update_meta_new=T, 
        # before running on cluster with update_meta_new=F
        meta_df = data.frame(read_sheet(meta_df_link, na = "NA"), 
                             stringsAsFactors = F) %>% 
            select(-starts_with("Publication"), -Submission_Pollendata)
        meta_df$Site_ID = as.integer(meta_df$Site_ID)
        meta_df$Dataset_ID = as.integer(meta_df$Dataset_ID)
        meta_df$Basin_Comment = as.character(meta_df$Basin_Comment)
        
        # add areas from GIS computations
        basin_areas = read.csv(basin_areas_file, sep="\t") %>% 
            select(-d_m)
        colnames(basin_areas) = c("Dataset_ID", "Basin_Area")
        basin_areas = basin_areas %>% filter(Dataset_ID %in%  meta_df$Dataset_ID)
        
        # add basin areas from this file to meta table and 
        meta_df = meta_df %>% 
            left_join(basin_areas, by = "Dataset_ID")
        if ("Basin_Area.x" %in% colnames(meta_df)){
            meta_df = meta_df %>%
                mutate(Basin_Area = ifelse((Basin_Area.x!=Basin_Area.y) & (!is.na(Basin_Area.y)),
                       Basin_Area.y, Basin_Area.x)) %>%
            select(-c("Basin_Area.x", "Basin_Area.y"))
        }
        
        # write to file. from here the REVEALS script reads the new meta_df when it cant be updated on the cluster
        meta_df = meta_df[, cols[cols %in% colnames(meta_df)]]
        cols[cols=="Archive_Type"] = "Basin_Type"; colnames(meta_df) = cols
        write.table(meta_df, file=meta_df_file, row.names = F, sep = "\t")
        
        return(meta_df)
    }
    
    modify_meta_df_for_reveals = function(meta_df){
        # merge eastern and western north america to continent North America
        meta_df$Continent[which(meta_df$Continent %in% c("Western North America",
                                                         "Eastern North America"))] = "North America"
        
        # make sure all basin types in the tables (different spellings etc) are defined 
        # in settings. Only then the translation to REVEALS standard basin types works
        basin_types_in_meta = unique(meta_df$Basin_Type)
        all(basin_types_in_meta %in% as.character(unlist(basin_types)))
        meta_df$Basin_Type[which(meta_df$Basin_Type %in% basin_types$Peatland)] = "peatland"
        meta_df$Basin_Type[which(meta_df$Basin_Type %in% basin_types$Lake)] = "lake"
        meta_df$Basin_Type[which(meta_df$Basin_Type %in% basin_types$Other)] = "other"
        
        # check for duplicated IDs
        meta_df = unique(meta_df)
        dupl_ids_in_meta = meta_df$Dataset_ID[duplicated(meta_df$Dataset_ID)]
        del_inds = which((meta_df$Dataset_ID %in% dupl_ids_in_meta) & (meta_df$Project == "LegacyPollen_1.0"))
        if (length(del_inds) > 0){meta_df = meta_df[-del_inds,]}
        dupl_ids_in_meta = meta_df$Dataset_ID[duplicated(meta_df$Dataset_ID)]
        if (length(dupl_ids_in_meta) > 0){stop("Duplicated Dataset_ID in meta data table.")}
        
        global_lake_area_median = median(meta_df$Basin_Area[which(meta_df$Basin_Type == "lake")], na.rm=T)
        #global_lake_radius_median = sqrt(global_lake_area_median / pi)
        # european_lake_area_median = median(meta_df$Basin_Area[which((meta_df$Basin_Type == "lake") &
        #                                                                 (meta_df$Continent == "Europe"))], na.rm=T)
        # european_lake_radius_median = sqrt(european_lake_area_median / pi)
        # asian_lake_area_median = median(meta_df$Basin_Area[which((meta_df$Basin_Type == "lake") & 
        #                                                              (meta_df$Continent == "Asia"))], na.rm=T)
        # asian_lake_radius_median = sqrt(asian_lake_area_median / pi)
        # north_american_lake_area_median = median(meta_df$Basin_Area[which((meta_df$Basin_Type == "lake") & (meta_df$Continent == "North America"))], na.rm=T)
        # north_american_lake_radius_median = sqrt(north_american_lake_area_median / pi)
        # cat("\nGlobal vs. continental median lake radius in meter for asia, north america, europe:\n\t",
        #     global_lake_radius_median, "vs.", asian_lake_radius_median, ", ", north_american_lake_radius_median, ", ", european_lake_radius_median)
        
        # replace all missing lake and other areas with global median  
        meta_df$Basin_Area[which(is.na(meta_df$Basin_Area) & 
                                     (meta_df$Basin_Type %in% c("lake", "other")))] = global_lake_area_median
        length(which(is.na((meta_df$Basin_Area) & (meta_df$Basin_Type == "lake")))) == 0
        length(which(is.na((meta_df$Basin_Area) & (meta_df$Basin_Type == "other")))) == 0
        
        # replace all missing peatland  areas with fixed 100 m^2 value  
        meta_df$Basin_Area[which(is.na(meta_df$Basin_Area) & (meta_df$Basin_Type == "peatland"))] = peat_area
        length(which(is.na((meta_df$Basin_Area) & (meta_df$Basin_Type == "peatland")))) == 0
        
        return(meta_df)
    }
    
    filter_subset = function(meta_df, subset_def){
        
        lib("dplyr")
        
        if (run_other_as_peatland){
            meta_df$Basin_Type[which(meta_df$Basin_Type == "other")] = "peatland"
            meta_df$Basin_Area[which(meta_df$Basin_Type == "peatland")] = peat_area
        }
        
        # filter for subset defined in settings
        inds_to_run = which(eval(subset_def))
        subset_IDs = meta_df[inds_to_run,]    
        subset_IDs = subset_IDs$Dataset_ID[which(subset_IDs$Continent %in% run_continents)]
        length(subset_IDs)
        
        # find Dataset_IDs in pollen data and defined in subset
        pollen_IDs = unique(unlist(lapply(pollen_dfs, function(u){ u$Dataset_ID })))
        dataset_IDs = intersect(pollen_IDs, subset_IDs)
        cat("\nFiltered", length(dataset_IDs), "dataset_IDs for subset", 
            paste0(output_folder, ".\n"))
        
        if (length(skip_sites)>=1){
            dataset_IDs = setdiff(dataset_IDs, skip_sites)
        }
        
        meta_df = meta_df %>% filter(Dataset_ID %in% dataset_IDs)
        
        # compute huge basins (with ~ 100 - 400 km basin radius) at last 
        # to get the majority of sites done
        if (large_basins_last == T){
            meta_df = meta_df %>% arrange(Basin_Area)
            dataset_IDs = meta_df$Dataset_ID
        }
        
        return(list(meta_df, dataset_IDs))
    }
}

# PARALLEL COMPUTATION
{   
    create_loggers = function(logger_list){
        
        lib(c("ParallelLogger"))
        
        clearLoggers()
        for (ii in 1:length(logger_list)){
            
            logger_name = names(logger_list)[ii]
            
            logfile_path = logger_list[[ii]]$path
            log_level = logger_list[[ii]]$log_level
            if (is.null(log_level)){ log_level = "INFO" }
            overwrite = logger_list[[ii]]$overwrite_log_file
            if (is.null(overwrite)){ overwrite = F }
            
            if (overwrite){ close(file(logfile_path, open="w")) }
            logger_list[[ii]] = createLogger(name = logger_name,
                                             threshold = log_level, # i.e. log infos and all messages with higher priorities, so don't log logTrace and logDebug
                                             appenders = list(createFileAppender(layout = layoutParallel,
                                                                                 fileName = logfile_path,
                                                                                 overwrite = overwrite)))
        }
        return(logger_list)
    }
    
    run_parallel = function(parallel_IDs, 
                            job_function, 
                            export_variables, 
                            list_of_loggers = NULL, 
                            nr_parallel_max = 2){
        
        # find nr. of parallel jobs
        {
            cores_avail = parallel::detectCores()
            cat(paste0("\n\n\t",cores_avail," cores available. ",
                       "Limit number of parallel jobs to ", cores_avail-1),".")
            cat("\n\tNr. of jobs: ", length(parallel_IDs))
            cat("\n\tMax. nr. of cores defined in settings: ", nr_parallel_max)
            nr_parallel = min(cores_avail-1, nr_parallel_max)
            nr_parallel = min(length(parallel_IDs), nr_parallel)
            ind = as.integer(length(parallel_IDs)>nr_parallel)+1
            cat(paste0("\n\tReserving ",nr_parallel," cores ", 
                       "since ",c("only ","more than ")[ind],
                       c(length(parallel_IDs), nr_parallel)[ind],
                       " job(s) to do.\n"))
        }
        
        if (nr_parallel == 1){
            
            # run sequentially if nr of parallel jobs is 1
            {
                cat("\n\tCompute sites sequentially (nr_parallel = 1):")
                
                sequential_results = list()
                for (ID in IDs){
                    cat("\t\tPreprocessing for ID", ID)
                    ID_ind = which(IDs == ID)
                    sequential_results[[ID_ind]] = job_function(ID)
                }
            }
            
            return(sequential_results)
            
        } else{
            
            # initialize cluster depending on packages installed
            {
                parallelLogger_available = require("ParallelLogger")
                if (parallelLogger_available){
                    
                    lib(c("ParallelLogger", "snow"))
                    
                    if (length(parallel_IDs) > 1){
                        
                        cluster = ParallelLogger::makeCluster(nr_parallel)
                        
                        # export variables to parallel nodes that are relevant to execute job_function
                        {
                            # variables created in global environment, i.e. created or sourced in main script
                            export_variables = unique(export_variables)
                            
                            if (!is.null(list_of_loggers)){
                                snow::clusterExport(cluster, 
                                                    c("list_of_loggers",
                                                      export_variables), 
                                                    envir = environment())
                            } else{
                                snow::clusterExport(cluster, export_variables)
                            }
                        }
                        
                        # parallel computation with load balancing
                        parallel_results = ParallelLogger::clusterApply(cluster = cluster,
                                                                        x = parallel_IDs,
                                                                        fun = job_function,
                                                                        stopOnError = T)
                        ParallelLogger::stopCluster(cluster)
                        
                    } else if (length(parallel_IDs) == 1){
                        
                        parallel_results = job_function(parallel_IDs) #, list_of_loggers)
                    } else { 
                        
                        stop("Parallel computation: No parallel IDs passed.")
                        
                    }   
                } else {
                    
                    lib(c("parallel", "doSNOW", "snow"))
                    
                    # run with parallel package (without logging)
                    doSNOW::registerDoSNOW(cluster)
                    cluster = makeCluster(nr_parallel, type = "SOCK")
                    snow::clusterExport(cluster, export_variables)
                    # parallel computation with load balancing
                    parallel_results = parallel::clusterApplyLB(cluster,
                                                                parallel_IDs,
                                                                job_function)
                    stopCluster(cluster)
                    
                }
            }
            
            return(parallel_results)
        }
    }
    
    log_msg = function(logger_name, msg){
        parallelLogger_available = require("ParallelLogger")
        if (parallelLogger_available){
            if (!is.null(list_of_loggers)){
                lib("ParallelLogger")
                
                logger_names = names(list_of_loggers)
                logger = list_of_loggers[[which(logger_names == logger_name)]]
                
                if (class(logger) == "Logger"){
                    registerLogger(logger)
                    logInfo(msg)
                    success = unregisterLogger(logger, silent = T)
                    
                    if (!success){ cat(paste0("\n\t\tLogging message: '", msg, 
                                              "' failed for logger: ", logger_name))}
                } else {
                    cat("\nCheck definitions of loggers. Class(logger) != 'Logger'")
                }
            }
        }
    }
}

# REVEALS FUNCTIONS (partially modified from REVEALSinR 0.9.7)
{
    reveals_analysis_per_site = function(ID){
        
        lib(c("tictoc","tidyr", "dplyr"))
        
        # prepare job
        {
            tic.clear(); tic()
            
            # log_msg(logger_name = "IDs_started",
            #         msg = paste0(" -- ", ID, " -- "))
            
            meta_names = c("Dataset_ID", "Site_Name", "Project",
                           "Continent", "Latitude", "Longitude",
                           "Basin_Type", "Basin_Area", "Basin_Radius",
                           "Data_Source", "Project")
            meta_site = meta_df %>% filter(Dataset_ID == ID) %>% 
                mutate(Basin_Radius = sqrt(Basin_Area / pi))

            # extract and print information on ID from meta data table
            {
                project = meta_site$Project
                site_name = meta_site$Site_Name
                lat = meta_site$Latitude
                lon = meta_site$Longitude
                basin_source = meta_site$Basin_Source
                basin_type = meta_site$Basin_Type
                basin_diameter = meta_site$Basin_Radius *2
                basin_area = meta_site$Basin_Area
                continent = meta_site$Continent
                
                cat("\n\n\t------- Dataset_ID:", ID, 
                    "\n\t\tContinent:", continent, 
                    "\n\t\tProject:", project, 
                    "\n\t\tDataset name:", site_name, 
                    "\n\t\t(Lat, Lon):", paste0("(", lat, ", ", lon, ")"), "\n")
            }
            
            # check if ID's results already saved on disk. for saving time 
            {
                #   when overwrite_results, compute ID and save it to disk
                #   when dont overwrite_results, only compute ID when its new (not on disk yet)
                
                output_path = file.path(subset_dir, "Tables", replace_spaces(continent),
                                        paste0("site_", ID, "_reconstr_rel_veg_cover_timeseries_",
                                               model_setting,".csv"))
                
                # if not overwriting existing results, check if file already on disk
                # (when run_site == T will overwrite the site's results on disk)
                if (overwrite_results==F){
                    table_exist = file.exists(output_path)
                    if (table_exist==F){
                        run_site = T
                    } else {
                        run_site = F
                    }
                } else {
                    run_site = T
                }
            }
        }
        
        if (run_site){
            
            # extract continental pollen data frame
            pollen_continent = pollen_dfs[[replace_spaces(continent)]]
            
            # prepare tables for site
            {
                # subset columns and drop Dataset_ID
                {
                    pollen_site = pollen_continent %>% 
                        filter(Dataset_ID == ID) %>% 
                        select(-any_of(meta_cols_pollen[!(meta_cols_pollen %in%
                                                              c("Sample_ID", age_col_new, "Continent"))]))
                    taxa_cols = return_taxa(pollen_site, meta_cols_pollen)
                    cols = c("Sample_ID", age_col_new, "Continent", taxa_cols)
                    pollen_site = pollen_site[,cols]
                    row.names(pollen_site) = NULL
                    
                    # if any spaces in taxa names, replace them with _ 
                    # since reveals does not allow spaces
                    colnames(pollen_site) = sapply(colnames(pollen_site),replace_spaces)                        
                }
                
                # replace NAs with 0 and remove all zero columns and rows
                {
                    # in samples (age slices) with all pollen counts equal zero
                    # reveals produces an error due to dividing by zero in
                    # randomizing error function. That's why we remove all-zero lines
                    pollen_site = pollen_site %>% 
                        mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% 
                        filter(rowSums(is.na(.)) == 0)
                    
                    # also remove all-zero cols (i.e. taxa not present at this ID)
                    pollen_site = cbind(pollen_site %>% select(where(function(x) !is.numeric(x))),
                                        pollen_site %>% select(where(is.numeric)) %>% 
                                            select(where(function(x) sum(x) != 0)))
                }
                
                # take the subset of parameters of taxa for which also pollen count data are available
                {
                    params_site = parameters %>% select(-ends_with("_Source")) %>% 
                        filter(Taxon %in% colnames(pollen_site)) %>%
                        filter(Continent %in% c(continent, rpp_set))
                    
                    if (length(unique(params_site$Taxon)) != 
                       length(params_site$Taxon)){
                        # use continental values at first, 
                        # then those filled up from Northern Hemisphere set
                        stop("Check filtering parameters for site.")
                    }
                }
                
                # prepare reveals input and add missing parameter to site's parameters frame
                {
                    # drop continent columns
                    {
                        pollen_reveals = pollen_site %>% select(-Continent)
                        params_reveals = params_site %>% select(-Continent)
                    }
                    
                    cols = colnames(pollen_site)
                    meta_cols = c("Sample_ID", "Continent","Age_BP", "Age_AD_BC")
                    taxa_cols = return_taxa(pollen_site, meta_cols)
                    
                    # add taxa with 1 +/- 0.25 to parameter frame, where we have missing values and mean fallspeed 
                    {
                        if (continent %in% c("Asia", "North America", "Europe")){
                            mean_fallspeed = parameters %>% filter(Continent %in% c(continent, "Northern Hemisphere")) %>%
                                summarize(x = mean(Fallspeed)) %>% pull(x)
                        } else {
                            mean_fallspeed = parameters %>% filter(Continent == for_southern_sites_use_params_from) %>%
                                summarize(x = mean(Fallspeed)) %>% pull(x)
                        }
                        
                        missing_params = setdiff(colnames(pollen_reveals[, taxa_cols]),
                                                 params_reveals$Taxon) 
    
                        params_reveals = rbind(params_reveals, 
                                               do.call(rbind, lapply(missing_params,
                                                                     function(u){
                                                                         data.frame("Taxon" = u,
                                                                                    "Fallspeed" = mean_fallspeed,
                                                                                    "PPE" = repl_val_NA_PPE,
                                                                                    "PPE_Error" = repl_val_NA_PPE_error,
                                                                                    stringsAsFactors = F)
                                                                     })))
                        taxa_cols = return_taxa(pollen_reveals, meta_cols)
                        if (!all(sort(params_reveals$Taxon) == sort(taxa_cols))){
                            stop("Preparing REVEALS input:\n",
                                 "\n\tIn PPE and fallspeed parameter table for ID ", ID, 
                                 "\n\tdifferent taxa defined than in pollen data for ID ", ID)
                        }
                    }
                    
                    # adapt to REVEALS conventions
                    {
                        tBasin = basin_type
                        dBasin = as.integer(basin_diameter)
                        verbose = F
                        n = n_ppe_variations
                        params_reveals = params_reveals %>% rename(species = Taxon,
                                                                   fallspeed = Fallspeed, 
                                                                   PPE = PPE, 
                                                                   PPE.error = PPE_Error)
                    }
                }
                
                # drop age and run REVEALS with sample_ID
                {
                    pollen_reveals = pollen_reveals %>% select(-any_of(c("Age_BP", "Age_AD_BC")))
                }
            }
            
            # REVEALS and PSA analysis per site
            {
                # run reveals and PSA computation
                log_msg(logger_name = "parallel",
                        msg = paste0("Dataset_ID ", ID,
                                     ": Collected pollen data, meta data, and PPE and fallspeed parameter.",
                                     " Start computations."))
                pollenfun = rmultinom_reveals
                ppefun = rnorm_reveals
                if (compute_PSA == T){
                    cat(paste0("\n\t\tCompute Pollen Source Area and REVEALS cover."))
                    log_msg(logger_name = "parallel",
                            msg = paste0("Dataset_ID ", ID, ": ", 
                                                "Compute Pollen Source Area with ",
                                                influx_thresh, "% of total pollen influx",
                                                " and REVEALS cover."))
                    
                    reveals_results = run_reveals_with_PSA(pollen=pollen_reveals,
                                                           params=params_reveals,
                                                           tBasin = tBasin, 
                                                           dBasin = dBasin, 
                                                           dwm = dwm, 
                                                           n = n, 
                                                           regionCutoff = regionCutoff, 
                                                           ppefun = ppefun, 
                                                           pollenfun = pollenfun,
                                                           verbose = verbose,
                                                           ID = ID, 
                                                           influx_thresh = influx_thresh, 
                                                           compute_psr = T)
                    
                    # extract the results
                    reveals_output = reveals_results$reveals_results
                    source_radius = reveals_results$source_radius
                    source_area = reveals_results$source_area
                    median_rel_influx = reveals_results$median_rel_influx
                    std_rel_influx = reveals_results$std_rel_influx
                    rm(reveals_results)
                    
                } else{
                    # orginal versions of REVEALSinR (v. 0.9.7)
                    
                    log_msg(logger_name = "parallel",
                            msg = paste0("Dataset_ID ", ID, ": ", 
                                         "Compute REVEALS cover."))
                    
                    ppefun = rnorm_reveals
                    reveals_results = run_reveals_with_PSA(pollen = pollen_reveals,
                                                           params = params_reveals,
                                                           tBasin, dBasin, 
                                                           dwm, 
                                                           n, 
                                                           regionCutoff, 
                                                           ppefun, 
                                                           pollenfun,
                                                           verbose,
                                                           ID = ID, 
                                                           influx_thresh=NULL, 
                                                           compute_psr = F)
                    reveals_output = reveals_results$reveals_results
                }
            }
            
            # add information to reveals output 
            {
                if (!all(reveals_output$Sample_ID == pollen_site$Sample_ID)){
                    stop("Order or number of samples changed in REVEALS run.")
                }
                reveals_output$Dataset_ID = ID
                reveals_output$Continent = continent
                reveals_output[,age_col_new] = pollen_site$Age_BP
                reveals_output$Basin_Area = basin_area
                reveals_output$Basin_Source = basin_source
                reveals_output$RPP_Set = rpp_set
                reveals_output$Basin_Diameter = basin_diameter
                meta_cols = c(meta_cols_pollen[which(meta_cols_pollen %in% colnames(reveals_output))],
                              "Distance.weighting",
                              "Basin.type",
                              "Basin_Area",
                              "Basin_Diameter",
                              "RPP_Set",
                              "Basin_Source")
                taxa_cols = return_taxa(reveals_output, meta_cols)
                cols = c(meta_cols, taxa_cols)
                reveals_output = reveals_output[,cols]
                meta_cols = gsub(".", "_", meta_cols, fixed = T)
                cols = c(meta_cols, taxa_cols)
                colnames(reveals_output) = cols
                
                if (compute_PSA){
                    reveals_output$Pollen_Source_Radius = source_radius
                    reveals_output$Pollen_Source_Area = source_area
                    reveals_output$Median_rel_influx_at_PSR = median_rel_influx
                    reveals_output$Std_rel_influx_at_PSR = std_rel_influx
                    
                } else {
                    reveals_output$Pollen_Source_Radius = NA
                    reveals_output$Pollen_Source_Area = NA
                    reveals_output$Median_rel_influx_at_PSR = NA
                    reveals_output$Std_rel_influx_at_PSR = NA
                }
                meta_cols = c(meta_cols, 
                              "Pollen_Source_Radius", "Pollen_Source_Area",
                              "Median_rel_influx_at_PSR", "Std_rel_influx_at_PSR")
                cols = c(meta_cols,taxa_cols)
                reveals_output = reveals_output[,cols]
            }    
            
            # save per ID results and finish job
            {
                write.table(reveals_output, file = output_path, 
                            sep="\t", row.names = F)
                
                timer = toc(quiet=T)
                duration = round((as.numeric(timer$toc) - as.numeric(timer$tic))/60,2)
                
                cat("\n\n\t\tDuration:", duration, "minutes.")
            
                log_msg(logger_name = "IDs_finished",
                        msg = paste0(" -- ", ID, " -- ", "Computed."))
                
                log_msg(logger_name = "parallel",
                        msg = paste0("Dataset_ID ", ID, ": ", 
                                     "Computation done. ", duration, " min elapsed."))
                comp_time_min = data.frame(ID = ID, Comp_Time_Min = duration, Including_PSA = compute_PSA) 
                save(comp_time_min,
                     file = file.path(subset_dir, "Tables", "Comp_times_per_site",
                                      paste0("site_", ID, "_computation_time_min_with_PSA_", compute_PSA, ".rda")))
            }
            
        } else {
            cat(paste0("\n\t\tDataset_ID ", ID, ": File exists already,\n\t\tskipped site."))
            log_msg(logger_name = "IDs_finished",
                    msg = paste0(" -- ", ID, " -- ", "Skipped."))
            
            stop = toc(quiet = T)
            log_msg(logger_name = "parallel",
                    msg = paste0("Dataset_ID ", ID, ": ", 
                                 "Site skipped."))
        }
    }
    
    checkREVEALS = function (pollen, params, tBasin, dBasin, dwm, n, regionCutoff, verbose){
        if ((dBasin/2) > regionCutoff) 
            stop("Stop: Basin larger than the size of the region (regionCutoff) - that does'nt work!")
        if (!all(colnames(pollen[, -1]) %in% rownames(params), na.rm = FALSE)) 
            stop("STOP: Not all pollen taxa found in the parameter list")
        if (!isTRUE((dwm == "lsm unstable") | (dwm == "gpm neutral") | 
                    (dwm == "gpm unstable") | (dwm == "1overd"))) 
            stop("distance weighting method not defined; should be 'LSM unstable', 'GPM neutral', 'GPM unstable' or '1oved' ")
        if (!isTRUE((tBasin == "peatland") | (tBasin == "lake"))) 
            stop("basin type ('tBasin') not defined; should be 'peatland' or 'lake'")
        if (!isTRUE(dBasin == floor(dBasin))) 
            stop("basin size ('dBasin') should be an integer value")
        if (dBasin < 10) 
            stop("basin diameter ('dBasin') should be at least 10 m")
        if (!isTRUE(n == floor(n))) 
            stop("'n' should be an integer value")
        # if (n < 1000) 
        #     message("Warning: for sensible error estimates, 'n' should be at least 1000")
        if (!"fallspeed" %in% tolower(names(params))) 
            stop("STOP: no 'fallspeed' column in the parameters ")
        if (max(params$fallspeed) > 0.15) 
            stop("fallspeed(s) too high (>0.15 m/s), please check")
        if (min(params$fallspeed) < 0.01) 
            stop("fallspeed(s) too low (<0.01 m/s), please check")
        if (!"ppe" %in% tolower(names(params))) 
            stop("STOP: no 'ppe' column in the parameters ")
        if (!"ppe.error" %in% tolower(names(params))) 
            stop("STOP: no 'ppe.error' column in the parameters ")
        return(TRUE)
    }
    
    gpm = function (fallSpeed, cz, n, u){
        x = seq(from = 0, to = 100, by = 0.1)
        x = x^3
        y = exp((-4 * fallSpeed * x^(n/2))/(n * u * sqrt(pi) * cz))
        return(loess(y ~ x, span = 0.015, model = FALSE, degree = 2))
    }
    
    LakeModel = function (x, dBasin, disModel, stepSize, regionCutoff){
        rBasin = dBasin/2
        r2 = seq(from = rBasin - x, to = 2 * rBasin, by = 2)
        r2Cut = r2
        r2Cut[r2Cut > (rBasin + x)] = rBasin + x
        alpha = acos((r2Cut^2 + x^2 - rBasin^2)/(2 * r2Cut * x))
        propOut = 1 - alpha/pi
        a = diff(propOut)/2
        propOut = propOut[-1] - a
        airborne = predict(disModel, r2)
        influxRing = abs(diff(airborne))
        influxRing = influxRing * propOut
        influx = sum(influxRing) + predict(disModel, 2 * rBasin) - 
            predict(disModel, regionCutoff)
        ringArea = pi * x^2 - pi * (x - stepSize)^2
        return(influx * ringArea)
    }
    
    DispersalFactorK = function (vg, tBasin, dBasin, dwm, regionCutoff, u = 3){
        rBasin = dBasin/2
        if (dwm == "gpm neutral") {
            disModel = gpm(fallSpeed = vg, cz = 0.12, n = 0.25, 
                            u = u)
        }
        else if (dwm == "gpm unstable") {
            disModel = gpm(fallSpeed = vg, cz = 0.21, n = 0.2, u = u)
        }
        else if (dwm == "1overd") {
            disModel = OneOverD.f()
        }
        else if (dwm == "lsm unstable low") {
            vgr = 200 * round(vg, 2) + 1
            vgr[vgr == 1] = 2
            disModel = lsm_unstable_low[[vgr]]
        }
        else if (dwm == "lsm unstable") {
            vgr = 100 * round(vg, 2) + 1
            disModel = lsm_unstable[[vgr]]
        }
        else stop("no valid distance weighting method selected")
        if (tBasin == "peatland") {
            if (regionCutoff <= 1e+06) 
                influx = predict(disModel, rBasin) - predict(disModel, 
                                                              regionCutoff)
            else influx = predict(disModel, rBasin)
        }
        else if (tBasin == "lake") {
            stepSize = 2
            if (dBasin >= 2000) 
                stepSize = 10
            rSteps = seq(from = stepSize, to = rBasin - 1, by = stepSize)
            lakeInflux = do.call(rbind, lapply(rSteps, LakeModel, 
                                                dBasin = dBasin, disModel = disModel, 
                                                stepSize = stepSize, 
                                                regionCutoff = regionCutoff))
            influx = sum(lakeInflux)/(pi * (max(rSteps))^2)
        }
        return(influx)
    }
    
    DispersalFactorK_nrSteps_fixed = function (vg, tBasin, dBasin, dwm, regionCutoff, u = 3,
                                                nrSteps){
        rBasin = dBasin/2
        if (dwm == "gpm neutral") {
            disModel = gpm(fallSpeed = vg, cz = 0.12, n = 0.25, 
                            u = u)
        }
        else if (dwm == "gpm unstable") {
            disModel = gpm(fallSpeed = vg, cz = 0.21, n = 0.2, u = u)
        }
        else if (dwm == "1overd") {
            disModel = OneOverD.f()
        }
        else if (dwm == "lsm unstable low") {
            vgr = 200 * round(vg, 2) + 1
            vgr[vgr == 1] = 2
            disModel = lsm_unstable_low[[vgr]]
        }
        else if (dwm == "lsm unstable") {
            vgr = 100 * round(vg, 2) + 1
            disModel = lsm_unstable[[vgr]]
        }
        else stop("no valid distance weighting method selected")
        if (tBasin == "peatland") {
            if (regionCutoff <= 1e+06) 
                influx = predict(disModel, rBasin) - predict(disModel, 
                                                              regionCutoff)
            else influx = predict(disModel, rBasin)
        }
        else if (tBasin == "lake") {
            stepSize = round(rBasin / nrSteps, 0)
            rSteps = seq(from = stepSize, to = (rBasin-1), by = stepSize)
            lakeInflux = do.call(rbind, lapply(rSteps, LakeModel, 
                                                dBasin = dBasin, disModel = disModel, 
                                                stepSize = stepSize, 
                                                regionCutoff = regionCutoff))
            influx = sum(lakeInflux)/(pi * (max(rSteps))^2)
        }
        return(influx)
    }
    
    rnorm_reveals = function (n, mean = rep(1, n), sds = rep(1, n)){
        x = 0 * mean - 1
        while (any(x < 0)) {
            x = rnorm(n = n, mean = mean, sd = sds)
        }
        return(x)
    }
    
    REVEALS = function(pollen, n = n, ppes, ppefun = disqover::rnorm_reveals, 
                       deposition, pollenfun = disqover::rmultinom_reveals, 
                       verbose = verbose){
        
        #require("ParallelLogger")
        
        if (verbose) {
            cat(paste0("\n\t\t\tComputing sample\t", pollen[1]))
            #logInfo(paste("\n\t\t\tComputing sample ", pollen[1]))
        }
        coverM = replicate(n = n, {
            ppes[, 3] = ppefun(length(ppes[, "ppe"]), ppes[, "ppe"], ppes[, "ppe.error"])
            pollenPer = pollenfun(n = 1, pollen = as.numeric(pollen[-1]))
            disp_ppe = deposition * ppes[, 3]
            s = sum(pollenPer/disp_ppe)
            round(100 * pollenPer/(disp_ppe * s), digits = 3)
        }, simplify = TRUE)
        rownames(coverM) = rownames(pollen[-1])
        coverMean = apply(coverM, 1, mean)
        coverMedian = apply(coverM, 1, median)
        coverSD = apply(coverM, 1, sd)
        coverq10 = apply(coverM, 1, function(u) quantile(probs = 0.1, 
                                                          x = u))
        coverq90 = apply(coverM, 1, function(u) quantile(probs = 0.9, 
                                                          x = u))
        return(invisible(list(meansim = coverMean, sdsim = coverSD, 
                              mediansim = coverMedian, q90sim = coverq90, q10sim = coverq10)))
    }
    
    rmultinom_reveals = function (n = 1, pollen, ...) {
        100 * rmultinom(n, sum(pollen), pollen/sum(pollen))/sum(pollen, na.rm = TRUE)
    }
    
    run_reveals_with_PSA = function(pollen, params, tBasin, dBasin, dwm = "lsm unstable", 
                                    n = 1000, regionCutoff = 1e+05, ppefun = rnorm_reveals, 
                                    pollenfun = rmultinom_reveals, verbose = TRUE, 
                                    ID, influx_thresh=80, compute_psr=F){
        
        # prepare input
        {
            dwm = tolower(dwm)
            tBasin = tolower(tBasin)
            names(params) = tolower(names(params))
            pollen = as.data.frame(pollen)
            params = as.data.frame(params)
            row.names(params) = params[, 1]
            params = params[, -1]
            if (!checkREVEALS(pollen = pollen, params = params, tBasin = tBasin, 
                              dBasin = dBasin, dwm = dwm, n = n, 
                              regionCutoff = regionCutoff, 
                              verbose = verbose)) 
                stop("Some parameter is not suitable - check and correct")
            params = params[colnames(pollen[, -1]), ]
            ppes = params[, c("ppe", "ppe.error")]
            vg = params[, "fallspeed"]
        }
        
        if (compute_psr){
            # compute pollen source radius in meter 
            # as distance from site center with influx_thresh percent of total pollen income
            psa_object = compute_PSR(vg, params,
                                     influx_thresh, regionCutoff, 
                                     dBasin, tBasin)
            psr = psa_object$psr
            psa = psa_object$psa
            median_rel_influx = psa_object$median_rel_influx
            std_rel_influx = psa_object$std_rel_influx
        }
        
        # examine how k-factors differ with only n_steps_lake_model in lake model
        {
            if (compare_k_factors == T){
                
                lib("tidyverse")
                lib("tictoc")
                tic.clear()
                
                # last step size servers as reference value, i.e. "true" k-factors from Martin Theuerkaufs in lakeModel of REVEALSinR 0.9.7
                rBasin = dBasin/2
                if (dBasin >= 2000){
                    ref_value = round(rBasin/ 10)
                } else {
                    ref_value = round(rBasin/ 2)
                }
                
                # check relative error of k-factors compared to those for ref_value
                test_steps = c(5, 10, 20, 50, seq(100, 500, 50))
                test_steps = test_steps[test_steps < ref_value]
                test_steps = c(test_steps, ref_value)
                
                k_runs = list(); ii=1
                for (test_step in test_steps){
                    tic(paste0("n_steps_", as.character(test_step)))
                    deposition = do.call(rbind, lapply(vg, DispersalFactorK_nrSteps_fixed, 
                                                       tBasin = tBasin, dBasin = dBasin,
                                                       dwm = dwm, regionCutoff = regionCutoff, 
                                                       nrSteps = test_step))
                    timer = toc()
                    duration = (as.numeric(timer[[2]]) - as.numeric(timer[[1]]))/60
                    if (ii==1){ times_min = duration } else { times_min = c(times_min, duration) }
                    k_runs[[ii]] = deposition
                    ii=ii+1
                }
                k_facs = data.frame(do.call(cbind, lapply(k_runs,
                                                          function(u){
                                                              u 
                                                          })), 
                                    stringsAsFactors = F)
                colnames(k_facs) = paste0("n_",as.character(test_steps))
                k_facs_rel = data.frame("taxon" = colnames(pollen)[2:ncol(pollen)], 
                                        apply(k_facs, 2, 
                                              function(u){ 
                                                  ref_values = k_facs[,ncol(k_facs)]
                                                  if (is.character(u)){u} else {(u - ref_values) * 100 / ref_values}
                                              }), stringsAsFactors = F)
                # plotting for only means over taxa
                {
                    k_facs_rel = apply(k_facs_rel[,2:ncol(k_facs_rel)], 2, mean)
                    
                    lib("ggrepel")
                    lib("dplyr")
                    lib("ggplot2")
                    lib("viridis")
                    lib("colorspace")
                    #hcl_palettes(plot = TRUE)
                    lib("shinyjs")
                    #hclwizard()
                    
                    plot_df = data.frame("n_steps" = test_steps,
                                         "comp_time" = factor(round(times_min,2)),
                                         "rel_error" = as.numeric(k_facs_rel), 
                                         stringsAsFactors = F)
                    red_beige = sequential_hcl(n = nrow(plot_df), h = c(-10, 31), 
                                               c = c(95, 75, 80), l = c(35, 69), 
                                               power = c(1.2, 1))
                    p = ggplot(plot_df, aes(x=n_steps, y=rel_error, color=comp_time)) + 
                        geom_line(color="dimgray", size=0.6) +
                        geom_point(size=3.5) +
                        geom_hline(yintercept=0, linetype="dashed", color = "dimgray", size=.5) +
                        geom_label_repel(aes(label = paste0(n_steps, " steps | ",
                                                            "error = ", round(rel_error,1)," %"),
                                             segment.color=comp_time),
                                         box.padding   = .6,
                                         point.padding = 2,
                                         show.legend=F,
                                         max.overlaps = 50) +
                        scale_color_manual(values = red_beige) + 
                        ggtitle(paste0("Relative error of mean K-factor compared to REVEALSinR 0.9.7",
                                       "and nr. of steps for mixing in lake model", 
                                       "\nDataset ID ", ID,
                                       " | Basin radius ", round(rBasin), " m",
                                       " | Nr. of taxa ", (ncol(pollen)-1),
                                       " | Nr. of samples ", nrow(pollen))) +
                        ylab("Relative deviation to mean reference K-value (%)") + 
                        xlab("Nr. of steps withi basin radius ( )") +
                        labs(color = "Computation time\nfor K-factors (min)") + 
                        theme_light() + 
                        guides(color = guide_legend(override.aes = list(size = 4)))+
                        theme(legend.box.background = element_rect(color="grey90", size=2.5),
                              legend.position = c(.925, .925),
                              legend.justification = c("right", "top"),
                              legend.box.just = "right",
                              legend.margin = margin(6, 6, 6, 6))
                    
                    comp_steps_plot_path = file.path(subset_dir, "Plots", 
                                                     "compare_step_sizes_lake_model",
                                                     paste0("rel_error_k_factors_site_", ID,".png"))
                    
                    clear()
                    png(comp_steps_plot_path, width=1067, height=600)
                    print(p)
                    dev.off()
                    
                }
                # {
                #     # for all taxa, not ready yet
                #     plot_df = k_facs_rel %>%
                #                 pivot_longer(!taxon, names_to = "nr_steps", values_to = "rel_error")
                #     plot_df$nr_steps = factor(plot_df$nr_steps)
                #     # plot_df$taxon = factor(plot_df$taxon)
                # }
                rm(deposition)
            }
        }
        
        # compute pollen deposition (i.e. pollen influx from regioncutoff)
        {
            cat("\n\n\t\tCompute reconstr. cover with REVEALS.")
            rBasin = dBasin/2
            if (rBasin >= 1000){
                nr_steps_orig = round(rBasin/10)
            } else {
                nr_steps_orig = round(rBasin/2)
            }
            
            if (use_fixed_steps_lake_model){
                # modified version of lake model
                # use minimum nr of steps (or max step size) as in REVEALS 0.9.7
                n_steps_use = min(n_steps_lake_model, nr_steps_orig)
                cat("\n\t\t\tCompute deposition in",  n_steps_use, 
                    paste0("steps\n\t\t\t(orig. in 0.9.7: ", nr_steps_orig," steps)."))
                deposition = do.call(rbind, lapply(vg, DispersalFactorK_nrSteps_fixed, 
                                                   tBasin = tBasin, dBasin = dBasin,
                                                   dwm = dwm, regionCutoff = regionCutoff, 
                                                   nrSteps = n_steps_use))
            } else {
                # original version of lake model from REVEALSinR 0.9.7
                cat("\n\t\tCompute deposition with", nr_steps_orig,"steps (REVEALSinR 0.9.7).")
                deposition = do.call(rbind, lapply(vg, DispersalFactorK, 
                                                    tBasin = tBasin, dBasin = dBasin,
                                                    dwm = dwm, regionCutoff = regionCutoff))
            }
        }
        
        # compare magnitudes of relative K factors with that of PPE as rough estimate
        # which has greater impact on reconstructed cover
        {
            relative_ks = expand.grid(deposition, deposition)
            relative_ks = relative_ks[,1] / relative_ks[,2]
            min_ratio_ks = min(relative_ks)
            max_ratio_ks = max(relative_ks)
            mean_ratio_ks = mean(relative_ks)
            sd_ratio_ks = sd(relative_ks)
            

            relative_ppes = expand.grid(params$ppe, params$ppe)
            relative_ppes = relative_ppes[,1] / relative_ppes[,2]
            min_ratio_ppes = min(relative_ppes)
            max_ratio_ppes = max(relative_ppes)
            mean_ratio_ppes = mean(relative_ppes)
            sd_ratio_ppes = sd(relative_ppes)
        }
        
        # compute cover per sample
        {
            all = apply(pollen, 1, REVEALS, n = n, ppes = ppes, ppefun = ppefun, 
                        deposition = deposition, pollenfun = pollenfun,
                        verbose = verbose)
            results_df = data.frame(dwm, tBasin, 
                                     pollen[1], t(sapply(all,
                                                         FUN = function(x) x$meansim)),
                                     t(sapply(all, FUN = function(x) x$mediansim)), 
                                     t(sapply(all, FUN = function(x) x$q90sim)), 
                                     t(sapply(all, FUN = function(x) x$q10sim)), 
                                     t(sapply(all, FUN = function(x) x$sdsim)))
            colnames(results_df) = c("Distance.weighting", 
                                      "Basin.type", 
                                      names(pollen)[1], 
                                      paste(names(pollen[-1]),
                                            rep(c("mean",
                                                  "median",
                                                  "q90",
                                                  "q10",
                                                  "sd"),
                                                each = length(pollen[1, ]) - 1),
                                            sep = "."))
        }
        
        # return cover and pollen source information
        {
            if (compute_psr){
                return(list("reveals_results" = results_df, 
                            "source_radius" = psr,
                            "source_area" = psa,
                            "median_rel_influx" = median_rel_influx,
                            "std_rel_influx" = std_rel_influx
                            ))
            } else {
                return(list("reveals_results" = results_df))
            }
        }
    }
    
    compute_PSR = function(vg, params,
                           influx_thresh, regionCutoff, 
                           dBasin, tBasin){
        
        # define steps with different region_cutoffs
        {
            rBasin = dBasin/2
            dist = (regionCutoff - rBasin) / 2
            
            # many shorter steps closer to basin where influx changes more rapidly
            inner_steps_nr = as.integer(region_steps*3/4)
            inner_steps_length = (dist) / (inner_steps_nr-1)
            dist = (regionCutoff - rBasin) / 2
            inner_steps = seq(from=rBasin, to=rBasin + dist,
                              by=inner_steps_length)
            
            # few larger steps for outer half
            outer_steps_nr = as.integer(region_steps/4)
            outer_steps_length = (dist - inner_steps_length) / (outer_steps_nr-1)
            outer_steps = seq(from=rBasin + dist + inner_steps_length, to=regionCutoff,
                              by = outer_steps_length) 
            
            cutoff_steps = c(inner_steps, outer_steps)
            # diff(cutoff_steps)
            # plot(cutoff_steps)
        }
        
        # find step lengths for lake model
        {
            if (rBasin >= 1000){
            nr_steps_orig = round(rBasin/10)
            } else {
                nr_steps_orig = round(rBasin/2)
            }
            n_steps_use = min(n_steps_lake_model, nr_steps_orig)
        }
        
        compute_influx_per_taxon = function(vg, cutoff, nr_steps, dBasin, tBasin){
            influx_per_taxon = do.call(rbind, 
                                       lapply(vg, DispersalFactorK_nrSteps_fixed, 
                                              tBasin = tBasin, dBasin = dBasin,
                                              dwm = dwm, regionCutoff = cutoff, 
                                              nrSteps = nr_steps))
        }
        # test that compute_influx_per_taxon function does the same as original with cutoff=regionCutoff 
        # influx_at_regionCutoff = compute_influx_per_taxon(vg, cutoff=cutoff_steps[length(cutoff_steps)], n_steps_use, dBasin, tBasin)
        # all(deposition == influx_at_regionCutoff)
        
        influx_cutoffs = data.frame(matrix(NA, nrow = length(vg), ncol=0))
        ref_influx = compute_influx_per_taxon(vg, cutoff = regionCutoff, 
                                              nr_steps = n_steps_use, 
                                              dBasin, tBasin)
        
        cat(paste0("\n\t\t\tVary region cutoff to find\n\t\t\t", influx_thresh, 
                   "% of total pollen influx:"))
        # logInfo(paste0("Vary region cutoff to find ", influx_thresh, 
        #                "% of total pollen influx."))
        for (cutoff in cutoff_steps[1:(length(cutoff_steps)-1)]){
            step_name = paste0("Step_", 
                               as.character(which(cutoff_steps == cutoff)))
            
            cur_influx = compute_influx_per_taxon(vg, cutoff = cutoff, 
                                                  nr_steps = n_steps_use, 
                                                  dBasin, tBasin)
            rel_influx = cur_influx / ref_influx
            median_rel_influx = median(rel_influx, na.rm=T) * 100
            std_rel_influx = sd(rel_influx, na.rm = T) * 100
            
            cat(paste0("\n\t\t\t\t", step_name, ". Within ", 
                       round(cutoff/1000,1),"km ", 
                round(median_rel_influx,2),"%"))
            
            if (median_rel_influx >= influx_thresh){
                
                # add current cutoff and reference value to table 
                influx_cutoffs[,step_name] = cur_influx
                influx_cutoffs[,"region_cutoff"] = ref_influx
                
                psr = cutoff
                
                cat(paste0("\n\t\t\tStop iteration and use pollen source radius",
                           "\n\t\t\tof ", round(cutoff/1000,1), 
                           "km with ", paste0(round(median_rel_influx,2), "% of total influx.")))
                psa = pi*cutoff^2
                break
            }
            influx_cutoffs[,replace_spaces(step_name)] = cur_influx
            
            influx_cutoffs_rel = apply(influx_cutoffs[,-ncol(influx_cutoffs)], 2, 
                                       function(u){ 
                                           u / influx_cutoffs[,ncol(influx_cutoffs)]
                                       }) * 100
        }
        
        # plot pollen influx radii
        {
            # library(dplyr, tidyr, ggplot2)
            # steps_psr =  c(cutoff_steps[cutoff_steps <= psr])
            # 
            # plot_df = data.frame(Taxon = taxa_cols, 
            #                      influx_cutoffs_rel)
            # #plot_df$region_cutoff = 100
            # plot_df$Taxon = factor(plot_df$Taxon)
            # colnames(plot_df)[grep("Step_", colnames(plot_df))] = steps_psr
            # #colnames(plot_df)[which(colnames(plot_df) == "region_cutoff")] = regionCutoff 
            # plot_df = plot_df %>% pivot_longer(cols = colnames(plot_df)[-1],
            #                                    names_to = "Radius", values_to = "Influx")
            # plot_df$Radius = as.numeric(plot_df$Radius)
            # plt = ggplot(plot_df, aes(x = Radius, y = Influx)) +
            #     geom_line(aes(colour = Taxon)) + 
            #     geom_vline(xintercept = psr, linetype = "dashed", 
            #                color = "dimgray", linewidth = 0.8) + 
            #     geom_hline(yintercept = influx_thresh, linetype = "dashed", 
            #                color = "dimgray", linewidth = 0.8) +
            #     theme_minimal() + 
            #     labs(x = "Radius from site center (m)",
            #          y = paste0("Relative Pollen influx compared ",
            #                     "to Max distance", regionCutoff, " m (%)"),
            #          title = paste0("Pollen Source Radius is distance with ", influx_thresh, 
            #                         "% of pollen influx (median over taxa)"),
            #          subtitle = paste0("Basin type: ", tBasin, 
            #                            " | Basin raidus: ", round(dBasin/2,0), " m"))
            # png(file.path(output_dir, output_folder, "Plots", "pollen_source_radii",
            #               paste0("dataset_ID_", ID, "_psr_and_influx_over_radius.png")),
            #     width = 1200, height = 800)
            # print(plt)
            # dev.off()
        }
        
        
        # compute psr and std of psr accuratly after all steps
        # (mean and std of individual cutoff per taxon
        # instead of mean relative influx over all taxa >= influx_thresh)
        # if implementing, remove break condition above. 
        # Will take a lot longer though. Maybe decrease region_steps then.
        {
            # # compute psr and std of psr
            # influx_cutoffs_rel = influx_cutoffs * 100 / influx_cutoffs[,"region_cutoff"]
            # cutoffs_per_taxon = cutoff_steps[apply(influx_cutoffs_rel, 1, 
            #                                        function(u){
            #                                            which(u >= influx_thresh)[1]
            #                                        })]
            # psr = mean(cutoffs_per_taxon)
            # psr = sd(cutoffs_per_taxon)
        }
        return(list("psr"=psr, "psa"=psa,
                    median_rel_influx = median_rel_influx, 
                    std_rel_influx = std_rel_influx,
                    "influx_threshold"=influx_thresh))
    }
    
    plotREVEALS_mod = function(d, p=NULL, p_meta = NULL, 
                               title, ytitle = "Age BP", reverse = T, 
                               error = "sd", taxa_to_plot="All"){
        
        lib("ggplot2")
        
        plot_pollen = F
        if ((is.null(p)==F) & (is.null(p_meta)==F)) { 
            plot_pollen=T 
        }
        
        if (!requireNamespace("forcats", quietly = T)) {
            stop("Package forcats needed for this function to work. Please install it.", 
                 call. = FALSE)
        }
        
        # -----------------------
        # TMP
        if ("Sample_ID" %in% colnames(d)){
            age = suppressWarnings(as.integer(d[,"Sample_ID"]))
            if (any(is.na(age))){ age  = 1:length(age)}
            name_age = "Sample_ID"
            if ("Age_BP" %in% colnames(d)){
                age_labels = d[,"Age_BP"]
            } else if ("Age_AD_BC" %in% colnames(d)){
                age_labels = d[,"Age_AD_BC"]
            } else {
                stop("Plot REVEALS: no age column found for y labels.")
            }
        } else if ("Age_BP" %in% colnames(d)){
            age = d[,"Age_BP"]
            name_age = "Age_BP"
        } else {
            age = d[,"Age_AD_BC"]
            name_age = "Age_AD_BC"
        }
        # -----------------------
        
        d.median = d[grep("median", names(d), value = T)]
        names(d.median) = sub(".median", "", names(d.median))
        d.mean = d[grep("mean", names(d), value = T)]
        names(d.mean) = sub(".mean", "", names(d.mean))
        d.10 = d[grep("q10", names(d), value = T)]
        names(d.10) = sub(".q10", "", names(d.10))
        d.90 = d[grep("q90", names(d), value = T)]
        names(d.90) = sub(".q90", "", names(d.90))
        d.sd = d[grep("sd", names(d), value = T)]
        names(d.sd) = sub(".sd", "", names(d.sd))
        
        # keep only taxa which have been passed 
        # (here, from pollen plot, taxa with pollen larger than 0 in more than x time slices)
        if (taxa_to_plot == "All"){
            taxa_to_plot = colnames(d.mean)
        } else if ((is.vector(taxa_to_plot)==T) & (is.character(taxa_to_plot) == T)){
            taxa_to_plot = taxa_to_plot[which(taxa_to_plot %in% colnames(d.mean)==T)]
        } else if ((is.numeric(taxa_to_plot) == T) & (is.scalar(taxa_to_plot) == T)){
            taxa_to_plot = colnames(d.mean)[which(apply(d.mean, 2, mean) >= taxa_to_plot)]
        } 
        # if only 1 taxon remains due to high min mean percentage, repeat with 0.5 %
        if (length(taxa_to_plot) <=1 ){
            taxa_to_plot = colnames(d.mean)[which(apply(d.mean, 2, mean) >= 0.05)]
        }
        
        d.median = d.median[,taxa_to_plot]
        d.mean = d.mean[,taxa_to_plot]
        d.10 = d.10[,taxa_to_plot]
        d.90 = d.90[,taxa_to_plot]
        d.sd = d.sd[,taxa_to_plot]
        
        yr = per = low = up = NULL
        if (plot_pollen==T){
            # -----------------------
            # TMP
            if ("Sample_ID" %in% colnames(p)){
                age_p = suppressWarnings(as.integer(p[,"Sample_ID"]))
                if (any(is.na(age_p))){ age_p  = 1:length(age)}
                name_age_p = "Sample_ID"
            } else if ("Age_BP" %in% colnames(p)){
                age_p = p[,"Age_BP"]
                name_age_p = "Age_BP"
            } else {
                age_p = p[,"Age_AD_BC"]
                name_age_p = "Age_AD_BC"
            }
            # -----------------------
            
            if (suppressWarnings(all(age_p == age, na.rm=T))){
                # resume both frames are in line
            } else if (all(age %in% age_p)){
                # drop samples which have been dropped in REVEALS analysis
                age_p = age_p[which(age_p %in% age == T)]
                p = p[which(as.integer(p[,name_age_p]) %in% age_p),]
            } else {
                stop("plotREVEALS_mod: Ages in pollen data and REVEALS results are different.")
            }
            
            taxa_cols = return_taxa(p, p_meta)
            #p_cols = c(name_age_p, taxa_cols)
            #p = p[,p_cols]
            p = p[,taxa_cols]
            # p = as.data.frame(cbind(age_p,
            #                   t(apply(p[,2:ncol(p)], 1, 
            #                      function(u){
            #                          u/sum(u)*100
            #                      }))))
            p = as.data.frame(t(apply(p, 1,
                                      function(u){
                                        u/sum(u)*100
                                      })))
            #p = cbind(age, p)
            #colnames(p) = p_cols
            #colnames(p) = taxa_cols
            p = p[,taxa_to_plot]
            
            df = data.frame(yr = rep(age, ncol(d.mean)), per = as.vector(as.matrix(d.mean)),
                            pollen = as.vector(as.matrix(p)),
                            taxa = as.factor(rep(colnames(d.mean), each = nrow(d.mean))))
            
        } else {
            df = data.frame(yr = rep(age, ncol(d.mean)), per = as.vector(as.matrix(d.mean)), 
                            taxa = as.factor(rep(colnames(d.mean), each = nrow(d.mean))))
        }
        
        # build frames to plot depending on sd or quantiles
        {
            df.mean = data.frame(yr = rep(age, ncol(d.mean)), per = as.vector(as.matrix(d.mean)), 
                                 taxa = as.factor(rep(colnames(d.mean), each = nrow(d.mean))))
            df.10 = data.frame(yr = rep(age, ncol(d.10)), per = as.vector(as.matrix(d.10)), 
                               taxa = as.factor(rep(colnames(d.10), each = nrow(d.10))))
            df.90 = data.frame(yr = rep(age, ncol(d.90)), per = as.vector(as.matrix(d.90)), 
                               taxa = as.factor(rep(colnames(d.90), each = nrow(d.90))))
            df.sd = data.frame(yr = rep(age, ncol(d.sd)), per = as.vector(as.matrix(d.sd)), 
                               taxa = as.factor(rep(colnames(d.sd), each = nrow(d.sd))))
            
            if (error != "sd") {
                #cat(paste0("\n\t\t\tPlot 10an 90% quantiles of reconstructed cover as error."))
                df$mean = df$per
                df$low = df.10$per
                df$up = df.90$per
            } else {
                #cat(paste0("\n\t\t\tPlot Std. deviation of reconstructed cover as error."))
                df$mean = df.mean$per
                df$low = df.mean$per - df.sd$per
                df$up = df.mean$per + df.sd$per
            }
        }
        
        theme_plot = theme(panel.grid.major.x = element_line(colour = "gray70",
                                                             size = 0.4,
                                                             linetype = "dotted"),
                           panel.grid.major.y = element_line(colour = "gray70",
                                                             size = 0.4,
                                                             linetype = "dotted"),
                           panel.grid.minor.x = element_line(colour = "gray70",
                                                             size = 0.2,
                                                             linetype = "dotted"),
                           panel.grid.minor.y = element_line(colour = "gray70",
                                                             size = 0.2,
                                                             linetype = "dotted"),
                           panel.background = element_rect(fill="white"),
                           plot.title = element_text(size = 18,
                                                     hjust = 0),
                           plot.subtitle = element_text(size = 16),
                           axis.text = element_text(size = 12.5),
                           axis.title = element_text(size = 14),
                           plot.margin = unit(rep(1.25, 4), "cm"))
        
       
        theme_new = theme(legend.position = "bottom",
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(), 
                          panel.grid.major.x = element_line(colour="grey85", size=0.075,
                                                            linetype="solid"), 
                          panel.grid.major.y = element_line(colour="grey85", size=0.075,
                                                            linetype="solid"),
                          panel.background = element_rect(color=NA,
                                                          fill="grey97"), 
                          axis.line = element_line(colour = "grey20"), 
                          strip.text.x = element_text(face = "bold", 
                                                      size = 12, 
                                                      angle = 0, 
                                                      vjust = 0.5, 
                                                      hjust = 0.5), 
                          strip.background = element_rect(color=NA,
                                                          fill="grey97"),
                          strip.text.y = element_text(angle = 0), 
                          panel.border = element_blank(), 
                          axis.text.x = element_text(angle = 0, 
                                                     hjust = 0.5))
        # plot
        {
            scale_max = round(max(df$per), -1)
            step_size = round(scale_max/5,-1)
            
            scale_max_x = max(age); scale_min_x = min(age)
            step_size_x = round(scale_max_x/20,0)
            if (is.na(step_size_x) | step_size_x < 1){ step_size_x = 1}
            inds_ages = seq(scale_min_x, scale_max_x, step_size_x)
            
            plt = ggplot2::ggplot(df) +
                # geom_ribbon(aes(yr, ymin = low, ymax = up), fill = "coral2", 
                #             colour = "coral2") + 
                geom_line(aes(yr, pollen), col = "darkred", size = 1) + 
                scale_y_continuous(breaks = seq(0, scale_max, step_size)) + 
                aes(ymin = step_size + 1) + 
                #labs(title = title) + 
                xlab(ytitle) + ylab("Relative Abundance (%)") + 
                coord_flip() + theme_plot + theme_new +
                facet_grid(~forcats::fct_inorder(df$taxa), scales = "fixed", 
                           space = "fixed")
            if (plot_pollen==T){
                plt = plt + geom_line(aes(yr, pollen), col = "darkgreen", size = 1)
            }
            if (name_age == "Sample_ID"){
                plt = suppressMessages(plt + 
                                       # scale_x_continuous(name=ytitle, 
                                       #                    breaks = inds_ages, 
                                       #                    labels = round(age_labels[inds_ages],-2)) +
                                       scale_x_reverse(name=ytitle,
                                                       breaks = inds_ages, 
                                                       labels = round(age_labels[inds_ages],-2)))
            }
        }
        
        return(plt)
    }
}

