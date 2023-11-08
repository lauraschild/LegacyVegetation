#get new parameters
#load parameters for continent and edit
parameter_file <- "input/pollen_meta/parameters.xlsx"

sheet <- "CF_orig"
if(continent %in% c("South_America","Indopacific","Africa")) sheet <- "NH_orig"

parameters <- readxl::read_xlsx(parameter_file,
                                sheet = sheet)%>%
  #rename continents to fit conventions
  mutate(Continent = ifelse(sheet == "CF_orig" & Continent == "China",
                            "Asia",
                            ifelse(Continent == "North America",
                                   "North_America",
                                   Continent)),
         #change minimum and maximum fallspeed values
         Fallspeed = ifelse(Fallspeed < 0.01,
                            0.01,
                            ifelse(Fallspeed >0.15,
                                   0.15,
                                   Fallspeed)) )

#filter parameters if CF_orig
if(sheet == "CF_orig"){
  parameters <- parameters %>%
    filter(Continent == continent)
}

#expand parameters table with missing taxa
taxa <- names(pollen_df[,-(1:13)])
parameters <- parameters[match(taxa, parameters$Taxon),] 
mean_fallspeed <- mean(parameters$Fallspeed, na.rm = TRUE)
#fill NAs with fill values
parameters[is.na(parameters$PPE), 3:5] <- matrix(data = c(1,0.25,mean_fallspeed),
                                                 byrow = TRUE,
                                                 ncol = 3,
                                                 nrow = sum(is.na(parameters$PPE)))
parameters$Taxon <- taxa

#load optimized PPE values
optis <- read.csv(paste0("output/optimized_RPP_",continent,".csv"))

#replace old parameters
parameters[match(optis$Taxa, parameters$Taxon),'PPE'] <- optis$PPE
