#add missing Ks
#prepare K factor look up table
rm(list = ls())

source("scripts/REVEALS/reveals_and_PSA_functions.R", echo=FALSE)
library(parallel)
library(dplyr)

#get all possible fallspeeds
fallspeeds <- c()
for(continent in c("Europe")){
  source("scripts/supportive/load_pollen.R")
  source("scripts/supportive/prep_pollen_parameters.R")
  fallspeeds <- c(fallspeeds, parameters$Fallspeed)
  source("scripts/reconstruct/load_pollen_complete.R")
  source("scripts/supportive/get_new_parameters.R")
  fallspeeds <- c(fallspeeds, parameters$Fallspeed)
  rm(parameters)
}


vgs <- unique(fallspeeds)
K <- read.csv("input/K_factors.csv")

vgs <- vgs[!vgs %in% unique(K$fallspeed)]
#get unique combinations of basin types and basin diameters
source("scripts/supportive/edit_meta.R", echo=FALSE)
distinct_basins <- meta_df %>%
  filter(Dataset_ID %in% pollen_df$Dataset_ID) %>% 
  mutate(Basin_diameter = round(sqrt((4*Basin_Area)/pi)),
         Basin_Type = ifelse(Basin_Type == "other",
                             "peatland",
                             Basin_Type))%>%
  distinct(Basin_Type, Basin_diameter)



add_vg <- function(x){
  y <- distinct_basins
  y$fallspeed <- x
  return(y)
} 
distinct_variables <- lapply(vgs,
                             add_vg)%>%
  bind_rows()

get_K <- function(row){
  vg <- distinct_variables$fallspeed[row]
  tBasin <- distinct_variables$Basin_Type[row]
  dBasin <- distinct_variables$Basin_diameter[row]
  influx <- DispersalFactorK(vg = vg,
                             tBasin = tBasin,
                             dBasin = dBasin,
                             dwm = "gpm unstable",
                             regionCutoff = 1e+06,
                             u = 3)
  
  df <- data.frame(fallspeed = vg,
                   Basin_type = tBasin,
                   Basin_diameter = dBasin,
                   K = influx)
  return(df)
}

KFactors <- lapply(1:3,
                   get_K) %>% 
  bind_rows()
# cl <- makeCluster(n_cores)
# clusterExport(cl,
#               list("distinct_variables",
#                    "get_K",
#                    "DispersalFactorK",
#                    "gpm",
#                    "LakeModel"))
# 
# KFactors <- clusterApplyLB(cl,
#                            1:nrow(distinct_variables),
#                            get_K)%>%
#   bind_rows()
# 
# stopCluster(cl)

KFactors <- rbind(K,
                  KFactors)
write.csv(KFactors,
          "input/K_factors.csv",
          row.names = FALSE)
