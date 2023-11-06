#script to load and edit meta data to run in the optimization
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

#load meta data frame
meta_df <- read.csv("data/pollen_meta/metadata.txt",
                    sep = "\t")%>%
  dplyr::select(Dataset_ID, Basin_Type,Basin_Area)%>%
  #rename basin types to fit three categories
  mutate(Basin_Type = ifelse(Basin_Type %in% basin_types$Peatland,
                             "peatland",
                             ifelse(Basin_Type %in% basin_types$Lake,
                                    "lake",
                                    "other")),
         #fill NA basin area for peatlands with 100
         Basin_Area = ifelse(Basin_Type %in% c("peatland","other") & is.na(Basin_Area),
                             100,
                             Basin_Area))

#get global lake area median as fill value
global_median_lake_area <- median(meta_df$Basin_Area[meta_df$Basin_Type == "lake"],
                                  na.rm = TRUE)

#fill lake and other NA basin areas with global median
meta_df$Basin_Area[meta_df$Basin_Type %in% c("lake","other")&is.na(meta_df$Basin_Area)] <- global_median_lake_area

