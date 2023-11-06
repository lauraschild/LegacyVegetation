#function to reconstruct forests from compositional data
#input: dataframe with compositional information + Dataset_IDs
#output: average tree cover per site

reconstruct_forest <- function(REVEALS_df){
  trees <- class %>%
    filter(cat == "tree")%>%
    pull(taxa)
  
  REVEALS_df$forest <- rowSums(REVEALS_df[names(REVEALS_df) %in% trees],
                              na.rm = TRUE)
  
  arctic_forest <- REVEALS_df$forest - rowSums(as.data.frame(REVEALS_df[,names(REVEALS_df)%in% c("Betula","Betulaceae")]))
  REVEALS_df$forest[REVEALS_df$Latitude >= 70] <- arctic_forest[REVEALS_df$Latitude >= 70]
  
  forest_df <- REVEALS_df %>%
    group_by(Dataset_ID)%>%
    summarize(Longitude = mean(Longitude),
              Latitude = mean(Latitude),
              forest = mean(forest))
  
  return(forest_df)
  
}

reconstruct_forest_past <- function(composition){
  trees <- class %>%
    filter(cat == "tree") %>%
    pull(taxa)
  
  composition$forest <- rowSums(composition[names(composition) %in% trees],
                                na.rm = TRUE)
  arctic_forest <- composition$forest - rowSums(as.data.frame(composition[,names(composition)%in% c("Betula","Betulaceae")]))
  composition$forest[composition$Latitude >= 70] <- arctic_forest[composition$Latitude >= 70]
  
  return(composition)
}

# function to correct for unvegetated urban areas in reconstructed tree cover
# input: tree cover per site
# output: corrected tree cover per site

correct_urban <- function(forest_df){
  correct_forest_df <- merge(forest_df,
                             urban[,c("Dataset_ID","openness")],
                             all.x = TRUE) %>%
    mutate(forest = forest * (1- openness))%>%
    dplyr::select(-openness)
  
  return(correct_forest_df)
}


# function to validate reconstructed tree cover with landsat remote sensing tree cover
# input: reconstructed tree cover per site
# output: RSS and MAE

validate_forest <- function(correct_forest_df){
  valid <- merge(correct_forest_df,
                 RS)%>%
    rename(forest_RS = mean,
           forest_pollen = forest)%>%
    filter(!(is.na(forest_RS)))
  
  RSS <- sum((valid$forest_RS - valid$forest_pollen)^2, na.rm = TRUE)
  MAE <- mean(abs(valid$forest_RS - valid$forest_pollen), na.rm = TRUE)
  
  validation_result <- list(valid_df = valid,
                            RSS = RSS,
                            MAE = MAE)
  
  return(validation_result)
  
}

#function to plot validation results
plot_valid <- function(validation_result,
                       convergence,
                       type = "opti"){
  library(ggplot2)
  
  df <- validation_result$valid_df
  df$residual <- df$forest_pollen -df$forest_RS
  #plot scatter plot with RSS and MAe
  scatter <- ggplot(df,
                    aes(forest_RS,
                        forest_pollen))+
    xlim(c(0,100))+
    ylim(c(0,100))+
    geom_point(size= 1)+
    geom_abline(slope =1,
                col = ggsci::pal_futurama()(1))+
    theme_bw(base_size = 8)+
    annotate("label",
             x = rep(75,2),
             y = c(25,15),
             label = c(paste0("RSS = ",round(validation_result$RSS,2)),
                       paste0("MAE = ", round(validation_result$MAE,2))),
             alpha = 0.6,
             hjust  =0,
             size = 2,
             color = ggsci::pal_futurama()(1))+
    labs(x = "Remote sensing tree cover",
         y = "Reconstructed tree cover",
         subtitle = paste0("convergence = ",convergence)
         #axis.title = element_text(size = 2)
         )
  
  #plot map
  world <- rnaturalearth::ne_countries(scale = "small",
                                       returnclass = "sf")
  if(continent %in% c("Asia","Indopacific")){
    world <- sf::st_as_sf(maps::map("world2", plot = FALSE, fill = TRUE))
    df$Longitude[df$Longitude < 0] <- 360 + df$Longitude[df$Longitude < 0]
  }
  
  map <- ggplot(data =world)+
    geom_sf()+
    coord_sf(xlim = c(min(df$Longitude)-2,max(df$Longitude)+2),
             ylim = c(min(df$Latitude)-2, max(df$Latitude)+2))+
    geom_point(data = df,
               aes(Longitude,
                   Latitude,
                   col = residual),
               size = 1)+
    theme_minimal(base_size = 7)+
    labs(color = "Residuals")+
    scale_color_viridis_c(option = "A")+
    theme(legend.position = "right",
          #legend.title = element_text(size = 5),
          axis.title = element_blank(),
          #axis.text = element_text(size = 5),
          #legend.text = element_text(size = 4)
          )
  
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      paste0("Validation results for ",gsub("_"," ",continent)),
      #fontface = 'bold',
      x = 0,
      hjust = 0,
      size = 13
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 40)
    )
  
  
  #make cowplot::plot_grid
  plot_row <- cowplot::plot_grid(scatter, map,
                                 labels = "auto")
  total <- cowplot::plot_grid(title,plot_row,
                              ncol = 1,
                              rel_heights = c(0.1,1))+
    theme(plot.background = element_rect(fill = "white",color =NA))
  
  #save 
  file_name <- paste0("output/figures/valid_",continent,".png")
  if(type != "opti") file_name <- paste0("output/figures/valid_",type,"_",continent,".png")
  ggsave(plot = total,
         filename = file_name,
         width = 10,
         height = 3,
         dpi = 600)
  }

#function to summarize tree cover in 1000 yrs age slices
slice_forest <- function(forest_df,
                         slice_width = 1000,
                         last_slice = 20000){
  forest_df %>%
    as.data.frame()%>%
    filter(`Age_mean [yrs BP]` <= last_slice)%>%
    mutate(slice_end = cut(`Age_mean [yrs BP]`,
                       breaks = seq(-slice_width,last_slice,slice_width),
                       labels = seq(-slice_width, last_slice-slice_width,slice_width)
                       ))%>%
    group_by(Dataset_ID, Longitude, Latitude, slice_end)%>%
    summarize(forest = mean(forest))%>%
    filter(!is.na(slice_end))%>%
    return()
}

#function to summarize tree cover in 1000 yrs age slices
slice_taxon <- function(compo_df,
                        taxon = "Poaceae",
                        slice_width = 1000,
                        last_slice = 20000){
  
  age_col <-"Age_mean [yrs BP]"
  if("meanAgeBP" %in% names(compo_df)) age_col <- "meanAgeBP"
  slices <- compo_df %>%
    as.data.frame()%>%
    filter(get(age_col) <= last_slice)%>%
    mutate(slice_end = cut(get(age_col),
                           breaks = seq(-slice_width,last_slice,slice_width),
                           labels = seq(-slice_width, last_slice-slice_width,slice_width)
    ))%>%
    group_by(Dataset_ID, Longitude, Latitude, slice_end)%>%
    summarize(Taxon = mean(get(taxon)))%>%
    filter(!is.na(slice_end))
  
  names(slices)[names(slices) == "Taxon"] <- taxon
    
  return(slices)
}

