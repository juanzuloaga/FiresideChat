##############################################################

# MODEL EVALUTION SYSTEM - FIRESIDE CHAT - MATERIALS
  # By 

##############################################################

# 1. GENERAL SETTINGS -------------------

# 1.1 Installing ans loading packages ----------
  
Packages <- c(library(dplyr),
              library(reshape2),
              library(tidyr),
              library(kableExtra),
              library(sf),
              library(terra),
              library(tidyterra),
              library(tibble),
              library(ggplot2),
              library(leaflet),
              library(DT),
              library(gbm),
              library(wrapr),
              library(lubridate),
              library(openxlsx),
              library(shiny),
              library(widgetframe),
              library(leafem),
              library(tidyverse),
              library(htmltools),
              library(magrittr),
              library(quarto),
              library(leaflet.extras)
)

lapply(Packages, library, character.only =TRUE) 

# 1.2. Setting working directory --------------
setwd("C:/FiresideChat/")

# Setting geographic projection for interactive map in leaflet

crs_wgs84 <- "epsg:4326"


# 2. LOADING DATA --------------

# 2.1. List of species -----------
species_list <- read.csv("./Data/Focal Species _Fireside_Chat.csv")%>%
  dplyr::filter(To_run == "Y")

# 2.2. List of regions (BCRs) ----------------
regions <- c("can81")

# 2.3 Covariates table names ----------------
cov_label <- read.xlsx("./data/covariates_label_insert.xlsx")

# 2.4. Maps --------------

# Boundaries regions
bcr_boundary <- terra::vect("./data/BAM_BCR_NationalModel.shp")

# Boundary geographic for interactive leaflet
bcr_boundary_wgs84 <- terra::project(bcr_boundary, "epsg:4326")


# 2.3. Loading MODEL master file ------------------

load("./data/04_NM5.0_data_stratify_7181.R")

#Transform bird object (from matrix to df)
bird_d <- as.data.frame(as.matrix(bird))
names(bird_d)[1] <- "id"

# 3. LOOP (SPECIES + REGIONS)

for(i in 1:length(regions)){
  
# Selecting BCRs that are TRUE from the birdlist object
  
bcrlist_region <- dplyr::select(bcrlist, id, regions[i])

bcrlist_region <- bcrlist_region %>%
  dplyr::filter(bcrlist_region[2] == "TRUE")

for(j in 1:7){
  cat(paste0("Running....... SPECIES: ", species_list$Species_code[j], " (", species_list$Species_name[j], ") & REGION: ", regions[i]))  
  
  species_region <- paste0(species_list$Species_code[j], " (", species_list$Species_name[j], "): Region ",  regions[i])
  
# MASTER DATA FROM MODELS  

# Merging bcrlist_sel and bird_d by id, indicate if species was detected or not
  varsp <- species_list$Species_code[j]

  birds_merged_region <- merge(bcrlist_region, bird_d, by="id")%>%
    data.table::setnames(old=varsp, new="count") |> 
    subset(select = c("id", "count"))%>%

    dplyr::mutate(count=ifelse(count == 0, NA, as.integer(count)),
                  detection = ifelse(is.na(count), 0, 1))

 # adding the DATE, TIME, and COORDS
visit_detections_region <- merge(bcrlist_region, visit, by="id")%>%
  inner_join(birds_merged_region) |> 
  subset(select = c("id", "tssr", "jday", "count", "detection", "lat", "lon")) |> 
  pivot_longer(c(jday, tssr), names_to="metric", values_to="value") |> 
  mutate(metric = factor(metric, levels = c("tssr", "jday"),
                         labels = c("Time since sunrise", "Ordinal day")),
         detection = factor(detection, levels=c(0, 1), labels = c("Undetected", "Detected")))

# OBSERVATION IN REGIONS

locations_region <- st_as_sf(visit_detections_region, coords=c("lon", "lat"), crs=4326) |> 
  terra::vect() 

locations_region_detected_wgs84 <- locations_region[!is.na(locations_region$count), ]
locations_region_Undetected_wgs84 <- locations_region[is.na(locations_region$count), ]

# COVARIATES

# List of covariates

covlist_sel <- dplyr::filter(covlist, bcr %in% c(regions[i]))%>%
  dplyr::select(where(~ last(.x) != "FALSE"))%>%
  reshape2::melt(id.vars =  "bcr")
names(covlist_sel)[2] <- "Covariate.label"

merge_cov <- merge(covlist_sel, cov_label, by="Covariate.label")

# Filter list using BCR
cov_region <- filter(merge_cov, bcr == regions[i])%>%
  dplyr::select(Covariate.label, Category, Description)

colnames(cov_region)[1] <- "var"

# Covariates importance
if(!dir.exists(paste0("./data/",
                      species_list$Species_code[j],
                      "/covariates_bootstrps/",
                      regions[i], "/"))){
  next
}else{

list_bootstraps <- list.files(paste0("./data/",
                                     species_list$Species_code[j],
                                     "/covariates_bootstrps/",
                                     regions[i], "/"),
                              pattern=".R",
                              full.names = T)
}
df_bootstrap <- do.call(rbind,
                        lapply(list_bootstraps, function(x) {
                          
                          load(file = x)
                          gbm::summary.gbm(b.i, plotit = F)
                        }))

df_bootstrap_mean <- dplyr::summarise(df_bootstrap, .by=var, mean=mean(rel.inf))

df_bootstrap_sd <- dplyr::summarise(df_bootstrap, .by=var, sd=sd(rel.inf))

n = length(df_bootstrap$rel.inf)

df_bootstrap_sd$se <- df_bootstrap_sd$sd / sqrt(n)

df_bootstrap_merge <- merge(df_bootstrap_mean, df_bootstrap_sd, by="var")


# Merging list and importance
cov_region_merge <- merge(cov_region, df_bootstrap_merge, by="var")%>%
  dplyr::select(var, mean, Category, Description)%>%
  mutate_at(vars(mean), list(~ round(., 1)))
colnames(cov_region_merge)[1] <- "Covariate"
colnames(cov_region_merge)[2] <- "Relative importance"


                                
# INTERACTIVE MAP
print("Creating interactive map.........")
# Raster model Predictions(mean), Uncertainty(cv), Extrapolation
if(!file.exists(paste0("./data/", species_list$Species_code[j], "/",  species_list$Species_code[j], "_", regions[i], "_2020.tiff"))){
  next
}else{

raster_model_wgs84 <- terra::rast(paste0("./data/", species_list$Species_code[j], "/",  species_list$Species_code[j], "_", regions[i], "_2020.tiff"))%>%
                     terra::project(crs_wgs84)
}
# Palletes
pal <- colorFactor("viridis",locations_region_detected_wgs84$count,
                   na.color = "transparent")
pal2 <- colorNumeric("inferno", values(raster_model_wgs84$mean), 
                     na.color = "transparent", reverse = T)
pal3 <- colorQuantile("plasma", values(raster_model_wgs84$cv),
                      na.color = "transparent", reverse = T)
pal4 <- colorNumeric("Reds", values(raster_model_wgs84$extrapolation),
                     na.color = "transparent")

# Leaflet map
training_data_map <- leaflet() %>%
  
  # Setting view Canada
  setView(lng = -82, lat = 50, zoom = 4)%>%
  
  # Adding background map
  addProviderTiles(
    "Esri.WorldImagery",'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>%
  
  # Adding BCR map
  addPolygons(data=bcr_boundary_wgs84,
              weight=1,
              fillColor = '0',
              color = "white",
              label= paste0("BCR_", bcr_boundary_wgs84$subUnit),
              #  popup = bcr_boundary_wgs84$subUnit,
              group= "Regions")%>%
  # Adding species data DETECTED
  addCircleMarkers(data=locations_region_detected_wgs84,
                   fillColor = ~pal(locations_region_detected_wgs84$count),
                   weight = 0.9,
                   fillOpacity = 1,
                   radius = 3,
                   color = ~pal(locations_region_detected_wgs84$count),
                   label = ~htmlEscape(locations_region_detected_wgs84$count),
                   group="Observations_Detected")%>%
  
  
  
  # Adding legend species data detected
  addLegend('bottomright', pal = pal, values = locations_region_detected_wgs84$count,
            title = "Count", opacity = 1, group="Observations_Detected")%>%
  
  
  # Adding species data UNDETECTED
  addCircleMarkers(data=locations_region_Undetected_wgs84,
                   fillColor = "gray",
                   weight = 0.5,
                   radius = 3,
                   color = "gray",
                   group="Observations_Undetected")%>%
  
  
  
  # Adding raster PREDICTIONS
  addRasterImage(raster_model_wgs84[[i]]$mean,
                 colors=pal2,opacity=0.5,
                 group="Predictions",
                 layerId = "Predictions")%>%
 
  addLegend(pal = pal2,
            values = values(raster_model_wgs84[[i]]$mean),
            title = "Predictions (males/ha)",
            position = "bottomright",
            group = "Predictions")%>%
  
  
 
  
  # Adding raster UNCERTAINTY
  addRasterImage(raster_model_wgs84$cv,
                 colors=pal3, opacity=0.6, group="Uncertainty")%>%
  
  addLegend(pal = pal3, values = values(raster_model_wgs84$cv),
            title = "Uncertainty(CV-Quantiles)",
            position = "bottomright",
            layerId = "Uncertainty",
            group = "Uncertainty")%>%
  
  # Adding raster EXTRAPOLATION
  addRasterImage(raster_model_wgs84$extrapolation,
                 colors=pal4,opacity=0.6, group="Extrapolation")%>%
  
  addLegend(pal = pal4, values = values(raster_model_wgs84$extrapolation),
            title = "Extrapolation",
            position = "bottomright",
            layerId = "Extrapolation",
            group = "Extrapolation")%>%
  
  
  # add measurement tool
  addMeasure(primaryLengthUnit = "kilometers",
             secondaryLengthUnit = 'miles',
             primaryAreaUnit = "hectares",
             secondaryAreaUnit="acres",
             position = 'topleft')%>%
  
  # Adding layer control
  addLayersControl(options = layersControlOptions(collapsed = T),
                   position ="topleft",
                   overlayGroups = c( "Regions",
                                      "Observations_Detected",
                                      "Observations_Undetected",
                                      "Predictions",
                                      "Uncertainty",
                                      "Extrapolation"))%>%
  hideGroup(c("Predictions",
              "Uncertainty", "Extrapolation",
              "Observations_Undetected"))%>%

  addResetMapButton()


                                
# Creating HTML files

output_doc <- paste0(species_list$Species_code[j],  "_",  regions[i])

print("Rendering .html.........")
rmarkdown::render(input = "./Fireside_template.Rmd",
                  output_format = "html_document",
                  output_file =  output_doc,
                  output_yaml = "./_quarto.yml",
                  output_dir = "./docs/")     



}

 

  
  
}


