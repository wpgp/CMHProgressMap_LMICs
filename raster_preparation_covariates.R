library(tidyverse)
library(terra)
library(tidyterra)

##File to prepare rasters and create covariates at cluster locations from rasters


##Creation of a cropped mastergrid to remove areas where inland water is 100%

mastergrid <- rast("raster/mastergrid_1km.tif")

inland_water_pct <-
  rast("raster/inland_water_pct_100m.tif") %>%
  resample(mastergrid, method = "average")


mastergrid_water_subtract <- mastergrid

values(mastergrid_water_subtract)[values(inland_water_pct) == 100] <-
  NA


district_boundaries <-
  vect("shapes/round2/sdr_subnational_boundaries.shp") %>%
  project(mastergrid) %>%
  arrange(REGCODE) %>%
  mutate(id = 1:nrow(.))



cluster_locations1 <-
  vect("shapes/round1/cluster_locations.shp") %>%
  filter(SOURCE != "MIS")

cluster_locations2 <- vect("shapes/round2/cluster_locations.shp") %>%
  filter(SOURCE != "MIS")

all_cluster_locations <-
  rbind(cluster_locations1, cluster_locations2)

##Remove polygons (islands) with no cluster locations
district_boundaries <-
  disagg(district_boundaries)[all_cluster_locations]

mastergrid_water_subtract_district_crop <-
  mastergrid_water_subtract %>% crop(district_boundaries, mask = T)

writeRaster(mastergrid_water_subtract_district_crop,"raster/mastergrid_crop.tif",overwrite=T)



## function to create covariates at cluster level from raster data


##Function to remove values where inland water percentage is 100% and then
##resample to the 1km mastergrid

cut_and_resample <- function(raster) {
  inland_water_pct <- rast("raster/inland_100m.tif")
  mastergrid_water_subtract_district_crop <- rast("raster/mastergrid_crop.tif")
  
  
  ## If unprepared raster is at >100m resample to 100m
  if (res(raster)[1] > 0.0009) {
    raster <- resample(raster, inland_water_pct)
  }
  
  
  raster_1km <- terra::resample(raster, mastergrid_water_subtract_district_crop, method =
                                  "average")
  
  values(raster_1km)[is.na(values(mastergrid_water_subtract_district_crop))] <-
    NA
  
  return(raster_1km)
  
  
  
  
}



##Function to extract covariates from the separate urban and rural buffers


buffer_extraction <- function(urban_buffer, rural_buffer, raster) {
  ##Extract covariates from raster using 2km urban buffers
  urban_covariates <-
    terra::extract(
      raster,
      urban_buffer,
      fun = "mean",
      exact = T,
      na.rm = T,
      bind = T
    )
  
  ##Extract covariates from raster using 5km rural buffers
  rural_covariates <-
    terra::extract(
      raster,
      rural_buffer,
      fun = "mean",
      exact = T,
      na.rm = T,
      bind = T
    )
  
  
  ##Join urban and rural extractions
  
  covariates <-
    rbind(urban_covariates, rural_covariates) %>%
    arrange(DHSCLUST) %>%
    values()
  
  
  
  return(covariates)
  
}


raster_covariates_preparation <-
  function(raster_file_list,
           cluster_location_file,
           covariates_filename) {
    ##Read in cluster locations from file and remove those with missing coords
    cluster_locations <- vect(cluster_location_file) %>%
      filter(SOURCE != "MIS")
    
    ##Create 2km urban buffer around urban cluster locations
    urban_buffer <-
      cluster_locations %>% filter(URBAN_RURA == "U") %>%
      buffer(width = 2000)
    
    ##Create 5km rural buffer around rural cluster locations
    rural_buffer <-
      cluster_locations %>% filter(URBAN_RURA == "R") %>%
      buffer(width = 5000)
    
    ##Initialise covariate object
    
    full_covariates <- cluster_locations
    
    ##for each raster extract covariates and then cut and resample the raster
    
    for (raster_file in raster_file_list) {
      print(raster_file)
      
      raster <- rast(raster_file)
      
      
      covariates <- buffer_extraction(urban_buffer, rural_buffer, raster)
      colnames(covariates)[21] <- basename(raster_file) %>% str_remove(".tif")
      
      full_covariates <-
        full_covariates %>%
        left_join(covariates)
      
      
      raster_1km <- cut_and_resample(raster)
      raster_1km_file <- raster_file %>% str_remove("/unprepared")
      print(raster_1km_file)
      
      writeRaster(raster_1km, raster_1km_file, overwrite = T)
      
    }
    
    write_csv(values(full_covariates), covariates_filename)
    
  }



raster_file_list_round1 <- list.files("raster/unprepared/round1",
                                      full.names = T)

cluster_location_file_round1 <- "shapes/round1/cluster_locations.shp"

covariates_filename_round1 <- "covariates/covariates_round1.csv"

raster_covariates_preparation(
  raster_file_list_round1,
  cluster_location_file_round1,
  covariates_filename_round1
)


raster_file_list_round2 <- list.files("raster/unprepared/round2",
                                      full.names =T)

cluster_location_file_round2 <- "shapes/round2/cluster_locations.shp"

covariates_filename_round2 <- "covariates/covariates_round2.csv"

raster_covariates_preparation(
  raster_file_list_round2,
  cluster_location_file_round2,
  covariates_filename_round2
)


