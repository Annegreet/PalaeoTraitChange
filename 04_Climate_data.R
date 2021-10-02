# Extracting climate data from CHELSA-TraCE21k dataset
# This code is based on: https://github.com/JesperBorrePedersen/Hamburgian_Culture_Climate_Analysis
# Annegreet Veeken, last changed 28-7-2021
 
## Libraries
library(sf)
library(raster)
library(tidyverse)
library(cowplot)
library(rasterVis)
library(colorspace)

## Prepare climate data
# Raster options
rasterOptions(progress = "text")

# Load CHELSA files -- too big to upload
# download link: https://www.envidat.ch/#/metadata/chelsa_trace 
temp <- list.files(
  "file-path",
  pattern = "*.tif",
  full.names = T)

# convert to stack
temp <- stack(temp)

# Rename layers
names_temp <- gsub(".*bio01_(.*)_.*", "temp_\\1", names(temp))
names_temp <- gsub("temp_\\.(.*)", "temp_\\1_BCE", names_temp)
names_temp <- gsub("temp_([0-9]*)$", "temp_\\1_CE", names_temp)
names(temp) <- names_temp

# Create look up tables for years using BP 
temp_look_up <- data.frame(
  names_temp = names_temp,
  year_CE = as.numeric(gsub(".*_([0-9]*)_.*", "\\1", names_temp)) * 100,
  epoch = gsub(".*_([BCE]*)$", "\\1", names_temp),
  Time.BP = NA)

temp_look_up$Time.BP <- sapply(
  names_temp, 
  function(x){
    Time.BP <- temp_look_up[temp_look_up$names_temp == x,]$year_CE
    if (temp_look_up[temp_look_up$names_temp == x,]$epoch == "BCE"){
      Time.BP <- Time.BP + 2000
    } else {
      Time.BP <- 2000 - Time.BP 
    }
    return(Time.BP)
  }
)

temp_look_up <- arrange(temp_look_up, Time.BP)

## Load pollen site locations
lPOL <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds")

# Extract pollen location
lat <- lPOL %>% 
  purrr::map(., ~pull(.,lat) %>% unique) %>% unlist() 
long <- lPOL %>% 
  purrr::map(., ~pull(.,long) %>% unique) %>% unlist() 
sitenames <- names(lPOL)

pollenloc <- tibble(site.name = sitenames,
                    lat = lat,
                    long = long)

pollenloc <- st_as_sf(pollenloc, 
                      coords = c("long", "lat"),
                      crs = 4326,
                      remove = F)

## Extract climate variables for pollen locations
# Extract data from layers
temp_df <- as.data.frame(raster::extract(temp, 
                                         as_Spatial(pollenloc)))

temp_df$site.name <- pollenloc$site.name

temp_df <- pivot_longer(temp_df, 1:(ncol(temp_df)-1),
                        names_to = "names_temp", 
                        values_to = "temp") %>%  
  full_join(., temp_look_up) %>% 
  arrange(Time.BP)

# add time bin, average for 500 yr bins
timecat <- seq(from = 0, to = 11000, by = 500)
timecat.lab <- labels(timecat)
dftimecat <- tibble(time.bin = as.numeric(timecat.lab), Time.BP = timecat)
dftimecat <- map2_dfr(length(pollenloc$site.name), pollenloc$site.name, 
                      ~dftimecat %>% mutate(site.name = .y))

# climate data, binned in 500 years
temp_df_binned <- temp_df %>% 
  mutate(time.bin = as.numeric(cut(Time.BP, breaks = c(-Inf, timecat, Inf)), 
                               labels = timecat.lab)) %>% 
  group_by(site.name, time.bin) %>% 
  summarise(Temperature = mean(temp, na.rm = TRUE)) %>% 
  full_join(dftimecat, by = c("site.name", "time.bin"))

# save final files
saveRDS(temp_df_binned, file = "RDS_files/01_bio1_binned.rds") # 500 yr resolution
saveRDS(temp_df, file = "RDS_files/01_bio1_unbinned.rds")

