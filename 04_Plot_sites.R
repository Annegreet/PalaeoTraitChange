# This script plots a map of the study site location and time span of the record
# and the temperature map (supplementary information)
# Annegreet Veeken

library(tidyverse)
library(ggplot2)
library(ggmap)
library(leaflet)
library(htmltools)

# Load data
lPOL <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds")

# Create tibble with lat/long and time.span
lat <- lPOL %>% 
  purrr::map(., ~pull(.,lat) %>% unique) %>% unlist 
long <- lPOL %>% 
  purrr::map(., ~pull(.,long) %>% unique) %>% unlist 
sitenames <- names(lPOL)
age.old <- lPOL %>% 
  purrr::map(., ~pull(.,age) %>% max(.)) %>% unlist
age.young <- lPOL %>% 
  purrr::map(., ~pull(.,age) %>% min(.)) %>% unlist

sites <- tibble(site.name = sitenames,
                    lat = lat,
                    long = long,
                    age.old = round(age.old,0), 
                    age.young = round(age.young,0),
                    time.span = age.old - age.young)

# Europe map
eu <- ggplot2::map_data("world", 
                        region = c("UK", "Ireland","Netherlands",
                                   "Germany", "France", "Denmark",
                                   "Sweden", "Switzerland",
                                   "Czech republic",
                                   "Luxembourg", "Belgium"))
eumap <- ggplot() + 
  geom_polygon(data = eu, aes(x = long, y = lat, group = group),
                              colour = "#999999", fill = "#999999", 
               show.legend = FALSE) +
  coord_map("sinusoidal")

# map of study sites (figure 1)
pubmap <- eumap +
  geom_point(data = sites, aes(x = long, y = lat,
                               colour = time.span), size = 2) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  scale_color_steps(name = "Time span of pollen record",
                    high = "#132B43",
                    low = "#56B1F7",
                    breaks = c(2500,5000,7500),
                    labels = c("< 2500 years","","> 10 000 years"),
                    guide = guide_coloursteps(title.position = "top",
                                              barheight = 0.5, barwidth = 10,
                                              label.hjust = c(0.75,0,-0.1))) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
pubmap

ggsave("Figures/Fig1-MapStudySites.png", pubmap, dpi = 600,
       width = 170, height = 114, units = "mm")


# Temperature map (Supplementary information 6)
dfTEMP <- readRDS("RDS_files/01_bio1_binned.rds")

df <- left_join(dfTEMP, sites, by = c("site.name" = "site.name")) %>% 
  filter(Time.BP %in% c("1000", "2000", "3000", "4000",
         "5000", "6000", "7000", "8000", "9000"))

timelab <- c("1000 BP", "2000 BP", "3000 BP", "4000 BP", "5000 BP",
             "7000 BP", "6000 BP", "8000 BP", "9000 BP")

timebp <- df$Time.BP %>% unique %>% sort 
names(timelab) <- timebp

tempmap <- eumap + 
  geom_point(data = df, aes(x = long, y = lat, colour = Temperature), 
             size = 1.5) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  ggtitle("") +
  scale_color_viridis_c("Temperature (°C)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  facet_wrap(vars(Time.BP), ncol = 3, labeller = labeller(Time.BP = timelab))
tempmap

ggsave("Figures/SI6-TemperatureMap.png", tempmap, dpi = 600, height = 7, width = 7)

# leaflet map
lPOL_meta <- readRDS("C:/Users/lgxgv/OneDrive - The University of Nottingham/Data analysis/PhD/Meta-analysis/Analysis/RDS_files/01_Pollen_sites_meta.rds")
lPOL_meta <- lPOL_meta %>% 
  filter(site.name %in% sitenames) %>% 
  as_tibble()

popup <- paste0("<b>",lPOL_meta$site.name,"</b><br/>",
                "<b>Elevation</b> ", lPOL_meta$elev,"<br/>",
                "<b>Description</b> ", lPOL_meta$desc,"<br/>",
                "<b>Citation</b> ", lPOL_meta$citation) %>% 
  lapply(htmltools::HTML)
leaflet(lPOL_meta) %>% fitBounds(-7.68,71.27, 19.63, 42.22) %>% 
  addTiles() %>% 
  addMarkers(~long, ~lat, popup = popup)
