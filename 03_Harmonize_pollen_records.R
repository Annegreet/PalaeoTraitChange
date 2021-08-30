# This script harmonizes pollen nomenclature and adds the new chronologies to the records
# Annegreet Veeken, last changed 28-7-2021

library(neotoma)
library(tidyverse)
library(readxl)
library(magrittr) # set colnames function

# load pollen data sets and conversion table
pollen.equiv <- read_xlsx("Data/HarmonizationTablePollen.xlsx")
ds_neotoma <- readRDS("RDS_files/01-NeotomaSites-W-Europe.rds")
ds_davies <- readRDS("RDS_files/01_Pollen-Adavies.rds") %>% 
  purrr::map(., ~mutate(.,site.name = as.character(site.name)))

lPOL <- 
  ds_neotoma %>% 
  purrr::map(., ~compile_downloads(.)) 
# Append pollen from Davies
lPOL <- append(lPOL, ds_davies)

# name list
sitenames <- 
  lPOL %>% 
  purrr::map(., ~pull(., "site.name") %>% as.character() %>% unique()) %>% 
  unlist %>% as.character() 
names(lPOL) <- sitenames

#create categorical variable for time period
timecat <- seq(from = 0, to = 10500, by = 500)
timecat.lab <- labels(timecat)
dftimecat <- bind_cols(Time.BP = timecat, time.bin = timecat.lab)
timecat.lab[23] <- "23"

lPOL <- lPOL %>% 
  purrr::map(., 
             ~mutate(., time.bin = cut(as.numeric(age), 
                                       breaks = c(-Inf, timecat, Inf), 
                                       labels = timecat.lab)) %>% 
               left_join(dftimecat, by = "time.bin"))
# add calibrated chronologies to the dataframes

# split pollen data from rest of the variables
lPOlvar <- lPOL %>% 
  purrr::map(., ~select_if(., names(.) %in% 
                            c(".id", "site.name", "depth", 
                              "age", "age.old", "age.young", 
                              "date.type", "lat", "long", 
                              "dataset", "time.bin", "Time.BP")))
lPOL <- lPOL %>% 
  purrr::map(., ~select_if(., ! names(.) %in% 
                             c(".id", "site.name", "depth", 
                               "age", "age.old", "age.young", 
                               "date.type", "lat", "long", 
                               "dataset", "time.bin", "Time.BP")))
# harmonize taxa
harmonize.taxa <- function(x){
taxon.matches <- match(colnames(x), 
                       pollen.equiv$VarName)
used.taxa <- pollen.equiv[taxon.matches, ]
x <- x %>% set_colnames(used.taxa$AccVarName2)
x <- t(rowsum(t(x), group = colnames(x), na.rm = T))
}

# add missing taxa to HamonizePOllenTaxa file
# tax <- lPOL %>% purrr::map(., ~colnames(.)) %>% unlist %>% unique(.)
# x <- tax[!tax %in% pollen.equiv$VarName]
# write.csv(x, "temp.csv")
# 

lPOL <- lPOL %>% 
  purrr::map(., ~harmonize.taxa(.))

lPOLharm <- purrr::map2(lPOlvar, lPOL, ~cbind.data.frame(.x, .y)) 

# Check for missing column names
# natax <- lPOLfinal %>% purrr::map(., ~colnames(.) %>% anyNA) 

# add new chronologies from bchron
lCAL <- readRDS("RDS_files/02_Calibrated_chronologies.rds") %>% 
  bind_rows() %>% 
  group_split(site.name)

calnames <- 
  lCAL %>% 
  purrr::map(., ~pull(., "site.name") %>% as.character() %>% unique()) %>% unlist
names(lCAL) <- calnames

addchron <- function(site.name){
  if(site.name %in% calnames) {
lPOLharm[[site.name]] <- lPOLharm[[site.name]] %>% select(-age:-date.type) %>% 
  left_join(lCAL[[site.name]], by = c("site.name","depth",".id" = "dataset")) %>% 
  select(.id, site.name, depth, age, age.old, age.young, date.type, everything())
  } 
  else{
    lPOLharm[[site.name]]
  }
}

lPOLcal <- purrr::map(sitenames, ~addchron(.x))

# check for missing chronologies
lPOLcal %>% purrr::map(., ~pull(., date.type) %>% unique) %>% unlist
lPOLcal[[75]] <- lPOLcal[[75]] %>% mutate(site.name = "Vladar") %>% 
    select(-age:-date.type) %>% 
    left_join(lCAL[["Vladar"]], by = c("site.name","depth",".id" = "dataset")) %>% 
    select(.id, site.name, depth, age, age.old, age.young, date.type, everything()) # Vladar didnt join because of spelling
lPOLcal[[7]] <- lPOLcal[[7]] %>% 
  filter(!is.na(age)) # Carquefou has 2 dated and 1 undated core
lPOLcal[[23]] <- lPOLcal[[23]] %>% 
  filter(!is.na(age)) # Holzmaar has 1 dated and 1 undated core
lPOLcal[[26]] <- lPOLcal[[26]] %>% 
  filter(!is.na(age)) # il fuorn has 1 dated and 1 core  of which chronological info is missing
lPOLcal[[37]] <- lPOLcal[[37]] %>% 
  filter(!is.na(age)) # Lago di Origlio has 1 radiocarbon dated and 1 pb dated core 
lPOLcal[[38]] <- lPOLcal[[38]] %>% 
  filter(!is.na(age)) # Lake of Annency has 1 dated and 1 undated core

# name list
polcalnames <- lPOLcal %>% 
  purrr::map(., ~pull(., site.name) %>% unique) %>% unlist
names(lPOLcal) <- polcalnames

# delete sites without proper chronology
lPOLcal[["Krageholmssjön"]] <- NULL # only biostratigraphy
lPOLcal[["Puy des Gouttes"]] <- NULL # only 1 radiocarbon date
lPOLcal[["Lac de Pavin"]] <- NULL # only biostratigraphy

saveRDS(lPOLcal, "RDS_files/03_PollenWEU-Harmonised.rds")

