# This script harmonizes pollen nomenclature and adds the new chronologies to the records
# Annegreet Veeken

# Load libraries
library(neotoma) # for compile_downloads function
library(tidyverse)
library(readxl)
library(magrittr) # set colnames function

# load pollen data sets and conversion table
harm <- read_xlsx("Data/HarmonizationTablePollen_EPD_mottl.xlsx") %>% 
  select(VarName, Crop, ppe.taxon, GroupID, AccVarName) %>% 
  distinct
ds_neotoma <- readRDS("RDS_files/01-NeotomaSites-W-Europe.rds")
ds_davies <- readRDS("RDS_files/01_Pollen-Adavies.rds") %>% 
  purrr::map(., ~mutate(., site.name = as.character(site.name)))

# Harmonize nomenclature----
# compile neotoma datasets to one df per site
lPOL <- 
  ds_neotoma %>% 
  purrr::map(., ~compile_downloads(.)) 
# change "Vladař" to Vladar to avoid complications with matching
lPOL[[75]] <- lPOL[[75]] %>% 
  mutate(site.name = "Vladar")

# name list
sitenames <- 
  lPOL %>% 
  purrr::map(., ~pull(., "site.name") %>% unique()) %>% 
  unlist %>% as.character() 
names(lPOL) <- sitenames

# convert to long format
POLlong <- lPOL %>% 
  purrr::map(.,~pivot_longer(.,!.id:dataset, names_to = "taxa", values_to = "count") %>% 
               select(-age:-date.type) %>% 
               left_join(harm, by = c("taxa" = "VarName"))) %>% 
  bind_rows() %>% 
  select(-.id) %>% 
  mutate(dataset = as.numeric(dataset))

# add new chronologies from bchron ----
CAL <- readRDS("RDS_files/02_Calibrated_chronologies.rds") %>% 
  bind_rows() %>% 
  mutate(dataset = as.numeric(dataset))

POLcal <- 
  left_join(POLlong, CAL, by = c("site.name", "dataset", "depth"))

# Append pollen from Davies
DAVIES <- ds_davies %>% 
  purrr::map(.,~pivot_longer(.,!.id:dataset, names_to = "taxa", values_to = "count") %>% 
               left_join(harm, by = c("taxa" = "VarName"))) %>% 
  bind_rows() %>% 
  select(-.id) %>% 
  mutate(dataset = as.numeric(dataset))

lPOLcal <- bind_rows(POLcal, DAVIES) %>% 
  # replace NA in count by 0
  mutate(count = replace_na(count, 0)) %>% 
  # remove sites without proper chronology
  filter(!is.na(date.type)) %>% 
  filter(age < 10000) %>% 
  group_split(site.name) 
  
# Carquefou has 2 dated and 1 undated core
# Holzmaar has 1 dated and 1 undated core
# il fuorn has 1 dated and 1 core  of which chronological info is missing
# Lago di Origlio has 1 radiocarbon dated and 1 pb dated core 
# Lake of Annency has 1 dated and 1 undated core
# Krageholmssjön # only biostratigraphy
# Puy des Gouttes # only 1 radiocarbon date
# Lac de Pavin # only biostratigraphy

saveRDS(lPOLcal, "RDS_files/03_PollenWEU-Harmonised.rds")
 
