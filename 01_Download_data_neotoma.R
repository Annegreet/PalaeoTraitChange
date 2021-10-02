# This script downloads pollen records from the Neotoma database and creates
# a meta-data xlsx published in supplementary information 2
# Annegreet Veeken

# libraries
library(neotoma)
library(tidyverse)
library(readxl)
library(xlsx)

# load file that contains site names
neotomasites_full <- 
  read_xlsx("Data/SelectedPollenRecords.xlsx")

# create vector with names of sites with pollen data in neotoma
neotomasites <-
  neotomasites_full %>% 
  filter(pollen == "y") %>% # only pollen data
  filter(Country %in% c("Ireland", "Great-Britain", 
                     "France","Germany", "Sweden", 
                     "Czech", "Slovakia;Czech",
                     "Switzerland")) %>% 
  filter(`radiocarbon dating` == "y") %>% # only dated records
  filter(str_detect(source, "Neotoma")) %>% # only sites in neotoma
  pull(site.name) %>% sort() %>% strsplit(",") %>%  # create vector
  unlist() %>% unique()

# add wildcards for search
wildcard <- rep("*", length(neotomasites))
neotomasites_wild <-
  neotomasites %>% str_replace("bog", "") %>%  str_trim() %>% 
  str_c(wildcard, ., wildcard) 

# get site info
ds_neotoma_site <-
  neotomasites_wild %>% purrr::map( ~neotoma::get_site(.)) 

# get meta info of pollen data from sites
# This line could run in to this error: 
# Timeout was reached: [api.neotomadb.org] Connection timed out 
# (run line again till it works)
ds_neotoma_meta <-
  ds_neotoma_site %>% 
  purrr::map( ~neotoma::get_dataset(., datasettype = "pollen"))

# download the datasets
ds_neotoma <-
  ds_neotoma_meta %>%
  purrr::map( ~neotoma::get_download(.))

# save dataset
saveRDS(ds_neotoma, "RDS_files/01-NeotomaSites-W-Europe.rds")

# Create dataframe with meta info of pollen sites
# get publication info
publ <- ds_neotoma %>% 
  purrr::map(~neotoma::get_publication(.))
neotoma_publ <- 
  purrr::map_depth(publ, 3, pluck, "meta")  %>% 
  flatten() %>% 
  purrr::map(., ~bind_rows(.)) %>% 
  purrr::map(., ~dplyr::select(., citation) %>% 
               summarize(citation = paste(citation, collapse="/"))) %>% 
  purrr::map2(., neotoma_names, ~mutate(.x, site.name = .y)) %>% 
  bind_rows()

# get site info
neotoma_long <-
  purrr::map_depth(ds_neotoma_meta, 3, pluck, "long") %>% unlist()
neotoma_lat <-
  purrr::map_depth(ds_neotoma_meta, 3, pluck, "lat") %>% unlist()
neotoma_names <-
  purrr::map_depth(ds_neotoma_meta, 3, pluck, "site.name") %>% unlist()
neotoma_elev <-
  purrr::map_depth(ds_neotoma_meta, 3, pluck, "elev") %>% unlist()
neotoma_desc <-
  purrr::map_depth(ds_neotoma_meta, 3, pluck, "description") %>% unlist()
  
neotomameta <-
  cbind.data.frame(site.name = neotoma_names,
        long = neotoma_long,
        lat = neotoma_lat,
        elev = neotoma_elev,
        desc = neotoma_desc) %>% 
  bind_cols(neotoma_publ) %>% 
  distinct()

saveRDS(neotomameta, "RDS_files/01_Pollen_sites_meta.rds")
write.xlsx(neotomameta, "Output/01_PollenSites-meta.xlsx", row.names = FALSE)
