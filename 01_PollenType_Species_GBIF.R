# This script creates a conversion table for pollen type to species based on the current 
# distribution of terrestrial species belonging to the pollen taxon in NW-europe
# Annegreet Veeken, last changed 28-7-2021

# libraries
library(rgbif)
library(tidyverse)
library(Taxonstand)
library(CoordinateCleaner)
library(readxl)

# Taxa with avalable ppe's in https://doi.pangaea.de/10.1594/PANGAEA.922661
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx") 
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = ppe$`RPP v1 (Northern Hemisphere)`)

# create taxon list based on ppe data
taxa <- ppe$ppe.taxon
taxa <- taxa[!taxa %in% c("wild.herbs")]
taxa[taxa == "Sambucus.nigra"] <- "Sambucus nigra"
taxa[taxa == "Ericales"] <- "Ericaceae"
cerealia <- c("Secale cereale", "Avena", "Hordeum vulgare", "Panicum miliaceum", "Triticum",
              "Hordeum", "Zea")
taxa <- c(taxa, cerealia)

# divide in genera, families and subfamily, to facilitate GBIF search
sp <- taxa %>% str_subset("[:space:]") %>% str_replace("\\.", " ")
gn <- taxa %>% str_subset("^(?!.*eae)") 
gn <- gn[!str_detect(gn, "[:space:]")]
fam <- taxa %>% str_subset("ceae")

# find gbif id's
# for species
splist <-
  purrr::map(sp, ~ name_suggest(., rank = 'SPECIES', limit = 15))
spkeys <- splist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% sp) 

# for genera
gnlist <-
  purrr::map(gn, ~ name_suggest(., rank = 'GENUS', limit = 15))
gnkeys <- gnlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% gn) 

# for family
famlist <-
  purrr::map(fam, ~ name_suggest(., rank = 'FAMILY', limit = 15))
famkeys <- famlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% fam) 

# concatenate all taxon keys
taxonkeys <- c(spkeys$key, gnkeys$key, famkeys$key) 

# get country keys
countries <- c("Ireland", "United Kingdom",
               "Netherlands",
               "Germany", "France", "Denmark",
               "Switzerland", "Sweden",
               "Belgium", "Austria",
               "Czech Republic", "Luxembourg")
countrykeys <- isocodes %>% filter(name %in% countries) %>% pull(code)

## create function  that downloads, cleans, saves and creates conversion table for pollen type to species

gbifsaveclean <- function(taxonkey, countrykey){
  ## perform gbif search
  occsearch <- occ_search(taxonKey = taxonkey,
                        country = countrykey,
                        limit = 1000,  # evaluate limit setting
                        year = "1990, 2020",
                        hasCoordinate = T)     
  
# proceed to next steps of the function (cleaning) if the search retrieved data  
  if(!is.null(occsearch$data)) {
  # save raw data
  #saveRDS(occsearch, paste("RDS_files/RawGBIF/GBIFoccurances_", countrykey, taxonkey, ".rds", sep = ""))
  
  # only select data element, discard meta
  occdat <- occsearch$data 
  
  # clear empty elements (key didn't retrieve data)
  occdat <- occdat %>% compact()
  
  ## clean occurrences
  occdat <- occdat %>% #purrr::map(. 
    data.frame() %>%
      cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with empty coordinates
      cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with equal latitude and longitude coordinates
      cc_cap(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records within a certain radius around country capitals
      cc_cen(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records within a radius around the geographic centroids of political countries and provinces.
      cc_inst(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records assigned to the location of zoos, botanical gardens, herbaria, universities and museums
      cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with zero  long or lat
      cc_outl(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes geographical outliers
      cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude") # removes duplicates
     
  # save cleaned data frame
  #saveRDS(occdat, paste("RDS_files/CleanGBIF/CleanedGBIFOccurances_", countrykey, taxonkey, ".rds", sep = ""))

# create species per pollen taxon list
  spec <-
    occdat %>% 
    #purrr::map_dfr(. %>% 
    dplyr::select_if(colnames(.) %in% c("family", "genus", "species")) %>% 
    drop_na() %>% 
    arrange(species) 
  
  spec_num <-
    spec %>% 
    group_by(species) %>% 
    summarise(count = n()) %>% 
    drop_na() %>% 
    arrange(species)
  spec <- spec %>% distinct()
  spec <- cbind(spec, spec_num)
  spec <- spec[,-4] # remove duplicate species column
  
  # add column with pollen taxon the species are belonging to
  spec <-
    spec %>% mutate(pollentaxon = ifelse(genus %in% gn, genus, 
                                         ifelse(family %in% fam, family, "NA"))) 
  # standardise species names according to the Plant List
  spec.stand <-
    TPL(spec$species) 
  spec <-
    spec %>% mutate(stand.spec = str_c(spec.stand$New.Genus, 
                                       spec.stand$New.Species, sep = " "))
  
  saveRDS(spec, paste("RDS_files/PollenType_Species_tables/Species_PollenType_", countrykey,taxonkey, ".rds", sep = ""))
  }
}

# Apply function for all combinations of taxa and countries
combo <- crossing(taxonkeys, countrykeys) %>% 
  mutate(ID = paste(countrykeys, taxonkeys, sep = "")) 

purrr::map2(combo$taxonkeys, combo$countrykeys, 
           ~gbifsaveclean(taxonkey = .x, countrykey = .y))

# compile to dataframe
files <- 
  list.files("RDS_files/PollenType_Species_tables/") %>% 
  str_subset(pattern = "Species_PollenType_")

folderpath.fun <- function(x)
{paste("RDS_files/PollenType_Species_tables", x, sep = "/")}

specpol <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(., ~readRDS(.)) %>% 
  keep(., ~is.data.frame(.)) %>% 
  bind_rows() %>% 
  mutate(stand.spec = str_c(spec.stand$New.Genus, 
                                     spec.stand$New.Species, sep = " ")) %>% 
  dplyr::select(family, genus, species, stand.spec, pollentaxon)


specpol$pollentaxon[specpol$pollentaxon == "Ericaceae"] <- "Ericales"
specpol$pollentaxon[specpol$pollentaxon == "Hordeum vulgare"] <- "Hordeum"
specpol$pollentaxon[specpol$pollentaxon == "Zea"] <- "Cerealia"
specpol$pollentaxon[specpol$pollentaxon == "Secale cereale"] <- "Cerealia"
specpol$pollentaxon[specpol$pollentaxon == "Avena"] <- "Cerealia"
specpol$pollentaxon[specpol$pollentaxon == "Panicum miliaceum"] <- "Cerealia"
specpol$pollentaxon[specpol$pollentaxon == "Triticum"] <- "Cerealia"
specpol$pollentaxon[specpol$pollentaxon == "Thymelaeaceae"] <- "Thymelaceae"


# Read in TRY Ellenberg data
EVdata <- 
  vroom::vroom("Data/10082020_TRYdata_ellenberg/11356.txt", 
               col_select = c(ObservationID, AccSpeciesName, 
                              TraitName, DataName, OrigValueStr, StdValue))
EVdata <- EVdata  %>% 
  filter(TraitName == "Species environmental indicator value according to Ellenberg: moisture") %>% 
  select(AccSpeciesName, StdValue) 

# bind conversion table with ellenberg data
specpol <- specpol %>% 
  left_join(EVdata, by = c("stand.spec" = "AccSpeciesName" )) 

saveRDS(polspec, "RDS_files/01_PollenType_species.rds")
