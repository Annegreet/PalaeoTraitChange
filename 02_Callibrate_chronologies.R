# This script creates new chronologies for the pollen records from neotoma
# Annegreet Veeken

# code source used as reference:
# https://cran.r-project.org/web/packages/Bchron/vignettes/Bchron.html#running-the-bchronology-age-depth-model
# https://fishandwhistle.net/post/2018/comparing-approaches-to-age-depth-modelling-in-r/ 
# http://open.neotomadb.org/workbooks/AgeModels.htm

# Libraries
library(readxl)
library(tidyverse)
library(Bchron)
library(neotoma)

# Load data ----
ds_neotoma <- readRDS("RDS_files/01-NeotomaSites-W-Europe.rds")

depth <- ds_neotoma %>% purrr::map_depth(., 2, ~.$sample.meta %>%  pull(depth))
sitenames <- ds_neotoma %>% purrr::map_depth(., 2, ~.$dataset$site.data$site.name) %>% 
  unlist %>%  unique

# get chronological info
ds_neotoma_chron <-
  ds_neotoma %>%
  purrr::map(., ~neotoma::get_chroncontrol(.))

# filter out undated records and records that are too short
datasetid <- ds_neotoma_chron %>%
  purrr::map(., ~names(.)) %>%
  unlist
datasetidname <- ds_neotoma_chron %>%
  purrr::map_depth(., 2, ~.$parent$dataset.name) %>%
  unlist

chron <- ds_neotoma_chron %>%
  purrr::map_depth(., 2, pluck, "meta") %>%
  purrr::map(., ~bind_rows(.)) %>%
  purrr::map2(., sitenames, ~mutate(.x, sitenames = .y)) %>%
  bind_rows() %>%
  mutate(., datasetid = datasetid)

control <- ds_neotoma_chron %>%
  purrr::map_depth(., 2, pluck, "chron.control") %>%
  flatten %>% 
  purrr::map(., ~pull(., control.type)) %>% 
  unlist %>% 
  unique()

# only include sites that have a chronology and cover more than 500 yr
selectid <- chron %>%
  filter(!is.na(age.type)) %>%
  filter(age.older - age.younger > 500) %>%
  pull(datasetid)

# created df suitable for Bchron
uncalib <- ds_neotoma_chron %>% 
  purrr::map_depth(., 2, ~.$chron.control) %>% 
  purrr::map_depth(., 2, ~filter(., !is.na(control.type)) %>% 
                     transmute(., ages = age, 
                               ageSds = ifelse(age < 71, 250, 
                                             ifelse(is.na(age.old), 250, age.old - age)), 
                               positions = depth, 
                               calCurves = ifelse(control.type %in% c("Core top", "Section top", "Core top, estimated",
                                                                      "Biostratigraphic, pollen", "Firbas pollen-zone boundary",
                                                                      "Biostratigraphic, pollen, regional", "Guess"), "normal",
                                                  ifelse(age < 71, "normal", "intcal13")), # radiocarbon not suitable for top sediments younger than 71 BP
                               ids = chron.control.id,
                               positionThicknesses = ifelse(is.na(thickness), 1, thickness)) %>% 
                     filter(., !ageSds < 0 )) %>% # remove if SD is negative as probably something wrong 
  # remove data sets without chronology or when the records are too short
  purrr::map(., ~keep(., names(.) %in% selectid)) 

ncores <- purrr::map(uncalib, ~length(.))
sitenamesn <- purrr::map2(sitenames, ncores, ~rep(.x, times = .y)) %>% unlist()

# flatten and create site.name column
uncalib <- uncalib %>% 
  flatten() %>% 
  purrr::map2(., sitenamesn, ~mutate(., site.name = .y)) %>%
  # discard empty elements
  discard(~length(.) == 0) 
  
# get relevant ids to select the depths for predictions
ids <- names(uncalib)

depth <- depth %>% 
  flatten %>% 
  keep(., names(.) %in% ids) 

# Run Bchron----
for(i in 1:length(uncalib)){
  calib <- with(uncalib[[i]], Bchronology(ages = ages, ageSds = ageSds, 
                                          calCurves = calCurves,
                                          positions = positions,
                                          ids = paste(site.name,ids),
                                          positionThicknesses = positionThicknesses,
                                          predictPositions = depth[[i]]))
  saveRDS(calib, paste0("RDS_files/Calibrations/02_Calibrated_chronology_",
                        uncalib[[i]]$site.name[[1]],"_",names(uncalib[i]),".rds"))
}

# compile to dataframe
files <- 
  list.files("RDS_files/Calibrations/") %>% 
  str_subset(pattern = "02_Calibrated_chronology_")

folderpath.fun <- function(x)
{paste("RDS_files/Calibrations/", x, sep = "/")}

chronid <- files %>% 
  str_remove(pattern = "02_Calibrated_chronology_") %>% 
  str_remove(pattern = ".rds") %>% 
  str_split("_", simplify = TRUE) %>% 
  as.data.frame(col.names = c("site.name", "dataset")) %>% 
  split(., seq(nrow(.)))

l <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(., ~readRDS(.))

calage <- purrr::map2(l, chronid,  ~tibble(age = apply(.x$thetaPredict, 2, "quantile", 
                                                       probs = 0.5, na.rm = TRUE), 
                                           age.old = apply(.x$thetaPredict, 2, "quantile", 
                                                           probs = 1 - (1 - 0.95)/2, na.rm = TRUE), 
                                           age.young = apply(.x$thetaPredict, 2, "quantile", 
                                                             probs = (1 - 0.95)/2, na.rm = TRUE),
                                           depth =  .x$predictPositions,
                                           site.name = .y[[1]],
                                           dataset = .y[[2]],
                                           date.type = "Bchron")
) 

saveRDS(calage, "RDS_files/02_Calibrated_chronologies.rds")

# Testing age uncertainties----
# Run Bchron
for(i in length(uncalib)){
  calib <- with(uncalib[[i]][1:85,], Bchronology(ages = ages, ageSds = ageSds, 
                                          calCurves = calCurves,
                                          positions = positions,
                                          ids = paste(site.name,ids),
                                          positionThicknesses = positionThicknesses,
                                          predictPositions = depth[[i]]))
  age_uncertainties <- 
    Bchron:::predict.BchronologyRun(
      calib,
      newPositions = depth[[i]])
  
  saveRDS(age_uncertainties, paste0("RDS_files/02_Calibrated_chronology_age_uncertainties_",
                        uncalib[[i]]$site.name[[1]],"_", names(uncalib[i]),".rds"))
}

# compile to one list
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "02_Calibrated_chronology_age_uncertainties_")

corename <- files %>% 
  str_remove(., "02_Calibrated_chronology_age_uncertainties_") %>% 
  str_remove(., ".rds")

corename <- files %>% 
  str_remove(., "02_Calibrated_chronology_age_uncertainties_") %>% 
  str_remove(., ".rds")
  
folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

age_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(., corename, ~readRDS(.) %>%
                t %>% 
                as.data.frame %>% 
                # add depth column (rownames of the dataframe)
                rownames_to_column(var = "position") %>% 
                mutate(depth = str_remove(position, "Pos") %>% as.numeric(),
                       # add site.name and dataset id
                  corename = .y, .before = V1) %>% 
                separate(corename, into = c("site.name", "dataset.id"), 
                         sep = "_", convert = TRUE)) 

saveRDS(age_list, "RDS_files/02_Calibrated_chronologies_uncertainties.rds")