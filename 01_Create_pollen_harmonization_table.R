# This script merges tables from the EPD 2017 access distribution: http://www.europeanpollendatabase.net/data/downloads/ 
# to create a harmonization table
# Annegreet Veeken 

# Load libraries
library(tidyverse)

# Create 1 df with accepted pollen names and ecological groupings from epd data tables ----
epddirc <- "C:/Users/lgxgv/OneDrive - The University of Nottingham/Data analysis/PhD/Validation-study/Data/EPD data tables/"

# Accepted pollen types
p_vars <- read.csv(paste0(epddirc,"P_VARS.csv"), sep = ";")
# Ecological group names
p_group <- read.csv(paste0(epddirc, "P_GROUP.csv"), sep = ";")
# Group codes
group_code <- read.csv(paste0(epddirc, "GROUP_CODES.csv"), sep = ";")

# Filter for accepted pollen types names
acc_varname <- p_vars %>% 
  select(AccVar., VarName) %>% 
  distinct(AccVar., .keep_all = TRUE) %>%
  rename(AccVarName = VarName) %>% 
  # add original pollen types names
  left_join(p_vars, by = "AccVar.") %>% 
  select(AccVar. , AccVarName, Var., VarName) 

# Add ecological groups
harm_epd <- p_group %>% 
  left_join(acc_varname, by = "Var.") %>% 
  left_join(group_code, by = "GroupId") %>% 
  # select relevant columns and rename
  select(VarID = Var., GroupID = GroupId, AccVarID = AccVar.,
         AccVarName, VarName, GroupName) %>% 
  # replace punctuation and spaces by .
  mutate(VarName = str_replace(VarName, "[:punct:]", ".") %>% 
           str_replace("[:space:]", ".") %>% 
           str_replace("[:space:]", "."))

# Load harmonization table from Mottl et al. (2021) https://doi.org/10.1126/science.abg1685
# downloadable from https://doi.org/10.6084/m9.figshare.13049735.v2
harm_mottl <- read_xlsx("Data/Mottl_etal_EU_HarmonizationTable.xlsx",
                        sheet = 2) %>% 
  rename(VarName = NeotomaTaxonName, AccVarName = MHVar.2) %>% 
  # replace punctuation and spaces by .
  mutate(VarName = str_replace(VarName, "[:punct:]", ".") %>% 
           str_replace("[:space:]", ".") %>% 
           str_replace("[:space:]", "."))


# Check availability of taxa in harmonization table ----
ds_neotoma <- readRDS("RDS_files/01-NeotomaSites-W-Europe.rds")
ds_davies <- readRDS("RDS_files/01_Pollen-Adavies.rds") %>% 
  purrr::map(., ~mutate(.,site.name = as.character(site.name)))

# compile neotoma datasets to one df per site
lPOL <- 
  ds_neotoma %>% 
  purrr::map(., ~compile_downloads(.)) 

# Append pollen from Davies
lPOL <- append(lPOL, ds_davies)

neotomataxa <- lPOL %>% 
  purrr::map(., ~colnames(.)) %>% 
  unlist %>% 
  # calculate frequency of occurrence in pollen records
  table() %>% 
  as.data.frame() %>% 
  rename(VarName = 1, Freq = 2)

# Create one common harmonization table ----
harm <- neotomataxa %>% 
  left_join(harm_mottl, by = "VarName") %>% 
  left_join(harm_epd, by = "VarName", suffix = c("_mottl", "_epd")) %>% 
  # when harmonization of mottl is available, use this
  mutate(., AccVarName = case_when(!is.na(AccVarName_mottl) ~ AccVarName_mottl,
                                  # when harmonization of mottl is not available, use epd                                
                                  is.na(AccVarName_mottl) ~ AccVarName_epd,
                                  # when frequency of occurrence in the pollen records = 1 set to indeterminable
                                  Freq <= 1 ~ "Not applicable"))
# add ppe taxa 
pollen.equiv <- read_xlsx("Data/HarmonizationTablePollen.xlsx") %>% 
  select(VarName, Crop, species, genus, family, ppe.taxon)
harm <- harm %>% 
  left_join(pollen.equiv, by = "VarName") 

write.xlsx(harm, "Data/HarmonizationTablePollen_EPD_R.xlsx",
           showNA = FALSE, row.names = FALSE)
# Missing taxa in harmonisation list were added manually and saved in 
# HarmonizationTablePollen_mottl_EPD.xlsx