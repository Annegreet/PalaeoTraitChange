# This script filters the pollen data, applies the ppe correction and calculates 
# pollen percentages of trees and shrubs, for use in the GAMs
# Annegreet Veeken

# Libraries
library(tidyverse)
library(readxl)

## Prepare pollen data: 
lPOL <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds")
sitenames <- lPOL %>% 
  purrr::map(., ~pull(., site.name)) %>% 
  unlist() %>% 
  unique()
names(lPOL) <- sitenames
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = ppe$`RPP v1 (Northern Hemisphere)`)

# add ppe taxon and group id 
lPOL <- lPOL %>% 
  purrr::map(. %>% 
               # filter out non-pollen
               filter(., !ppe.taxon %in%  "not applicable") %>% 
               # filter out ages past 10000
               filter(., age <= 10000))

# calculate percentage ppe corrected
lPOLcorrected <-  lPOL %>%
  purrr::map(. %>%
               filter(!ppe.taxon == "no ppe") %>% 
               # bind ppe data
               left_join(., ppe, by = "ppe.taxon") %>% 
               # adjust count with ppe
               mutate(adjustedcount = count*PPE) %>%
               group_by(site.name, age, depth) %>% 
               # calculate percentages
               mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
                      percent = count/sum(count, na.rm = TRUE)) %>%
               group_by(site.name, age, depth, GroupID, .drop = FALSE) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>% 
               # sort
               arrange(depth) %>% 
               filter(GroupID == "TRSH")) 

# check
lPOLcorrected[[1]]

saveRDS(lPOLcorrected, "RDS_files/03_TRSH_percentage.rds")
