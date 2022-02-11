# This script filters the pollen data, applies the ppe correction and calculates 
# pollen percentages
# Annegreet Veeken

# Libraries
library(tidyverse)
library(readxl)

# Load data
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

lPOL %>% 
  purrr::map(. %>% 
               group_by(age) %>% 
               # calculate percentages
               mutate(percent = count/sum(count, na.rm = TRUE)) %>%
               group_by(site.name, age, depth, ppe.taxon) %>%
               summarise(percent = sum(percent)) %>% 
               # sort
               arrange(age,depth, ppe.taxon) %>% 
               filter(ppe.taxon == "no ppe")) %>% 
  bind_rows() %>% 
  summary()
  

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
               group_by(site.name, age, depth, ppe.taxon) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>% 
               # sort
               arrange(depth, ppe.taxon))

# check if adjusted.percent > 0, otherwise PPE is probably missing
# x <- lPOLcorrected %>% purrr::map(., ~group_by(., age) %>%
#                                     summarise(., sum = sum(.$adjustedpercent)) %>%
#                                     filter(sum == 0))

saveRDS(lPOLcorrected, "RDS_files/04_PollenWEU-PPEcorrected.rds")


lPOLcorrected_herbs <-  lPOL %>%
  purrr::map(. %>%
               filter(!ppe.taxon == "no ppe") %>% 
               # bind ppe data
               left_join(., ppe, by = "ppe.taxon") %>% 
               # filter for herbs 
               filter(GroupID == "HERB") %>% 
               # adjust count with ppe
               mutate(adjustedcount = count*PPE) %>%
               group_by(site.name, age, depth) %>% 
               # calculate percentages
               mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
                      percent = count/sum(count, na.rm = TRUE)) %>%
               group_by(site.name, age, depth, ppe.taxon) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>% 
               # sort
               arrange(depth, ppe.taxon))

saveRDS(lPOLcorrected_herbs, "RDS_files/04_PollenWEU-PPEcorrected-HERBS.rds")

lPOLcorrected_trees <-  lPOL %>%
  purrr::map(. %>%
               filter(!ppe.taxon == "no ppe") %>% 
               # bind ppe data
               left_join(., ppe, by = "ppe.taxon") %>% 
               # filter for trees 
               filter(GroupID %in% c("TRSH","DWAR")) %>% 
               # adjust count with ppe
               mutate(adjustedcount = count*PPE) %>%
               group_by(site.name, age, depth) %>% 
               # calculate percentages
               mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
                      percent = count/sum(count, na.rm = TRUE)) %>%
               group_by(site.name, age, depth, ppe.taxon) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>% 
               # sort
               arrange(depth, ppe.taxon))

saveRDS(lPOLcorrected_trees, "RDS_files/04_PollenWEU-PPEcorrected-TRSH.rds")
