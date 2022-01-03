# This script filters the pollen data, applies the ppe correction and calculates 
# pollen percentages of trees and shrubs, binned, for use in the PCA
# Annegreet Veeken, last changed 5-8-2021

# libraries
library(tidyverse)
library(readxl)
library(xlsx)

## Prepare pollen data: 
lPOL <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds")
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = ppe$`RPP v1 (Northern Hemisphere)`)
harm <- read_xlsx("Data/HarmonizationTablePollen.xlsx")
harm <- harm %>% 
  dplyr::select(AccVarName2, GroupId, ppe.taxon) %>% 
  distinct() %>% 
  filter(!is.na(ppe.taxon))

# convert dataframe to long format
lPOLlong <- lPOL %>% 
  purrr::map(.,
             ~dplyr::select_if(., !names(.) %in%
                                 c(".id", "depth", "age.old", "age.young",
                                   "date.type", "lat", "long", "dataset"))) %>%
  purrr::map(., ~gather(., key = "taxa", value = "count", -c(1:4)))

# add ppe taxon and group id 
lPOLlong <- lPOLlong %>% 
  purrr::map(. %>% 
               # WATCH OUT FOR DUPLICATION (check length of each df in list before and after join)
               left_join(., harm, by = c("taxa" = "AccVarName2")) %>% 
               # filter out non-pollen
               filter(., !ppe.taxon %in%  "not applicable"))


# discard samples with low sample count 
lowcount <- 
  lPOLlong %>% purrr::map(. %>% 
                            # filter out no ppe
                            filter(., !ppe.taxon %in%  "no ppe") %>%
                            group_by(time.bin) %>% 
                            summarise(pollencount = sum(count)) %>%
                            filter(pollencount < 300) %>% 
                            pull(time.bin)) 

lPOLlong <- 
  purrr::map2(lPOLlong, 1:length(lowcount), 
              ~filter(.x, !time.bin %in% lowcount[[.y]])) 

# at what sites is pollen sites too low?
sitelowcount <- lPOLlong %>% 
  keep(.,~nrow(.) == 0) %>% 
  names
sitelowcount

lPOLlong <- lPOLlong %>% 
  # filter sites without data (due to low count)
  discard(., ~nrow(.) == 0) %>% 
  # filter no ppe data
  purrr::map(., ~filter(., !ppe.taxon %in%  c("no ppe")) %>% 
               # select relevant columns
               dplyr::select(site.name, age, time.bin, 
                             Time.BP, taxa, GroupId, ppe.taxon,
                             count))



# calculate percentage per time bin, ppe corrected
lPOLcorrected <-  lPOLlong %>%
  purrr::map(. %>%
               # bind ppe data
               left_join(., ppe, by = "ppe.taxon") %>% 
               # adjusted count with ppe
               mutate(adjustedcount = count*PPE) %>%
               group_by(time.bin) %>% #
               mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
                      percent = count/sum(count, na.rm = TRUE),
                      time.bin = ) %>%
               group_by(time.bin, GroupId) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>%
               filter(GroupId == "TRSH"))

lPOLcorrected[[1]]

saveRDS(lPOLcorrected, "RDS_files/03_TRSH_percentage.rds")
