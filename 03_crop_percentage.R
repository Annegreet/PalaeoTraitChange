# This script filters the pollen data, applies the ppe correction and calculates 
# pollen percentages of crops, unbinned, for use in the GAMS
# Annegreet Veeken

# Libraries
library(tidyverse)
library(readxl)
library(xlsx)

lPOL <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds")
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = ppe$`RPP v1 (Northern Hemisphere)`)
taxa <- lPOL %>% purrr::map(., ~colnames(.)) %>% unlist() %>% unique()
harm <- read_xlsx("Data/HarmonizationTablePollen.xlsx")
harm <- harm %>% 
  filter(AccVarName2 %in% taxa) %>% 
  dplyr::select(AccVarName2, GroupId, Crop, ppe.taxon) %>% 
  distinct() 
harm$Crop[is.na(harm$Crop)] <- "Non-crop"
harm$Crop <- as.factor(harm$Crop)

# convert data frames to long format
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
               filter(., !ppe.taxon %in%  "not applicable") %>% 
               # filter out ages past 10000
               filter(., age <= 10000)) 

lPOLlong <- lPOLlong %>% 
  # filter no ppe data
  purrr::map(., ~ filter(., !ppe.taxon %in%  c("no ppe")) %>% 
               # select relevant columns
               dplyr::select(site.name, age, time.bin, 
                             Time.BP, taxa, GroupId, ppe.taxon, Crop,
                             count))

# calculate percentage per time bin, ppe corrected
lPOLcorrected <-  lPOLlong %>%
  purrr::map(. %>%
               # bind ppe data
               left_join(., ppe, by = "ppe.taxon") %>% 
               # adjusted count with ppe
               mutate(adjustedcount = count*PPE) %>%
               group_by(age) %>% #
               mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
                      percent = count/sum(count, na.rm = TRUE),
                      time.bin = as.factor(time.bin)) %>%
               ungroup() %>% 
               group_by(age, Crop, .drop = FALSE) %>%
               summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
                         percent = sum(percent)) %>% 
               filter(Crop == "Crop")
  )

#check
lPOLcorrected[[1]]
lPOLcorrected[["Sluggan Moss"]]

saveRDS(lPOLcorrected, "RDS_files/03_Crop_percentage_unbinned.rds")
