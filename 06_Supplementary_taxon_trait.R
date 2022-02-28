# This script  produces the supplementary tables S3 and S7
# Annegreet Veeken

# Load packages
library(tidyverse)
library(readxl)
library(xlsx)

## S3 - Pollen taxa information -----
# Pollen data
lPOL <- readRDS("RDS_files/04_PollenWEU-PPEcorrected.rds")
pollentaxa <- lPOL %>% purrr::map(., ~.$ppe.taxon) %>% unlist %>% unique()

# Species conversion table
spec <- readRDS("RDS_files/01_PollenType_species.rds") %>% 
  as_tibble()

# trait data
gapdata <- read.csv("C:/Users/lgxgv/OneDrive - The University of Nottingham/Data analysis/PhD/Meta-analysis/Analysis/Data/GapfilledTraitDataFS.csv") 
gapdata <- gapdata %>% 
  mutate(Species = str_replace(Spp, "\\.", " ")) %>% 
  filter(Species %in% spec$stand.spec) # only the relevant species
trait <- gapdata %>% 
  inner_join(spec, by = c("Species" = "stand.spec")) %>% 
  arrange(pollentaxon) %>%  # sort alphabetically
  dplyr::select(pollentaxon,PlantHeight,SLA, LA = Area, LDMC, LeafC = CperDryMass, LeafP = P, 
                LeafN = N, Seed.length = SeedLength, Seed.count = SeedNr_reproductonUnit, 
                Seed.mass = SeedMass)

# PPE data
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = round(ppe$`RPP v1 (Northern Hemisphere)`,2))

# table
nspec <-  spec %>% filter(pollentaxon %in% pollentaxa) %>% group_by(pollentaxon) %>% summarise(nspec = n())
ntrait <- trait %>% filter(pollentaxon %in% pollentaxa) %>% group_by(pollentaxon) %>% summarise(ntrait = n())
ppe <- ppe %>% filter(ppe.taxon %in% pollentaxa)
polperc <- lPOL %>% bind_rows() %>% 
  ungroup() %>% 
  select(site.name, ppe.taxon) %>% 
  distinct() %>% 
  group_by(ppe.taxon) %>% 
  summarise(perc = round(n()/78*100,2))

s3 <- left_join(ppe, polperc, by = "ppe.taxon") %>% 
  left_join(ntrait, by = c("ppe.taxon" = "pollentaxon")) %>% 
  left_join(nspec, by = c("ppe.taxon" = "pollentaxon")) %>% 
  rename("Pollen taxon" = ppe.taxon, "Occurrence in pollen records" = perc,  
         "Number of trait observations" = ntrait, "Number of species in pollen taxon" = nspec) 

# write.xlsx(s3, file = "Output/06_Summary_table_pollentypes.xlsx")

## S7 Trait estimates----
pol <- tibble(ID = 1:47, pollentaxon = pollentaxa)
dfMEAN <- read_rds("RDS_files/05_multivariate_taxon_mean_100.rds") %>% 
  as.data.frame() %>% 
  rownames_to_column("parameter") %>%
  filter(str_detect(parameter, "mu")) %>%
  dplyr::select(parameter, Mean) %>% 
  mutate(ID = 
           str_extract(parameter, pattern = "\\[(.*?)\\,") %>% 
           str_remove_all(pattern = "[[:punct:]]") %>% 
           as.numeric(),
         trait = str_extract(parameter, pattern = "\\,(.*?)\\]") %>% 
           str_remove_all(pattern = "[[:punct:]]") %>% 
           as.factor() %>% 
           recode("1" = "SLA", 
                  "2" = "PlantHeight",
                  "3" = "LA", 
                  "4" = "LDMC",
                  "5" = "LeafC",
                  "6" = "LeafP",
                  "7" = "LeafN",
                  "8" = "Seed.length",
                  "9" = "Seed.count",
                  "10" = "Seed.mass")) %>% 
  left_join(pol, by = "ID") %>% 
  mutate(Mean = exp(Mean))

dfSD <- read_rds("RDS_files/05_multivariate_taxon_sd_100.rds")%>% 
  as.data.frame() %>% 
  rownames_to_column("parameter") %>%
  filter(str_detect(parameter, "sigma")) %>%
  dplyr::select(parameter, SD = Mean)  %>% 
  mutate(ID = 
           str_extract(parameter, pattern = "\\[(.*?)\\,") %>% 
           str_remove_all(pattern = "[[:punct:]]") %>% 
           as.numeric(),
         trait = str_extract(parameter, pattern = "\\,(.*?)\\]") %>% 
           str_remove_all(pattern = "[[:punct:]]") %>% 
           as.factor() %>% 
           recode("1" = "SLA", 
                  "2" = "PlantHeight",
                  "3" = "LA", 
                  "4" = "LDMC",
                  "5" = "LeafC",
                  "6" = "LeafP",
                  "7" = "LeafN",
                  "8" = "Seed.length",
                  "9" = "Seed.count",
                  "10" = "Seed.mass")) %>% 
  left_join(pol, by = "ID") %>% 
  mutate(SD = exp(SD))


df <- left_join(dfMEAN,dfSD, by = c("pollentaxon", "trait")) %>% 
  mutate(Mean = round(Mean, 2),
         SD = round(SD,2),
         MeanSD = paste(Mean, "±", SD)) %>% 
  select(pollentaxon,trait, MeanSD) %>% 
  pivot_wider(., names_from = trait, values_from = MeanSD) %>% 
  mutate( LA = str_remove(LA, "\\.(.*?)[:space:]"),
          Seed.count = str_remove(Seed.count, "\\.(.*?)[:space:]")) 

write.xlsx(df, "Output/06_Supplement_trait_values.xlsx")

