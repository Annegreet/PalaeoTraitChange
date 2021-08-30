library(xlsx)
library(tidyverse)

## Summary table 1 for manuscript
## Species per pollen type list
spec <- readRDS("RDS_files/01_PollenType_species.rds") %>% as_tibble()
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = round(ppe$`RPP v1 (Northern Hemisphere)`,2))

# filter out aquatic species
spec <- spec %>% 
  filter(StdValue < 10 | is.na(StdValue)) 
spec$pollentaxon[spec$pollentaxon == "Sambucus nigra"] <- "Sambucus.nigra"
spec <- spec %>% 
  dplyr::select(stand.spec, pollentaxon) %>% 
  filter(!pollentaxon == "NA" & !is.na(stand.spec)) %>% 
  distinct()

# Gap-filled trait data
gapdata <- readRDS("RDS_files/Gap_filled_traits_2020.rds") 
gapdata <- gapdata %>% filter(Species %in% spec$stand.spec) # only the relevant species

# merge trait data with spec/pol list
trait <- gapdata %>% 
  inner_join(spec, by = c("Species" = "stand.spec")) %>% 
  arrange(pollentaxon) # sort alphabetically

# Pollen data
lPOL <- readRDS("RDS_files/04_PollenWEU-PPEcorrected.rds")
pollentaxa <- lPOL %>% 
  purrr::map(., ~.$ppe.taxon) %>% unlist %>% unique
sitenames <- names(lPOL)
# occurrence of pollen taxa in records
polper <- lPOL %>% 
  purrr::map2(., sitenames, ~mutate(.x, site.name = .y)) %>% 
  bind_rows() %>% 
  group_by(site.name, ppe.taxon) %>% 
  filter(adjustedpercent > 0) %>% 
  summarise(n = n()) %>% 
  mutate(present = ifelse(n > 0, 1, 0)) %>% 
  group_by(ppe.taxon) %>% 
  summarise(occ = round(sum(present)/length(sitenames)*100, 0))

# unique species in pollen taxon
sumspec <-
  spec %>%
  group_by(pollentaxon) %>% 
  filter(pollentaxon %in% pollentaxa) %>% 
  summarise(unique.spec = length(unique(stand.spec))) 
# number of trait observations per pollen taxon
nobs <- trait %>% 
  group_by(pollentaxon) %>% 
  filter(pollentaxon %in% pollentaxa) %>% 
  summarise(nobs = n())

t1 <- left_join(sumspec, nobs, by = "pollentaxon") %>% 
  left_join(ppe, by = c("pollentaxon" = "ppe.taxon")) %>% 
  left_join(polper, by = c("pollentaxon" = "ppe.taxon"))


write.xlsx(t1,"Output/05_Summary_table_pollentypes.xlsx")
