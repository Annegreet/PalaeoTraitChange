library(tidyverse)
library(runjags)
library(mcmc)
library(coda)
library(xlsx)


## Pollen data
lPOL <- readRDS("RDS_files/04_PollenWEU-PPEcorrected.rds")
pollentaxa <- lPOL %>% purrr::map(., ~.$ppe.taxon) %>% unlist %>% unique()

files <- 
  list.files("RDS_files/Posterior_CWM/") %>% 
  str_subset(pattern = "05_CWM_estimates_")
files <- files[!str_detect(files, pattern = "Stats")]
files <- files[!str_detect(files, pattern = "herbs")]
files <- files[!str_detect(files, pattern = "trees")]

folderpath.fun <- function(x)
{paste("RDS_files/Posterior_CWM/", x, sep = "/")}

traitname <- files %>% 
  str_remove(., "05_CWM_estimates_") %>% 
  str_remove(., ".rds")

df <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(traitname, ~readRDS(.) %>% 
                rownames_to_column("parameter") %>% 
                filter(., str_detect(parameter, "mean.tax")) %>% 
                dplyr::select(Mean, SD) %>%  
                mutate(trait = .y, pollentaxon = pollentaxa,
                       value = paste(round(exp(Mean), digits = 1),"±",
                                     round(exp(SD), digits = 2)))) %>% 
  bind_rows() %>%
  select(-Mean, -SD) %>% 
  pivot_wider(., names_from = trait)

write.xlsx(df, "Output/06_Supplement_trait_values.xlsx")

x <- readRDS("RDS_files/Posterior_CWM/05_Stats_CWM_estimates_LA.rds")
y <- readRDS("RDS_files/Posterior_CWM/05_CWM_estimates_Seed.count.rds")

yy <- as.mcmc(y)

gelman.plot(yy)

plot(x, vars =  "sd.tax")

