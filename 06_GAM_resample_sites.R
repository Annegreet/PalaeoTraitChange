# load packages
library(tidyverse)
library(ggplot2)
library(rjags)
library(runjags)
library(mcmc)
library(coda)
library(mgcv) 
library(gridExtra) # for combining plots
library(cowplot) # for get_legend

## Prepare data ----
dfCWM <- readRDS("RDS_files/05_Posterior_alltraits_unbinned.rds")
dfAGRI <- readRDS("RDS_files/Archaeological_indicators.rds") %>% 
  filter(age <= 10000) %>% 
  select(site.name, age, PresenceAgri) %>% 
  distinct()

dfTEMP <- readRDS("RDS_files/01_bio1_unbinned.rds") %>% 
  # add time label for joining
  mutate(timecat = factor(Time.BP, labels = as.character(1:length(unique(Time.BP)))))

lTRSH <- readRDS("RDS_files/03_TRSH_percentage_unbinned.rds")
dfTRSH <- lTRSH %>% 
  purrr::map2(., names(lTRSH), 
              ~mutate(., site.name = .y, tree = round(adjustedpercent * 100, 1))) %>% 
  bind_rows() 
lCROP <- readRDS("RDS_files/03_Crop_percentage_unbinned.rds")
dfCROP <- lCROP %>% 
  purrr::map2(., names(lCROP), 
              ~mutate(., site.name = .y, crop = round(adjustedpercent * 100, 1))) %>% 
  bind_rows() 

# colour palette
cbf <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

selectedtrait <- "SLA"

units <- tibble(trait = c("SLA", "PlantHeight","LDMC", 
                          "LA", "LeafC", "LeafN", "LeafP", 
                          "Seed.count", "Seed.lenght", "Seed.mass"),
                units = c("Community SLA (mg/mm2)",
                          "Community height (m)",
                          "Community LDMC (g/g)",
                          "Community leaf area (mm2)",
                          "Community leaf C (mg/g)",
                          "Community leaf N (mg/g)",
                          "Community leaf P (mg/g)",
                          "Community seed count",
                          "Community seed length (mm)",
                          "Community seed mass (mg)"),
                title = c("Specific leaf area",
                          "Plant height",
                          "Leaf dry matter content",
                          "Leaf area",
                          "Leaf carbon content",
                          "Leaf nitrogen content",
                          "Leaf phosphorus content",
                          "Seed number",
                          "Seed length",
                          "Seed mass")
)

sitenames <- unique(dfCWM$site.name)

# merge data sets
CWMall <- dfCWM %>% 
  left_join(., dfAGRI, by = c("site.name", "age")) %>% 
  mutate(site.name = as.factor(site.name)) %>% 
  mutate(., timecat = cut(age, breaks = c(-Inf, unique(dfTEMP$Time.BP), Inf),
                          labels = c(1,1:length(unique(dfTEMP$Time.BP))))) %>% 
  left_join(dfTEMP, by = c("site.name", "timecat")) %>% 
  left_join(dfTRSH, by = c("site.name", "age")) %>%
  left_join(dfCROP, by = c("site.name", "age")) %>%
  dplyr::select(site.name, age, trait,
                Mean, SD, PresenceAgri, temp, tree, crop) %>%
  mutate(crop = replace_na(crop, 0)) %>% 
  as_tibble()

trait <- c("SLA","PlantHeight","LA", "LeafC", "LeafN", "LeafP", 
           "Seed.count", "Seed.lenght", "Seed.mass","LDMC")
selectedtrait <- "LA"

## GAM's
# GAM Mean ~ s(Time.BP) ----
gam_func <- function(selectedtrait, sitename){
  CWM <- CWMall %>% 
    filter(trait == selectedtrait) %>% 
    # sort time.bp in descending order
    arrange(site.name, desc(age)) %>% 
    # filter out one site at the time
    filter(!site.name == sitename) %>%  
    as_tibble()
  
  # Reference for code:
  # Simon N.Wood (2017). Generalized additive models. An introduction with R (2nd ed). page 374-376
  gam.model <- jagam(Mean ~ s(age), data = CWM, 
                     file = "gam_agri0_time.jags",
                     diagonalize = TRUE) # priors on smooths Gaussian
  
  # changes from generated files are saved in gam_agri_time.jags
  # use previously calculated SD on CWM mean
  gam.model$jags.data$sd <- CWM$SD
  
  # add site name and site number for random effect
  gam.model$jags.data$site <- as.numeric(as.factor(CWM$site.name))
  gam.model$jags.data$nsite <- length(unique(CWM$site.name))
  
  load.module("glm")
  jm <- jags.model("gam_agri_time.jags", data = gam.model$jags.data,
                   inits = gam.model$inits, n.chains = 2, n.adapt = 500)
  sam <- jags.samples(jm, c("b", "rho"),
                      n.iter = 10000, thin = 10) # set to n.iter = 10000 thin = 10
  jam <- sim2jam(sam, gam.model$pregam)
jam
}

jam_sites <- purrr::map(sitenames, ~gam_func(selectedtrait, sitename = .x))
saveRDS(jam_sites, 
        paste0("RDS_files/06_GAM_resample_sites_", selectedtrait, ".rds"))


# Plot  GAM's ----

# extract fitted model from jam object
jam_values <- purrr::map2(jam_sites, sitenames,
           ~plot(.) %>% flatten %>% keep(names(.) %in% c("x", "fit")) %>% 
             # create df and add column with name of site that was left out
             bind_cols %>% mutate(sitename = .y)) %>% 
  bind_rows()

# plot fitted models
p <- ggplot(data = jam_values, aes(x = x, y = fit, group = sitename)) +
  geom_line() +
  scale_y_continuous(paste(units$units[units$trait == selectedtrait], "(log)")) +
  scale_x_reverse("Time (calibrated years BP)") +
  ggtitle(units$title[units$trait == selectedtrait]) +
  theme_bw() 
ggsave(paste0("Figures/S5_Time-Resampled-sites_", selectedtrait, ".rds"), p) 

# GAM Mean ~ s(since.agri) + s(temp) ----
gam_func2 <- function(selectedtrait, sitename){
  CWM <- CWMall %>% 
    #filter sites with agriculture already present
    filter(!site.name %in% c("Bjärsjöholmssjön", "Changeon", "Claraghmore",
                             "Diheen", "Garry Bog", "Glen West", "Hares Down",
                             "Lough Henney", "Montbé", "Owenduff", "Paleochenal de Neublans",
                             "Port des Lamberts", "Silberhohl", "Sources de l'Yonne",
                             "Vladar", "Vrbka", "Windmill Rough", "Zahájí")) %>%

    filter(trait == selectedtrait) %>% 
    # filter out one site at the time
    filter(!site.name == sitename) %>%  
    # create column with time since arrival of agriculture
    group_by(site.name, PresenceAgri) %>% 
    mutate(lag = ifelse(PresenceAgri == "After agriculture", age - lead(age), age - lag(age)),
           years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),0)) %>% 
    select(-lag)
  
  # Reference for code:
  # Simon N.Wood (2017). Generalized additive models. An introduction with R (2nd ed). page 374-376
  gam.model <- jagam(Mean ~ s(years.since) + s(temp), data = CWM, 
                     file = "gam_agri0_arrival_temp.jags",
                     diagonalize = TRUE) # priors on smooths Gaussian
  
  # changes from generated files are saved in gam_agri_time.jags
  # use previously calculated SD on CWM mean
  gam.model$jags.data$sd <- dfCWM$SD
  
  # add site name and site number for random effect
  gam.model$jags.data$site  <- as.numeric(as.factor(CWM$site.name))
  gam.model$jags.data$nsite <- length(unique(CWM$site.name))
  
  load.module("glm")
  jm <- jags.model("gam_agri_arrival_temp.jags", data = gam.model$jags.data,
                   inits = gam.model$inits, n.chains = 2, n.adapt = 500)
  sam <- jags.samples(jm, c("b", "rho"),  
                      n.iter = 10000, thin = 10)
  jam <- sim2jam(sam, gam.model$pregam)
  jam
}

agrisites <- sitenames[!sitenames %in% c("Bjärsjöholmssjön", "Changeon", "Claraghmore",
                        "Diheen", "Garry Bog", "Glen West", "Hares Down",
                        "Lough Henney", "Montbé", "Owenduff", "Paleochenal de Neublans",
                        "Port des Lamberts", "Silberhohl", "Sources de l'Yonne",
                        "Vladar", "Vrbka", "Windmill Rough", "Zahájí")]

jam_sites2 <- purrr::map(agrisites, ~gam_func2(selectedtrait, sitename = .x))
saveRDS(jam_sites2, 
        paste0("RDS_files/06_GAM_agri_resample_sites_", selectedtrait, ".rds"))

# Plot  GAM's agriculture ----

# extract fitted model from jam object
jam_values2 <- purrr::map2(jam_sites, agrisites,
                          ~plot(., select(1)) %>% pluck(1) %>% keep(names(.) %in% c("x", "fit")) %>% 
                            # create df and add column with name of site that was left out
                            bind_cols %>% mutate(sitename = .y)) %>% 
  bind_rows()

# plot fitted models
p <- ggplot(data = jam_values2, aes(x = x, y = fit, group = sitename)) +
  geom_line() +
  scale_y_continuous(paste(units$units[units$trait == selectedtrait], "(log)")) +
  scale_x_continuous("Years since arrival of agriculture",
                     limits = c(-5000, 5000)) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  theme_bw() 
ggsave(paste0("Figures/S5_Agriculture-Resampled-sites_", selectedtrait, ".rds"), p) 
