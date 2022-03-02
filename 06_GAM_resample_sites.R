# This codes reruns the GAM models leaving out one site at the time, with the 
# goal of testing the sensivity to choice of sites.
# Annegreet Veeken

## load packages
library(tidyverse)
library(ggplot2)
library(rjags)
library(runjags)
library(mcmc)
library(coda)
library(mgcv) 
library(furrr)
library(gridExtra) # for combining plots
library(cowplot) # for get_legend

memory.limit(size = 9999999)

## Prepare data ----
dfCWM <- readRDS(paste0("RDS_files/05_multivariate_CWM_100.rds")) 
dfAGRI <- readRDS("RDS_files/Archaeological_indicators.rds") %>% 
  filter(age <= 10000) %>% 
  select(site.name, age, PresenceAgri) %>% 
  distinct()

dfTEMP <- readRDS("RDS_files/01_bio1_unbinned.rds") %>% 
  # add time label for joining
  mutate(timecat = factor(Time.BP, labels = as.character(1:length(unique(Time.BP)))))

lTRSH <- readRDS("RDS_files/03_TRSH_percentage.rds")
dfTRSH <- lTRSH %>% 
  purrr::map2(., names(lTRSH), 
              ~mutate(., site.name = .y, tree = round(adjustedpercent * 100, 1))) %>% 
  bind_rows() 
lCROP <- readRDS("RDS_files/03_Crop_percentage.rds")
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
                          "Seed.count", "Seed.length", "Seed.mass"),
                units = c("SLA (mg/mm2)",
                          "Height (m)",
                          "LDMC (g/g)",
                          "Leaf area (mm2)",
                          "Leaf C (mg/g)",
                          "Leaf N (mg/g)",
                          "Leaf P (mg/g)",
                          "Seed count",
                          "Seed length (mm)",
                          "Seed mass (mg)"),
                title = c("Specific leaf area",
                          "Plant height",
                          "Leaf dry matter content",
                          "Leaf area",
                          "Leaf carbon content",
                          "Leaf nitrogen content",
                          "Leaf phosphorus content",
                          "Seed number",
                          "Seed length",
                          "Seed mass"),
                ymin = c(5,1,0.2,1,380,10,1,1,1,1),
                ymax = c(30,30,0.5,4000,550,40,4,12000,15,1000)
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
           "Seed.count", "Seed.length", "Seed.mass","LDMC")
rm(dfAGRI,dfCROP,dfTEMP,lCROP,lTRSH,dfTRSH)

## GAM's
# GAM Mean ~ s(Time.BP) ----
if(0){
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
  saveRDS(jam, 
          paste0("RDS_files/06_GAM_resample_sites_", selectedtrait,"_",sitename, ".rds"))
  rm(jm, sam, jam)
}

plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[1], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[2], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[3], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[4], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[5], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[6], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[7], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[8], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[9], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(sitenames, ~gam_func(trait[10], sitename = .x), .options = furrr_options(seed = TRUE))
}

if(1){
# Plot  GAM1 ----
# extract fitted model from jam object
plot_gam <- function(selectedtrait) {
files <- list.files("RDS_files/") %>% 
  str_subset("06_GAM_resample_sites_") %>% 
  str_subset(paste0("_", selectedtrait))
files

folderpath.fun <- function(x)
{paste("RDS_files", x, sep = "/")}

jam_sites <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))

sitenames <- files %>% 
  str_remove(paste0("06_GAM_resample_sites_", selectedtrait, "_")) %>% 
  str_remove(".rds")

jam_int <- 
  purrr::map2(jam_sites, sitenames,
              ~tibble(intercept = coef(.x)[1],
                      sitename = .y)) %>% 
  bind_rows()

jam_values <- 
  purrr::map2(jam_sites, sitenames,
              ~plot(.) %>% 
                flatten %>% 
                keep(names(.) %in% c("x", "fit")) %>%
             # create df and add column with name of site that was left out
             bind_cols %>% 
              mutate(sitename = .y)) %>%
  bind_rows() %>% 
  left_join(jam_int, by = "sitename")

# plot fitted models
p <- ggplot(data = jam_values, aes(x = x, y = exp(fit+intercept), group = sitename)) +
  geom_line() +
  scale_y_continuous(paste(units$units[units$trait == selectedtrait]), trans = "log10",
                     limits = c(units$ymin[units$trait == selectedtrait],
                              units$ymax[units$trait == selectedtrait])) +
  scale_x_reverse("Time (calibrated years BP)") +
  ggtitle(units$title[units$trait == selectedtrait]) +
  theme_bw() +
  theme(axis.title.y =  element_text(size = 8),
        axis.title.x =  element_text(size = 8))
p
} 
windows()
q <- purrr::map(trait, ~plot_gam(.))
names(q) <- trait

qall <- grid.arrange(q$PlantHeight, q$SLA, q$LA, q$LDMC, q$LeafC, q$LeafN, q$LeafP,
                     q$Seed.count, q$Seed.length, q$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     widths = c(2,2),
                     nrow = 5, ncol = 2)

ggsave("Figures/SI5-GAM_time_resample_sites.png", qall, width = 174, 
       height = 247, units = "mm", dpi = 600)
}

# GAM Mean ~ s(since.agri) + s(temp) ----
if(0){
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
           years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),
                                -cumsum(replace_na(lag, 0)))) %>% 
    select(-lag)
  
  # Reference for code:
  # Simon N.Wood (2017). Generalized additive models. An introduction with R (2nd ed). page 374-376
  gam.model <- jagam(Mean ~ s(years.since) + s(temp), data = CWM, 
                     file = "gam_agri0_arrival_temp.jags",
                     diagonalize = TRUE) # priors on smooths Gaussian
  
  # changes from generated files are saved in gam_agri_time.jags
  # use previously calculated SD on CWM mean
  gam.model$jags.data$sd <- CWM$SD
  
  # add site name and site number for random effect
  gam.model$jags.data$site  <- as.numeric(as.factor(CWM$site.name))
  gam.model$jags.data$nsite <- length(unique(CWM$site.name))
  
  load.module("glm")
  jm <- jags.model("gam_agri_arrival_temp.jags", data = gam.model$jags.data,
                   inits = gam.model$inits, n.chains = 2, n.adapt = 500)
  sam <- jags.samples(jm, c("b", "rho"),  
                      n.iter = 10000, thin = 10)
  jam <- sim2jam(sam, gam.model$pregam)
  saveRDS(jam, 
          paste0("RDS_files/06_GAM2_resample_sites_", selectedtrait,"_",sitename, ".rds"))
  rm(jm, sam, jam)
}

agrisites <- sitenames[!sitenames %in% c("Bjärsjöholmssjön", "Changeon", "Claraghmore",
                        "Diheen", "Garry Bog", "Glen West", "Hares Down",
                        "Lough Henney", "Montbé", "Owenduff", "Paleochenal de Neublans",
                        "Port des Lamberts", "Silberhohl", "Sources de l'Yonne",
                        "Vladar", "Vrbka", "Windmill Rough", "Zahájí")]

plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[1], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[2], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[3], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[4], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[5], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[6], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[7], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[8], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[9], sitename = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agrisites, ~gam_func2(trait[10], sitename = .x), .options = furrr_options(seed = TRUE))
}

# Plot  GAM2 ----
if(1){
plot_gam2 <- function(selectedtrait){
files <- list.files("RDS_files/") %>% 
  str_subset("06_GAM2_resample_sites_") %>% 
  str_subset(paste0("_", selectedtrait))
files

sites <- files %>% 
  str_remove(paste0("06_GAM2_resample_sites_",selectedtrait, "_")) %>% 
  str_remove(".rds")

folderpath.fun <- function(x)
  {paste("RDS_files/", x, sep = "/")}

jam_sites <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))

jam_int2 <- 
  purrr::map2(jam_sites, sites,
              ~tibble(intercept = coef(.x)[1],
                      sitename = .y)) %>% 
  bind_rows()

# extract fitted model from jam object
jam_values2 <- purrr::map2(jam_sites, sites,
                          ~plot(., select = 1) %>% 
                            pluck(1) %>% 
                            keep(names(.) %in% c("x", "fit")) %>% 
                            # create df and add column with name of site that was left out
                            bind_cols %>% 
                            mutate(fit = as.numeric(fit),
                                   sitename = .y)) %>% 
  bind_rows() %>% 
  left_join(jam_int2, by = "sitename")

# plot fitted models
p <- ggplot(data = jam_values2, aes(x = x, y = exp(fit + intercept), group = sitename)) +
  geom_line() +
  scale_y_continuous(paste(units$units[units$trait == selectedtrait]),trans = "log10",
                     limits = c(units$ymin[units$trait == selectedtrait],
                                units$ymax[units$trait == selectedtrait])) +
  scale_x_continuous("Years since arrival of agriculture",
                     limits = c(-5000, 5000)) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  theme_bw() +
  theme(axis.title.y =  element_text(size = 8),
        axis.title.x =  element_text(size = 8))

p
}
windows()
q <- purrr::map(trait, ~plot_gam2(.))
names(q) <- trait
qall <- grid.arrange(q$PlantHeight, q$SLA, q$LA, q$LDMC, q$LeafC, q$LeafN, q$LeafP,
                     q$Seed.count, q$Seed.length, q$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     widths = c(2,2),
                     nrow = 5, ncol = 2)

ggsave("Figures/SI5-GAM_agri_resample_sites.png", qall, width = 174, 
       height = 247, units = "mm", dpi = 600)

# Temperature
plot_gam3 <- function(selectedtrait){
files <- list.files("RDS_files/") %>% 
  str_subset("06_GAM2_resample_sites_") %>% 
  str_subset(paste0("_", selectedtrait))
files

sites <- files %>% 
  str_remove(paste0("06_GAM2_resample_sites_",selectedtrait, "_")) %>% 
  str_remove(".rds")

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

jam_sites <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))

jam_int3 <- 
  purrr::map2(jam_sites, sites,
              ~tibble(intercept = coef(.x)[1],
                      sitename = .y)) %>% 
  bind_rows()

jam_values3 <- purrr::map2(jam_sites, sites,
                           ~plot(.,  select = 2) %>% pluck(2) %>% keep(names(.) %in% c("x", "fit")) %>%
                             # create df and add column with name of site that was left out
                             bind_cols %>% mutate(sitename = .y)) %>%
  bind_rows() %>% 
  left_join(jam_int3, by = "sitename")

# plot fitted models
p <- ggplot(data = jam_values3, aes(x = x, y = exp(fit+intercept), group = sitename)) +
  geom_line() +
  scale_y_continuous(paste(units$units[units$trait == selectedtrait]), trans = "log10",
                     limits = c(units$ymin[units$trait == selectedtrait],
                                units$ymax[units$trait == selectedtrait])) +
  scale_x_continuous("Temperature (°C)", limits = c(0, 13)) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  theme_bw() +
  theme(axis.title.y =  element_text(size = 8),
        axis.title.x =  element_text(size = 8))
p
}

windows()
p <- purrr::map(trait, ~plot_gam3(.))
names(p) <- trait
pall <- grid.arrange(p$PlantHeight, p$SLA, p$LA, p$LDMC, p$LeafC, p$LeafN, p$LeafP,
                     p$Seed.count, p$Seed.length, p$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     widths = c(2,2),
                     nrow = 5, ncol = 2)

ggsave("Figures/SI5-GAM_temp_resample_sites.png", pall, width = 174, 
       height = 247, units = "mm", dpi = 600)
}
graphics.off()