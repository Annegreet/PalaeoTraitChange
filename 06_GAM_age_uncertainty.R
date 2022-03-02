# load packages
library(tidyverse)
library(ggplot2)
library(rjags)
library(runjags)
library(mcmc)
library(coda)
library(mgcv) 
library(furrr)
library(gridExtra)

memory.limit(size = 99999999)

## Prepare data ----
dfCWM <- readRDS("RDS_files/05_multivariate_CWM_100.rds")  
dfAGRI <- readRDS("RDS_files/Archaeological_indicators.rds") %>% 
  filter(age <= 10000) %>% 
  select(site.name, age, PresenceAgri) %>% 
  distinct()
dfAGE_uncertainty <- readRDS("RDS_files/02_Calibrated_chronologies_uncertainties.rds") %>% 
  bind_rows() %>% 
  as_tibble()
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
  # create joining column for temperature data
  mutate(., timecat = cut(age, breaks = c(-Inf, unique(dfTEMP$Time.BP), Inf),
                          labels = c(1,1:length(unique(dfTEMP$Time.BP))))) %>% 
  left_join(dfTEMP, by = c("site.name", "timecat")) %>% 
  left_join(dfTRSH, by = c("site.name", "age", "depth")) %>%
  left_join(dfCROP, by = c("site.name", "age", "depth")) %>%
  dplyr::select(site.name, age, depth,trait,
                Mean, SD, PresenceAgri, temp, tree, crop) %>%
  mutate(crop = replace_na(crop, 0)) 

trait <- c("SLA","PlantHeight","LA", "LeafC", "LeafN", "LeafP", 
           "Seed.count", "Seed.length", "Seed.mass","LDMC")
selectedtrait <- "SLA"
rm(dfAGRI,dfCROP, dfCWM, dfTEMP, dfTRSH,lCROP, lTRSH)

## GAM's
# GAM Mean ~ s(Time.BP) ----
if(0){
gam_func <- function(agerandom, selectedtrait){
  # random select sample from age df
  age_sample <- dfAGE_uncertainty %>% 
    select(depth, site.name, dataset.id, age = all_of(agerandom))
  CWM <- CWMall %>% 
    filter(trait == selectedtrait) %>% 
    # remove median age
    select(-age) %>% 
    # add sampled age
    left_join(age_sample, by = c("site.name", "depth")) %>% 
    # sort time in descending order
    arrange(site.name, desc(age)) 
  
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
          paste0("RDS_files/06_GAM_age_uncertainty_", selectedtrait,"_",
                 agerandom,".rds"))

}
agerandom <- sample(5:1004, 50)
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[1], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[2], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[3], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[4], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[5], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[6], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[7], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[8], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[9], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func(trait[10], agerandom = .x), .options = furrr_options(seed = TRUE))
}

# Plot GAM1 ----
if(1){
plot_gam <- function(selectedtrait){
  files <- list.files("RDS_files/") %>%
    str_subset("06_GAM_age_uncertainty_") %>%
    str_subset(paste0("_", selectedtrait))

  ages <- files %>%
    str_remove(paste0("06_GAM_age_uncertainty_", selectedtrait, "_")) %>%
    str_remove(".rds")

  folderpath.fun <- function(x)
  {paste("RDS_files/", x, sep = "/")}

  jam_ages <- files %>%
    folderpath.fun(.) %>%
    purrr::map(~readRDS(.))
  
  jam_int <- 
    purrr::map2(jam_ages[1:50], ages[1:50],
                ~tibble(intercept = coef(.x)[1],
                        agedraw = .y)) %>% 
    bind_rows()

  # extract fitted model from jam object
  jam_values <- purrr::map2(jam_ages[1:50], ages[1:50],
                             ~plot(.) %>%
                               flatten(.) %>%
                               keep(names(.) %in% c("x", "fit")) %>%
                               # create df and add column with name of site that was left out
                               bind_cols %>%
                               mutate(fit = as.numeric(fit),
                                      agedraw = .y)) %>%
    bind_rows() %>% 
    left_join(jam_int, by = "agedraw")

  # plot fitted models
  p <- ggplot(data = jam_values, aes(x = x, y = exp(fit + intercept), group = agedraw)) +
    geom_line(col = "darkgreen") +
    scale_y_continuous(paste(units$units[units$trait == selectedtrait]), trans = "log10",
                       limits = c(units$ymin[units$trait == selectedtrait],
                                  units$ymax[units$trait == selectedtrait])) +
    scale_x_reverse("Time (calibrated years BP)",limits = c(10000,0)) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    theme_bw() +
    theme(axis.title.y =  element_text(size = 8),
          axis.title.x =  element_text(size = 8))
  p
}

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

ggsave("Figures/SI5-GAM_time_age_uncertainty.png", qall, width = 174, 
       height = 247, units = "mm", dpi = 600)
}

# GAM 2 Mean ~ s(since.agri) + s(temp) ----
if(0){
gam_func2 <- function(selectedtrait, agerandom){
  age_sample <- dfAGE_uncertainty %>% 
    select(depth, site.name, dataset.id, age = all_of(agerandom))
  CWM <- CWMall %>% 
    filter(trait == selectedtrait) %>% 
    # remove median age
    select(-age) %>% 
    # add sampled age
    left_join(age_sample, by = c("site.name", "depth")) %>%
    #filter sites with agriculture already present
    filter(!site.name %in% c("Bjärsjöholmssjön", "Changeon", "Claraghmore",
                             "Diheen", "Garry Bog", "Glen West", "Hares Down",
                             "Lough Henney", "Montbé", "Owenduff", "Paleochenal de Neublans",
                             "Port des Lamberts", "Silberhohl", "Sources de l'Yonne",
                             "Vladar", "Vrbka", "Windmill Rough", "Zahájí")) %>%
    filter(trait == selectedtrait) %>% 
    # create column with time since arrival of agriculture
    group_by(site.name, PresenceAgri) %>% 
    mutate(lag = ifelse(PresenceAgri == "After agriculture", age - lead(age), age - lag(age)),
           years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),
                                -cumsum(replace_na(lag, 0)))) 
  
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
saveRDS(jam, file = paste("RDS_files/06_GAM2_age_uncertainty_",
                          selectedtrait,"_",
                          agerandom, ".rds", sep = ""))
rm(jm, sam,jam)
}

agerandom <- sample(5:1004, 50)
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[1], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[2], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 2))
furrr::future_map(agerandom, ~gam_func2(trait[3], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[4], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[5], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[6], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[7], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[8], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[9], agerandom = .x), .options = furrr_options(seed = TRUE))
plan(multisession(workers = 4))
furrr::future_map(agerandom, ~gam_func2(trait[10], agerandom = .x), .options = furrr_options(seed = TRUE))
}

# Plot GAM2 -----
if(1){
plot_gam2 <- function(selectedtrait){
  files <- list.files("RDS_files/") %>%
    str_subset("06_GAM2_age_uncertainty_") %>%
    str_subset(paste0("_", selectedtrait))

 ages <- files %>%
    str_remove(paste0("06_GAM2_age_uncertainty_", selectedtrait, "_")) %>%
    str_remove(".rds")
  
  folderpath.fun <- function(x)
  {paste("RDS_files/", x, sep = "/")}
  
  jam_ages <- files %>%
    folderpath.fun(.) %>%
    purrr::map(~readRDS(.))
  
  jam_int2 <- 
    purrr::map2(jam_ages[1:50], ages[1:50],
                ~tibble(intercept = coef(.x)[1],
                        agedraw = .y)) %>% 
    bind_rows()
  
  # extract fitted model from jam object
  jam_values2 <- purrr::map2(jam_ages[1:50], ages[1:50],
                             ~plot(., select = 1) %>% 
                               pluck(1) %>% 
                               keep(names(.) %in% c("x", "fit")) %>% 
                               # create df and add column with name of site that was left out
                               bind_cols %>% 
                               mutate(fit = as.numeric(fit),
                                      agedraw = .y)) %>% 
    bind_rows()  %>% 
    left_join(jam_int2, by = "agedraw")
  
  # plot fitted models
  p <- ggplot(data = jam_values2, aes(x = x, y = exp(fit + intercept), group = agedraw)) +
    geom_line(col = "darkgreen") +
    scale_y_continuous(paste(units$units[units$trait == selectedtrait]), trans = "log10",
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

ggsave("Figures/SI5-GAM_agri_age_uncertainty.png", qall, width = 174, 
       height = 247, units = "mm", dpi = 600)

plot_gam3 <- function(selectedtrait){
  files <- list.files("RDS_files/") %>%
    str_subset("06_GAM2_age_uncertainty_") %>%
    str_subset(paste0("_", selectedtrait))
  
  ages <- files %>%
    str_remove(paste0("06_GAM2_age_uncertainty_", selectedtrait, "_")) %>%
    str_remove(".rds")
  
  folderpath.fun <- function(x)
  {paste("RDS_files/", x, sep = "/")}
  
  jam_ages <- files %>%
    folderpath.fun(.) %>%
    purrr::map(~readRDS(.))
  
  jam_int3 <- 
    purrr::map2(jam_ages[1:50], ages[1:50],
                ~tibble(intercept = coef(.x)[1],
                        agedraw = .y)) %>% 
    bind_rows()
  
  jam_values3 <- purrr::map2(jam_ages[1:50], ages[1:50],
                             ~plot(.,  select = 2) %>% pluck(2) %>% keep(names(.) %in% c("x", "fit")) %>%
                               # create df and add column with name of site that was left out
                               bind_cols %>% mutate(agedraw = .y)) %>%
    bind_rows() %>% 
    left_join(jam_int3, by = "agedraw")
  
  # plot fitted models
  p <- ggplot(data = jam_values3, aes(x = x, y = exp(fit+intercept), group = agedraw)) +
    geom_line(col = "darkgreen") +
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
l <- purrr::map(trait, ~plot_gam3(.))
names(l) <- trait
lall <- grid.arrange(l$PlantHeight, l$SLA, l$LA, l$LDMC, l$LeafC, l$LeafN, l$LeafP,
                     l$Seed.count, l$Seed.length, l$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     widths = c(2,2),
                     nrow = 5, ncol = 2)

ggsave("Figures/SI5-GAM_temp_age_uncertainty.png", lall, width = 174, 
       height = 247, units = "mm", dpi = 600)
}
graphics.off()
