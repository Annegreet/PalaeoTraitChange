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
## GAM's
# GAM Mean ~ s(Time.BP) ----
gam_func <- function(selectedtrait){
  CWM <- CWMall %>% 
    filter(trait == selectedtrait) %>% 
    # sort time.bp in descending order
    arrange(site.name, desc(age)) %>% 
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
                      n.iter = 10000, thin = 10)
  jam <- sim2jam(sam, gam.model$pregam)
  coda.sam <- coda.samples(jm, c("b", "rho"),
                           n.iter = 10000, thin = 10)
  
  saveRDS(jam, file = paste("RDS_files/06_GAM_time_unbinned_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(sam, file = paste("RDS_files/06_GAM_time_sam_unbinned_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(coda.sam, file = paste("RDS_files/06_GAM_time_coda_unbinned_",
                                 selectedtrait, ".rds", sep = ""))
}

purrr::map(trait, ~gam_func(.))

# Plot  GAM's ----
# create list of jam's
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_time_unbinned_") %>% 
  str_subset(pattern = "sam", negate = TRUE) %>% 
  str_subset(pattern = "coda", negate = TRUE) 
files

traitname <- files %>% 
  str_remove(., "06_GAM_time_unbinned_") %>% 
  str_remove(., ".rds")

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

jam_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(jam_list) <- traitname

# create list of sam objects
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_time_sam_unbinned_")  
files

sam_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(sam_list) <- traitname

# create list of coda objects
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_time_coda_unbinned_") 
files

coda_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(coda_list) <- traitname

# check convergence
x <- 10
gelman.plot(coda_list[[x]], ask = FALSE)
plot(coda_list[[x]])

# quick plot smooths for all traits
plot_func <- function(x,y) {
  plot(x,  shade = TRUE, xlab = "Time (cal BP)",
       ylab =  units$units[units$trait == y],
       xlim = c(10000, 0), 
       #ylim = c(0, max(CWMall$Mean[CWMall$trait == y])+quantile(CWMall$Mean[CWMall$trait == y],0.25)),
       main = units$title[units$trait == y],
       shift = coef(x)[1])
  # points(CWMall$Time.BP[CWMall$trait == y], 
  #        CWMall$Mean[CWMall$trait == y],
  #        pch = 20)
  # abline(h=0)
}

# png(filename = "Figures/06_GAM/gam_time.png", width = 1024, height = 544)
windows()
par(mfrow = c(2,5))
purrr::map2(jam_list, names(jam_list), 
            ~plot_func(.x,.y))
# dev.off()

# ggplot -----
gg_func <- function(jam, selectedtrait) {
  values <- plot(jam) %>% flatten
  intercept <- coef(jam)[1]
  df <- data.frame(age = values$x, Mean = intercept + values$fit, se = values$se)
  ggplot(data = CWMall[CWMall$trait == selectedtrait,], 
         aes(x = age, y = exp(Mean))) +
    geom_point(aes(colour = tree),  size = .5) +
    scale_color_gradient("Tree pollen\n(%)", low = cbf[5],
                         high = cbf[4], limits = c(0, 100)) +
    scale_y_continuous(units$units[units$trait == selectedtrait]) +
     scale_x_reverse("Time (calibrated years BP)", limits = c(10000,0)) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    geom_line(data = df, aes(x = age, y = exp(Mean)), colour = "#0072B2",
              size = 1) +
    geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
                aes(ymin = exp(Mean - se), ymax = exp(Mean + se))) + 
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 10))
}

windows()
p <- purrr::map2(jam_list, names(jam_list), 
                 ~gg_func(.x,.y))

selectedtrait <- "SLA"
traitindex <- which(names(jam_list) %in% selectedtrait)

windows()
values <- plot(jam_list[[traitindex]]) %>% flatten
intercept <- coef(jam_list[[traitindex]])[1]
df <- data.frame(age = values$x, Mean = intercept + values$fit, se = values$se)

p1 <- ggplot(data = CWMall[CWMall$trait == selectedtrait,], 
             aes(x = age, y = exp(Mean))) +
  geom_point(aes(colour = tree)) +
  scale_color_gradient("Tree pollen\n(%)", low = cbf[5],
                       high = cbf[4], limits = c(0, 100)) +
  scale_y_continuous(units$units[units$trait == selectedtrait]) +
  scale_x_reverse("Time (calibrated years BP)", limits = c(10000,0)) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  geom_line(data = df, aes(x = age, y = exp(Mean)), colour = "#0072B2",
            size = 1) +
  geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
              aes(ymin = exp(Mean - se), ymax = exp(Mean + se))) + 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 10))
p1

selectedtrait <- "LDMC"
traitindex <- which(names(jam_list) %in% selectedtrait)

windows()
values <- plot(jam_list[[traitindex]]) %>% flatten
intercept <- coef(jam_list[[traitindex]])[1]
df <- data.frame(age = values$x, Mean = intercept + values$fit, se = values$se)

p2 <- ggplot(data = CWMall[CWMall$trait == "LDMC",], 
             aes(x = age, y = Mean)) +
  geom_point(aes(colour = tree), size = .5) +
  scale_color_gradient("Tree pollen (%)", low = cbf[5], 
                       high = cbf[4], limits = c(0, 100)) +
  scale_y_continuous(units$units[units$trait == "LDMC"],
                     limits = c(0.2,0.5)) +
  scale_x_reverse("Time (calibrated years BP)") +
  ggtitle(units$title[units$trait == "LDMC"]) +
  geom_line(data = df, aes(x = age, y = Mean), colour = "#0072B2",
            size = 1) +
  geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
              aes(ymin = Mean - se, ymax = Mean + se)) + 
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 10))
p2

mylegend <- get_legend(p1)

# set limits to be the same across GAM's
p$LeafC <- p$LeafC +  scale_y_continuous(units$units[units$trait == "LeafC"],limits = c(425, 560)) 

# increase margins
p$PlantHeight <- p$PlantHeight + theme(plot.margin = grid::unit(c(0,0,0,1),"cm"))
p$SLA <- p$SLA + theme(plot.margin = grid::unit(c(0,0.5,0,0.5),"cm"))
p$LA <- p$LA + theme(plot.margin = grid::unit(c(0,0,0,1),"cm"))
p2 <- p2 + theme(plot.margin = grid::unit(c(0,0.5,0,0.5),"cm"))
p$LeafC <- p$LeafC + theme(plot.margin = grid::unit(c(0,0,0,1),"cm"))
p$LeafN <- p$LeafN + theme(plot.margin = grid::unit(c(0,0.5,0,0.5),"cm"))
p$LeafP <- p$LeafP + theme(plot.margin = grid::unit(c(0,0,0,1),"cm"))
p$Seed.count <- p$Seed.count + theme(plot.margin = grid::unit(c(0,0.5,0,0.5),"cm"))
p$Seed.lenght <- p$Seed.lenght + theme(plot.margin = grid::unit(c(0,0,0,1),"cm"))
p$Seed.mass <- p$Seed.mass + theme(plot.margin = grid::unit(c(0,0.5,0,0.5),"cm"))

windows()
pall <- grid.arrange(p$PlantHeight, p$SLA, p$LA, p2, p$LeafC,
                     p$LeafN, p$LeafP, p$Seed.count, p$Seed.lenght, p$Seed.mass,
                    layout_matrix = rbind(c(1,2,NA),
                                          c(3,4,NA),
                                          c(5,6,11),
                                          c(7,8,NA),
                                          c(9,10,NA)),
                    mylegend, nrow = 5, ncol = 3, widths = c(2,2,0.5))

ggsave("Figures/Fig5-GAM_time.pdf", pall, width = 174, 
       height = 247, units = "mm", dpi = 600)

# draws from posteriors for all traits----
gam_samples <- function(trait){
  pd <- data.frame(age = seq(0, 10000, by = 500))
  Xp <- predict(jam_list[[trait]], newdata = pd, type = "lpmatrix")
  ii <- 1:25 * 20 + 500
  
  for (i in ii) {
    fv <- Xp %*% sam_list[[trait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd$age, fv, type = "l", lty = 3, xlab = "Time (cal BP)", 
           ylab = paste0(units$units[units$trait == trait],"(log)"), xlim = c(10000, 0),
           main = units$title[units$trait == trait],
           ylim = c(ylim$ymin[units$trait == trait], ylim$ymax[ylim$trait == trait]))
    else
      lines(pd$age, fv, lty = 3)
  }
}
ylim <- tibble(trait = c("SLA", "PlantHeight","LDMC", 
                          "LA", "LeafC", "LeafN", "LeafP", 
                          "Seed.count", "Seed.lenght", "Seed.mass"),
                ymin = c(1.8,  0.3, 0.25, 5.4,  6, 2.5, 0.25, 6, 1, 0),
                ymax = c(2.75, 2, 0.5, 7.2, 6.3, 3.15, 0.6, 8, 1.9, 2.8))
# quick plot
png("Figures/SI4-Draws-time.png", width = 174 , height = 300,units = "mm", res = 300)
par(mfrow = c(5,2))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
              "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.lenght", "Seed.mass"), ~gam_samples(.))
dev.off()

# GAM Mean ~ s(since.agri) + s(temp) ----
gam_func2 <- function(selectedtrait){
  CWM <- CWMall %>% 
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
  coda.sam <- coda.samples(jm, c("b", "rho"),
                           n.iter = 10000, thin = 10)
  
  saveRDS(jam, file = paste("RDS_files/06_GAM_arrival_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(sam, file = paste("RDS_files/06_GAM_arrival_sam_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(coda.sam, file = paste("RDS_files/06_GAM_arrival_coda_",
                                 selectedtrait, ".rds", sep = ""))
}

purrr::map(trait, ~gam_func2(.))

# Compile results to list
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_arrival2_") %>% 
  str_subset(pattern = "sam", negate = TRUE) %>% 
  str_subset(pattern = "coda", negate = TRUE) 
files
traitname <- files %>% 
  str_remove(., "06_GAM_arrival2_") %>% 
  str_remove(., ".rds")

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

jam_list2 <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(jam_list2) <- traitname

files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_arrival2_sam") 
files

sam_list2 <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(sam_list2) <- traitname

files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_arrival2_coda")
files

coda_list2 <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(coda_list2) <- traitname

# check convergence
x <- 10
plot(coda_list2[[x]])
gelman.plot(coda_list2[[x]], ask = FALSE)

# plot smooths for all traits
# agriculture plots ----
# time since arrival of agriculture
plot_func_agri <- function(x,y) {
  plot(x,  shade = TRUE, select = 1,
       scale = 0, xlab = "Years since arrival of agriculture",
       xlim = c(0,5000),
       main = units$title[units$trait == y])
  abline(h=0)
}

windows()
par(mfrow = c(2,5))
purrr::map2(jam_list2, names(jam_list2), 
            ~plot_func_agri(.x,.y))

gg_func_agri <- function(jam, selectedtrait) {
  CWM <- CWMall %>% 
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
  
  # extract smooth for plotting
  intercept <- coef(jam)[1]
  values <-  plot(jam, select = 1) %>% pluck(1)
  df <- tibble(Time = values$x, Mean = values$fit + intercept,
               se = values$se)
  
  ggplot(data = CWM, aes(x = years.since, y = exp(Mean))) +
    # plot residuals
    geom_point(aes(colour = crop), size = 0.5) +
    scale_color_gradient("Crop pollen (%)", low = cbf[5], high = cbf[4],
                         limits = c(0.01,100), trans = "log",breaks = c(0.1,2.5,50),
                         labels = c(0.1,2.5,50), na.value = "#999999") +
    geom_line(data = df, aes(x = Time, y = exp(Mean)), size = 1.5, 
              colour = "#0072B2") +
    geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
                aes(x = Time, y = exp(Mean), ymin = exp(Mean - se), 
                    ymax = exp(Mean + se))) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    
    # plot smooth 
    scale_y_continuous(units$units[units$trait == selectedtrait]) +
    scale_x_continuous("Years since arrival of agriculture",
                       limits = c(-5000, 5000)) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 10))
}

windows()
p <- purrr::map2(jam_list2, names(jam_list2), 
                 ~gg_func_agri(.x,.y))

# plot for extracting legend
selectedtrait <- "SLA"
CWM <- CWMall %>% 
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

# extract smooth for plotting
intercept <- coef(jam_list2[[selectedtrait]])[1]
values <-  plot(jam_list2[[selectedtrait]], select = 1) %>% pluck(1)
df <- tibble(Time = values$x, Mean = values$fit + intercept,
             se = values$se)

p1 <- ggplot(data = CWM, aes(x = years.since, y = exp(Mean))) +
  # plot residuals
  geom_point(aes(colour = crop), size = 0.5) +
  scale_color_gradient("Crop pollen\n(%)", low = cbf[5], high = cbf[4],
                       limits = c(0.01,100), trans = "log",breaks = c(0.1,2.5,50),
                       labels = c(0.1,2.5,50), na.value = "#999999") +
  geom_line(data = df, aes(x = Time, y = exp(Mean)), size = 1.5, 
            colour = "#0072B2") +
  geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
              aes(x = Time, y = exp(Mean), ymin = exp(Mean - se), 
                  ymax = exp(Mean + se))) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  
  # plot smooth 
  scale_y_continuous(units$units[units$trait == selectedtrait]) +
  scale_x_continuous("Years since arrival of agriculture",
                     limits = c(-5000, 5000)) +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 10))

mylegend <- get_legend(p1)

# plot LDMC separately (no exp transformation)
CWM <- CWMall %>% 
  #filter sites with agriculture already present
  filter(!site.name %in% c("Bjärsjöholmssjön", "Changeon", "Claraghmore",
                           "Diheen", "Garry Bog", "Glen West", "Hares Down",
                           "Lough Henney", "Montbé", "Owenduff", "Paleochenal de Neublans",
                           "Port des Lamberts", "Silberhohl", "Sources de l'Yonne",
                           "Vladar", "Vrbka", "Windmill Rough", "Zahájí")) %>%
  filter(trait == "LDMC") %>% 
  # create column with time since arrival of agriculture
  group_by(site.name, PresenceAgri) %>% 
  mutate(lag = ifelse(PresenceAgri == "After agriculture", age - lead(age), age - lag(age)),
         years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),
                              -cumsum(replace_na(lag, 0)))) 

# extract smooth for plotting
intercept <- coef(jam_list2[["LDMC"]])[1]
values <-  plot(jam_list2[["LDMC"]], select = 1) %>% pluck(1)
df <- tibble(Time = values$x, Mean = values$fit + intercept,
             se = values$se)
df

p2 <- ggplot(data = CWM, aes(x = years.since, y = Mean)) +
  # plot residuals
  geom_point(aes(colour = crop), size = 0.5) +
  scale_color_gradient("Crop pollen (%)", low = cbf[5], high = cbf[4],
                       limits = c(0.01,100), trans = "log",breaks = c(0.1,2.5,50),
                       labels = c(0.1,2.5,50), na.value = "#999999") +
  geom_line(data = df, aes(x = Time, y = Mean), size = 1.5, 
            colour = "#0072B2") +
  geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
              aes(x = Time, y = Mean, ymin = Mean - se, 
                  ymax = Mean + se)) +
  ggtitle(units$title[units$trait == "LDMC"]) +
  
  # plot smooth 
  scale_y_continuous(units$units[units$trait == "LDMC"]) +
  scale_x_continuous("Years since arrival of agriculture",
                     limits = c(-5000,5000)) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 10))

p2

# Set common scale
p$LeafC <- p$LeafC +  scale_y_continuous(units$units[units$trait == "LeafC"],limits = c(425, 560)) 

# Increase margins
p$PlantHeight <- p$PlantHeight + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$SLA <- p$SLA + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LA <- p$LA + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p2 <- p2 + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LeafC <- p$LeafC + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$LeafN <- p$LeafN + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LeafP <- p$LeafP + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$Seed.count <- p$Seed.count + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$Seed.lenght <- p$Seed.lenght + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$Seed.mass <- p$Seed.mass + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))

windows()
pall <- grid.arrange(p$PlantHeight, p$SLA, p$LA, p2, p$LeafC,
                     p$LeafN, p$LeafP, p$Seed.count, p$Seed.lenght, p$Seed.mass,
                     layout_matrix = rbind(c(1,2,NA),
                                           c(3,4,NA),
                                           c(5,6,11),
                                           c(7,8,NA),
                                           c(9,10,NA)),
                     mylegend, nrow = 5, ncol = 3, widths = c(2,2,0.5))

ggsave("Figures/Fig6a_GAM_arrival_agri.pdf", pall,
       width = 174, 
       height = 247, units = "mm", dpi = 600)

# draws from posterior ---- 

gam_samples_agri <- function(trait){
  years.since <-  seq(-5000, 5000, by = 500)
  temperature <- seq(0, 12, by = 1)
  pd2 <- data.frame(years.since = years.since,
                    temp = years.since*0)
  Xp2 <- predict(jam_list2[[trait]], newdata = pd2, type = "lpmatrix")
  ii <- 1:25 * 20 + 500

  for (i in ii) {
    fv <- Xp2 %*% sam_list2[[trait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd2$years.since, fv, type = "l", lty = 3, xlab = "Years since arrival of agriculture", 
           main = units$title[units$trait == names(jam_list2[trait])], 
           ylab = paste0(units$units[units$trait == names(jam_list2[trait])], "(log)"),
           xlim = c(-5000, 5000),
           ylim = c(ylim$ymin[ylim$trait == names(jam_list2[trait])],
                    ylim$ymax[ylim$trait == names(jam_list2[trait])]))
    else
      lines(pd2$years.since, fv, lty = 3)
  }
}

png("Figures/SI4-Draws-agri.png", width = 174, 
    height = 300, units = "mm", res = 300)
par(mfrow = c(5,2))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
             "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.lenght", "Seed.mass"), ~gam_samples_agri(.))
dev.off()


# temperature plots ----
plot_func_temp <- function(x,y) {
  plot(x,  shade = TRUE, select = 2,
       scale = 0, xlab = "Temperature (°C)",
       xlim = c(-3, 13),
       main = units$title[units$trait == y])
  abline(h=0)
}

windows()
par(mfrow = c(2,5))
purrr::map2(jam_list2, names(jam_list2), 
            ~plot_func_temp(.x,.y))

gg_func_temp <- function(jam, selectedtrait) {
  CWM <- CWMall %>% 
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
           years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),0)) %>% 
    select(-lag)
  
  values <- plot(jam, select = 2) %>% pluck(2)
  intercept <- coef(jam)[1]
  df <- data.frame(temp = values$x, Mean = values$fit + intercept, 
                   se = values$se)
  
  ggplot(data = CWM, aes(x = temp, y = exp(Mean))) +
    geom_point(color = cbf[1], size = 0.5) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    geom_line(data = df, aes(x = temp, y = exp(Mean)),
              size = 1.5, colour = "#0072B2") +
    geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
                aes(x = temp, y = exp(Mean), ymin = exp(Mean - se), ymax = exp(Mean + se))) +
    scale_y_continuous(units$units[units$trait == selectedtrait]) +
    scale_x_continuous("Temperature (°C)",
                       limits = c(0,13)) +
    theme_bw() +
    theme(legend.position = "none")
}
windows()
q <- purrr::map2(jam_list2, names(jam_list2), 
                 ~gg_func_temp(.x,.y))

# plot LDMC separately (no log transformation)
values <- plot(jam_list2[["LDMC"]], select = 2) %>% pluck(2)
intercept <- coef(jam)[1]
df <- data.frame(temp = values$x, Mean = values$fit + intercept, 
                 se = values$se)

q2 <- ggplot(data = CWM, aes(x = temp, y = Mean)) +
  geom_point(color = cbf[1], size = 0.5) +
  ggtitle(units$title[units$trait == "LDMC"]) +
  geom_line(data = df, aes(x = temp, y = Mean), 
            size = 1.5, colour = "#0072B2") +
  geom_ribbon(data = df, alpha = 0, colour = "#0072B2", fill ="#0072B2",
              aes(x = temp, y = Mean, ymin = Mean - se, ymax = Mean + se)) +
  scale_y_continuous(units$units[units$trait == selectedtrait]) +
  scale_x_continuous("Temperature (°C)",
                     limits = c(0,13)) +
  theme_bw() +
  theme(legend.position = "none")

# set same axis limits
q$LeafC <- q$LeafC +  scale_y_continuous(units$units[units$trait == "LeafC"],limits = c(425, 560)) 

# Increase margins
p$PlantHeight <- p$PlantHeight + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$SLA <- p$SLA + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LA <- p$LA + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p2 <- p2 + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LeafC <- p$LeafC + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$LeafN <- p$LeafN + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$LeafP <- p$LeafP + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$Seed.count <- p$Seed.count + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))
p$Seed.lenght <- p$Seed.lenght + theme(plot.margin = grid::unit(c(0,0.1,0,1),"cm"))
p$Seed.mass <- p$Seed.mass + theme(plot.margin = grid::unit(c(0,0.5,0,0.6),"cm"))

windows()
pall <- grid.arrange(p$PlantHeight, p$SLA, p$LA, p2, p$LeafC,
                     p$LeafN, p$LeafP, p$Seed.count, p$Seed.lenght, p$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     mylegend, nrow = 5, ncol = 3, widths = c(2,2,0.5))
windows()
qall <- grid.arrange(q$PlantHeight, q$SLA, q$LA, q2, q$LeafC, q$LeafN, q$LeafP,
                     q$Seed.count, q$Seed.lenght, q$Seed.mass,
                     layout_matrix = rbind(c(1,2),
                                           c(3,4),
                                           c(5,6),
                                           c(7,8),
                                           c(9,10)),
                     widths = c(4,4),
                     nrow = 5, ncol = 2)

ggsave("Figures/Fig6b_GAM_temp.pdf", qall,  width = 174, 
       height = 247, units = "mm", dpi = 600)
# draws from posterior ----
gam_samples_temp <- function(trait){
  years.since <-  seq(-5000, 5000, by = 500)
  temperature <- seq(2, 12, by = 1)
  pd2 <- data.frame(years.since = temperature*0,
                    temp = temperature)
  Xp2 <- predict(jam_list2[[trait]], newdata = pd2, type = "lpmatrix")
  ii <- 1:25 * 20 + 500
  
  for (i in ii) {
    fv <- Xp2 %*% sam_list2[[trait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd2$temp, fv, type = "l", lty = 3, xlab = "Temperature (°C)", 
           main = units$title[units$trait == names(jam_list2[trait])], 
           ylab = paste0(units$units[units$trait == names(jam_list2[trait])], "(log)"),
           xlim = c(2,12), 
           ylim = c(ylim$ymin[ylim$trait == names(jam_list2[trait])],
                    ylim$ymax[ylim$trait == names(jam_list2[trait])]))
    else
      lines(pd2$temp, fv, lty = 3)
  }
}
png("Figures/SI4-Draws-temp.png", width = 174, height = 300, units = "mm", res = 300)
par(mfrow = c(5,2))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
             "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.lenght", "Seed.mass"), ~gam_samples_temp(.))
dev.off()

