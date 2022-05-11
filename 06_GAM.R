# This script produces the GAMs (figure 4, 5a and 5b) and the posterior simulations Appendix S5
# Annegreet Veeken

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
library(furrr)

# memory.limit(size = 999999)

## Prepare data ----
dfCWM <- readRDS("RDS_files/05_multivariate_CWM.rds") 
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
                          "Seed count ",
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
                          "Seed mass")
)

ylim <- tibble(trait = c("PlantHeight","SLA", "LA","LDMC", 
                         "LeafC", "LeafN", "LeafP", 
                         "Seed.count", "Seed.length", "Seed.mass"),
               ymin = c(0.1, 5, 1, 0.2,  380, 10, 1, 1, 1, 0.1),
               ymax = c(20, 30, 4000,0.5, 550, 40, 4, 12000, 15, 500))
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
rm(dfAGRI, dfTEMP, dfTRSH, lCROP, lTRSH, dfCROP, dfCWM)

## GAM's
# GAM Mean ~ s(Time.BP) ----
if(1){
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
  
  saveRDS(jam, file = paste("RDS_files/06_GAM_time_multivariate_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(sam, file = paste("RDS_files/06_GAM_time_sam_multivariate_",
                            selectedtrait, ".rds", sep = ""))
  saveRDS(coda.sam, file = paste("RDS_files/06_GAM_time_coda_multivariate_",
                                 selectedtrait, ".rds", sep = ""))
}
plan(multisession(workers = 2))
furrr::future_map(trait, ~gam_func(.), .options = furrr_options(seed = TRUE))

}

# Plot  GAM's ----
if(1){
# create list of jam's
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_time_multivariate_") %>% 
  str_subset(pattern = "sam", negate = TRUE) %>% 
  str_subset(pattern = "coda", negate = TRUE)
files

traitname <- files %>% 
  str_remove(., paste0("06_GAM_time_multivariate_")) %>% 
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
  str_subset(pattern = "06_GAM_time_sam_multivariate_") 
files

sam_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(sam_list) <- traitname

# create list of coda objects
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_time_coda_multivariate_")
files

coda_list <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(coda_list) <- traitname

# check convergence
pdf("Convergence_diagnostics_GAM1_gelman.pdf")
for(x in 1:10){
gelman.plot(coda_list[[x]], ask = FALSE)
  }
dev.off()

pdf("Convergence_diagnostics_GAM1.pdf")
for(x in 1:10){
plot(coda_list[[x]])
}
dev.off()

# ggplot -----
gg_func <- function(jam, selectedtrait) {
  values <- plot(jam, seWithMean = TRUE) %>% flatten
  intercept <- coef(jam)[1]
  df <- data.frame(age = values$x, Mean = intercept + values$fit, se = values$se) %>% 
    mutate(Min = Mean - se, Max = Mean + se) 
  
  p <- ggplot(data = CWMall[CWMall$trait == selectedtrait,], 
         aes(x = age, y = exp(Mean))) +
    geom_point(aes(colour = tree),  size = .3) +
    scale_color_gradient("Tree pollen\n(%)", low = cbf[5],
                         high = cbf[4], limits = c(0, 100)) +
    scale_y_continuous(units$units[units$trait == selectedtrait], trans = "log10") +
    scale_x_reverse("Time (calibrated years BP)", limits = c(10000,0)) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    geom_ribbon(data = df, alpha = .5,  fill ="black",
                aes(ymin = exp(Mean - se), ymax = exp(Mean + se))) + 
    geom_line(data = df, aes(x = age, y = exp(Mean)), colour = "black",
              size = 0.5) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 10))
  
  df <- df %>% 
    select(age, Mean, Min, Max) %>% 
    mutate(Mean = exp(Mean), Min = exp(Min), Max = exp(Max))
  
  list(values=df,plot=p)
}

windows()
p <- purrr::map2(jam_list, names(jam_list), 
                 ~gg_func(.x,.y))

## plot 1 trait for extracting legend 
selectedtrait <- "LeafP"
traitindex <- which(names(jam_list) %in% selectedtrait)
windows()
values <- plot(jam_list[[traitindex]]) %>% flatten
intercept <- coef(jam_list[[traitindex]])[1]
df <- data.frame(age = values$x, Mean = intercept + values$fit, se = values$se)

p1 <-  ggplot(data = CWMall[CWMall$trait == selectedtrait,], 
              aes(x = age, y = exp(Mean))) +
  geom_point(aes(colour = tree),  size = .3) +
  scale_color_gradient("Tree pollen\n(%)", low = cbf[5],
                       high = cbf[4], limits = c(0, 100)) +
  scale_y_continuous(units$units[units$trait == selectedtrait], trans = "log10") +
  scale_x_reverse("Time (calibrated years BP)", limits = c(10000,0)) +
  ggtitle(units$title[units$trait == selectedtrait]) +
  geom_ribbon(data = df, alpha = .5,  fill ="black",
              aes(ymin = exp(Mean - se), ymax = exp(Mean + se))) + 
  geom_line(data = df, aes(x = age, y = exp(Mean)), colour = "black",
            size = 0.5) +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 10))
mylegend <- get_legend(p1)

windows()
pall <- grid.arrange(p$PlantHeight$plot, p$SLA$plot, p$LA$plot, p$LDMC$plot, p$LeafC$plot,
                     p$LeafN$plot, p$LeafP$plot, p$Seed.count$plot, p$Seed.length$plot, p$Seed.mass$plot,
                    layout_matrix = rbind(c(NA,1,2,NA),
                                          c(NA,3,4,NA),
                                          c(NA,5,6,11),
                                          c(NA,7,8,NA),
                                          c(NA,9,10,NA)),
                    mylegend, nrow = 5, ncol = 4, 
                    widths = c(0.25,2,2,0.5))

ggsave(paste0("Figures/Fig4-GAM_time_multivariate.pdf"), pall, width = 174, 
       height = 300, units = "mm", dpi = 600)
graphics.off()
}

# draws from posteriors for all traits----
if(1){
options(scipen=10)  
gam_samples <- function(selectedtrait){
  pd <- data.frame(age = seq(0, 10000, by = 500))
  Xp <- predict(jam_list[[selectedtrait]], newdata = pd, type = "lpmatrix")
  ii <- 1:25 * 20 + 500
  
for (i in ii) {
    fv <- Xp %*% sam_list[[selectedtrait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd$age, exp(fv), type = "l", lty = 1, xlab = "Time (cal BP)", col = 4,log="y",
           ylab = units$units[units$trait == selectedtrait], xlim = c(10000, 0),
           main = units$title[units$trait == selectedtrait],
           ylim = c(ylim$ymin[ylim$trait == selectedtrait],ylim$ymax[ylim$trait == selectedtrait]),
           yaxt = "n")
    else
      lines(pd$age, exp(fv), lty = 1, col = 4)
      grid(nx = NULL, ny = NULL,
           lty = 1,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)
      axis(side = 2,
           ## Rotate the labels.
           las = 2)
  }
}


# plot
png("Figures/SI5-Draws-time-multivariate.png", width = 174 , height = 300, units = "mm", res = 600)
par(mfrow = c(5,2), mar = c(4.1, 4.1, 4.1, 1.1))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
              "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.length", "Seed.mass"), ~gam_samples(.))
dev.off()
}

# GAM Mean ~ s(since.agri) + s(temp) ----
if(1){
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
           years.since = ifelse(PresenceAgri == "After agriculture", -rev(cumsum(rev(replace_na(lag, 0)))),
                                -cumsum(replace_na(lag, 0)))) %>% 
    filter(years.since > -5000 & years.since < 5000)
  
  # Reference for code:
  # Simon N.Wood (2017). Generalized additive models. An introduction with R (2nd ed). page 374-376
  gam.model <- jagam(Mean ~ s(years.since) + s(temp), data = CWM, 
                     file = "gam_agri0_arrival_temp.jags",
                     diagonalize = TRUE) 
  
  # changes from generated files are saved in gam_agri_time.jags
  # use previously calculated SD on CWM mean
  gam.model$jags.data$sd <- CWM$SD
  
  # add site name and site number for random effect
  gam.model$jags.data$site  <- as.numeric(as.factor(CWM$site.name))
  gam.model$jags.data$nsite <- length(unique(CWM$site.name))
  
  load.module("glm")
  jm <- jags.model("gam_agri_arrival_temp.jags", data = gam.model$jags.data,
                   inits = gam.model$inits, n.chains = 2,  n.adapt = 500)
  sam <- jags.samples(jm, c("b", "rho"),
                      n.iter = 10000, thin = 10)
  saveRDS(sam, file = paste("RDS_files/06_GAM_arrival_multivariate_sam_",
                            selectedtrait, ".rds", sep = ""))
  jam <- sim2jam(sam, gam.model$pregam)
  saveRDS(jam, file = paste("RDS_files/06_GAM_arrival_multivariate_jam_",
                            selectedtrait, ".rds", sep = ""))
  coda.sam <- coda.samples(jm, c("b", "rho"),
                           n.iter = 10000, thin = 10)
  saveRDS(coda.sam, file = paste("RDS_files/06_GAM_arrival_multivariate_coda_",
                                 selectedtrait, ".rds", sep = ""))
}
plan(multisession(workers = 2))
furrr::future_map(trait[c(2:5,7:10)], ~gam_func2(.), .options = furrr_options(seed = TRUE))
}

if(1){
## Plot GAM2 ----
# Compile results to list
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_arrival_") %>% 
  str_subset(pattern = "jam") 
files
traitname <- files %>% 
  str_remove(., paste0("06_GAM_arrival_multivariate_jam_")) %>% 
  str_remove(., ".rds")
traitname
folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

jam_list2 <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(jam_list2) <- traitname

files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "06_GAM_arrival_multivariate_sam") 
files

sam_list2 <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map(~readRDS(.))
names(sam_list2) <- traitname

files <-
  list.files("RDS_files/") %>%
  str_subset(pattern = "06_GAM_arrival_multivariate_coda")
files

coda_list2 <- files %>%
  folderpath.fun(.) %>%
  purrr::map(~readRDS(.))
names(coda_list2) <- traitname

# check convergence
pdf("Convergence_diagnostics_GAM2_gelman.pdf")
for(x in 1:10){
  gelman.plot(coda_list2[[x]], ask = FALSE)
}
dev.off()

pdf("Convergence_diagnostics_GAM2.pdf")
for(x in 1:10){
  plot(coda_list2[[x]])
}
dev.off()

# ggplot -----
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
  values <-  plot(jam, seWithMean = TRUE, select = 1) %>% pluck(1)
  df <- data.frame(years.since = values$x, Mean = values$fit + intercept,
               se = values$se) %>% 
    mutate(Min = Mean - se, Max = Mean + se) 
  
  p <- ggplot(data = CWM, aes(x = years.since, y = exp(Mean))) +
    geom_point(aes(colour = crop), size = 0.3) +
    scale_color_gradient("Crop pollen (%)", low = cbf[5], high = cbf[4],
                         limits = c(0.01,100), trans = "log",breaks = c(0.1,2.5,50),
                         labels = c(0.1,2.5,50), na.value = "#999999") +
    ggtitle(units$title[units$trait == selectedtrait]) +
    scale_y_continuous(units$units[units$trait == selectedtrait], trans = "log10") +
    scale_x_continuous("Years since arrival of agriculture",
                       limits = c(-5000, 5000)) +
    geom_ribbon(data = df, alpha = 0.5, fill ="#0072B2",
                aes(x = years.since, y = exp(Mean), ymin = exp(Mean - se), 
                    ymax = exp(Mean + se))) +
    geom_line(data = df, aes(x = years.since, y = exp(Mean)), size = 0.5, 
              colour = "#0072B2") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 10))
  
  df <- df %>% 
    select(years.since, Mean, Min, Max) %>% 
    mutate(Mean = exp(Mean), Min = exp(Min), Max = exp(Max))
  
  list(values = df, plot = p)
  }

windows()
l <- purrr::map2(jam_list2, names(jam_list2), 
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
values <-  plot(jam_list2[[selectedtrait]], seWithMean = TRUE, select = 1) %>% pluck(1)
df <- data.frame(years.since = values$x, Mean = values$fit + intercept,
                 se = values$se)

l1 <-   ggplot(data = CWM, aes(x = years.since, y = exp(Mean))) +
  geom_point(aes(colour = crop), size = 0.3) +
  scale_color_gradient("Crop pollen\n(%)", low = cbf[5], high = cbf[4],
                       limits = c(0.01,100), trans = "log",breaks = c(0.1,2.5,50),
                       labels = c(0.1,2.5,50), na.value = "#999999") +
  ggtitle(units$title[units$trait == selectedtrait]) +
  scale_y_continuous(units$units[units$trait == selectedtrait], trans = "log10") +
  scale_x_continuous("Years since arrival of agriculture",
                     limits = c(-5000, 5000)) +
  geom_ribbon(data = df, alpha = 0.5, fill ="#0072B2",
              aes(x = years.since, y = exp(Mean), ymin = exp(Mean - se), 
                  ymax = exp(Mean + se))) +
  geom_line(data = df, aes(x = years.since, y = exp(Mean)), size = 0.5, 
            colour = "#0072B2") +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 10))

mylegend <- get_legend(l1)

windows()
lall <- grid.arrange(l$PlantHeight, l$SLA, l$LA, l$LDMC, l$LeafC,
                     l$LeafN, l$LeafP, l$Seed.count, l$Seed.length, l$Seed.mass,
                     layout_matrix = rbind(c(NA,1,2,NA),
                                           c(NA,3,4,NA),
                                           c(NA,5,6,11),
                                           c(NA,7,8,NA),
                                           c(NA,9,10,NA)),
                     mylegend, nrow = 5, ncol = 4, widths = c(0.25,2,2,0.5))

ggsave(paste0("Figures/Fig5a_GAM_arrival_multivariate_agri.pdf"), lall,
       width = 174, 
       height = 247, units = "mm", dpi = 600)
graphics.off()
}

# draws from posterior ---- 
if(1){
options(scipen=10)  
gam_samples_agri <- function(selectedtrait){
  years.since <-  seq(-5000, 5000, by = 500)
  temperature <- seq(0, 12, by = 1)
  pd2 <- data.frame(years.since = years.since,
                    temp = years.since*0)
  Xp2 <- predict(jam_list2[[selectedtrait]], newdata = pd2, type = "lpmatrix")
  ii <- 1:25 * 20 + 500

  for (i in ii) {
    fv <- Xp2 %*% sam_list2[[selectedtrait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd2$years.since, exp(fv), type = "l", lty = 1, xlab = "Years since arrival of agriculture", 
           main = units$title[units$trait == names(jam_list2[selectedtrait])], 
           ylab = paste0(units$units[units$trait == names(jam_list2[selectedtrait])]), col =4,
           xlim = c(-5000, 5000), log = "y",
           ylim = c(ylim$ymin[ylim$trait == names(jam_list2[selectedtrait])],
                    ylim$ymax[ylim$trait == names(jam_list2[selectedtrait])]),
           yaxt = "n"
            )
    else
      lines(pd2$years.since, exp(fv), lty = 1, col = 4)
      grid(nx = NULL, ny = NULL,
           lty = 1,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)
      grid(nx = NULL, ny = NULL,
           lty = 1,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)
      axis(side = 2,
           ## Rotate the labels.
           las = 2)
  }
}


# plot
png("Figures/SI5-Draws-agri-multivariate.png", width = 174, 
    height = 300, units = "mm", res = 300)

par(mfrow = c(5,2), mar = c(4.1, 4.1, 4.1, 1.1))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
             "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.length", "Seed.mass"), ~gam_samples_agri(.))
dev.off()
}
if(1){
# temperature plots ----
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
  
  values <- plot(jam, select = 2, seWithMean = TRUE) %>% pluck(2)
  intercept <- coef(jam)[1]
  df <- data.frame(temp = values$x, Mean = values$fit + intercept, 
                   se = values$se)
  
  ggplot(data = CWM, aes(x = temp, y = exp(Mean))) +
    geom_point(color = cbf[1], size = 0.5) +
    ggtitle(units$title[units$trait == selectedtrait]) +
    geom_line(data = df, aes(x = temp, y = exp(Mean)),
              size = 0.5, colour = "#0072B2") +
    geom_ribbon(data = df, alpha = 0.5, fill = "#0072B2",
                aes(ymin = exp(Mean - se), ymax = exp(Mean + se))) + 
    scale_y_continuous(units$units[units$trait == selectedtrait], trans = "log10") +
    scale_x_continuous("Temperature (°C)",
                       limits = c(0,13)) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 10))
}
windows()
q <- purrr::map2(jam_list2, names(jam_list2), 
                 ~gg_func_temp(.x,.y))

windows()
qall <- grid.arrange(q$PlantHeight, q$SLA, q$LA, q$LDMC, q$LeafC, q$LeafN, q$LeafP,
                     q$Seed.count, q$Seed.length, q$Seed.mass,
                     layout_matrix = rbind(c(NA,1,2),
                                           c(NA,3,4),
                                           c(NA,5,6),
                                           c(NA,7,8),
                                           c(NA,9,10)),
                     widths = c(0.25,2,2),
                     nrow = 5, ncol = 3)

ggsave(paste0("Figures/Fig5b_GAM__multivariate_temp.pdf"), qall,  width = 155, 
       height = 247, units = "mm", dpi = 600)
graphics.off()
}

if(1){
# draws from posterior ----
options(scipen=10)  
gam_samples_temp <- function(selectedtrait){
  years.since <-  seq(-5000, 5000, by = 500)
  temperature <- seq(2, 12, by = 1)
  pd2 <- data.frame(years.since = temperature*0,
                    temp = temperature)
  Xp2 <- predict(jam_list2[[selectedtrait]], newdata = pd2, type = "lpmatrix")
  ii <- 1:25 * 20 + 500
  
  for (i in ii) {
    fv <- Xp2 %*% sam_list2[[selectedtrait]]$b[, i, 1]
    if (i == ii[1])
      plot(pd2$temp, exp(fv), type = "l", lty = 1, col = 4, xlab = "Temperature (°C)", log = "y",
           main = units$title[units$trait == names(jam_list2[selectedtrait])], 
           ylab = paste0(units$units[units$trait == names(jam_list2[selectedtrait])]),
           xlim = c(2,12), 
           ylim = c(ylim$ymin[ylim$trait == names(jam_list2[selectedtrait])],
                    ylim$ymax[ylim$trait == names(jam_list2[selectedtrait])]),
           yaxt= "n")
    else
      lines(pd2$temp, exp(fv), lty = 1, col = 4)
      grid(nx = NULL, ny = NULL,
           lty = 1,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)      
      axis(side = 2,
           ## Rotate the labels.
           las = 2)
  }
}

png("Figures/SI5-Draws-temp-multivariate.png", width = 174, height = 300, units = "mm", res = 300)
par(mfrow = c(5,2), mar = c(4.1, 4.1, 4.1, 1.1))
purrr::map(c("PlantHeight","SLA", "LA","LDMC", 
             "LeafC", "LeafN", "LeafP", 
             "Seed.count", "Seed.length", "Seed.mass"), ~gam_samples_temp(.))
dev.off()
}

graphics.off()
