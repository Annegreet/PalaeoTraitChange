# This script produces the figure 2
# Annegreet Veeken

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggmap)
library(gridExtra)
library(sp)
library(rgdal)
library(gratia)

## Load and prepare data -----
dfCWM <- readRDS("RDS_files/05_multivariate_CWM.rds")

dfAGRI <- readRDS("RDS_files/Archaeological_indicators.rds") %>% 
  filter(age <= 10000) %>% 
  select(site.name, age, PresenceAgri) %>% 
  distinct()

dfCWM <- dfCWM %>% 
  left_join(., dfAGRI, by = c("site.name", "age")) %>%
  dplyr::select(site.name, age, trait, Mean, SD, PresenceAgri) %>% 
  filter(age <= 10000)

# re-order levels
dfCWM$PresenceAgri <- factor(dfCWM$PresenceAgri, levels = c("Before agriculture", 
                                                            "After agriculture"))

traitlab <-  c("SLA" = "Specific leaf area", "PH" ="Plant height","LDMC" = "Leaf dry matter content",
               "LA" = "Leaf area", "LeafC" = "Leaf carbon", "LeafN" =  "Leaf nitrogen", 
               "LeafP" =  "Leaf phosphorus", "Seed.count" = "Seed count", "Seed.lenght" = "Seed length",
               "Seed.mass" = "Seed mass")

sites <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds") %>% 
  purrr::map(., ~dplyr::select(., site.name, lat, long) %>% 
               unique) %>% 
  bind_rows()

## PCA ----
dfPCA <- dfCWM %>% 
  dplyr::select(-SD) %>% 
  spread(trait, Mean)

pca <- dfPCA %>% 
  dplyr::select(-site.name, -age, -PresenceAgri) %>%
  prcomp(., scale. = TRUE) 

pca_sum <- summary(pca)

pca_points <- 
  #  convert the pca results to a tibble
  as_tibble(pca$x) %>% 
  bind_cols(dfCWM[dfCWM$trait == "LA", c("site.name","age", "PresenceAgri")])

cbf <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Make clusters based on proximity of sites ----
x <- sites$long
y <- sites$lat

xy <- SpatialPointsDataFrame(matrix(c(x,y), ncol=2), data.frame(ID=seq(1:length(x))),
                             proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

xy <- spTransform(xy, CRS("+init=epsg:27700 +datum=WGS84"))
chc <- hclust(dist(data.frame(rownames=rownames(xy@data), x=coordinates(xy)[,1],
                              y=coordinates(xy)[,2])), method="complete")

# Make 6 groups  
chc.cut <- cutree(chc, k=6) 

# Join clust results to pca results
sites_clust <- sites %>% 
  mutate(., cluster = paste0("Clust",as.factor(chc.cut)))
pca_points_clust <- pca_points %>% 
  left_join(sites_clust) %>% 
  mutate(site.name = as.factor(site.name))

## Make plot ----
if(1){
  # get base map
  eu <- ggplot2::map_data("world", 
                          region = c("UK", "Ireland","Netherlands",
                                     "Germany", "France", "Denmark",
                                     "Sweden", "Switzerland",
                                     "Belgium", "Czech Republic",
                                     "Austria", "Luxembourg"))
  eumap <- ggplot() + 
    geom_polygon(data = eu, aes(x = long, y = lat, group = group),
                 colour = "#999999", fill = "#999999", 
                 show.legend = FALSE) +
    coord_map("bonne", lat0 = 45)
  
  # Create labeled colors
  clus_col <- cbf[2:7]
  clus_lab <- paste0("Clust" ,1:6)
  names(clus_col) <- clus_lab
  
  # plot clusters on map
  clustmap <- eumap +
    geom_point(data = sites_clust, aes(x = long, y = lat,
                                       colour = cluster), size = 1) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    scale_color_manual(values = clus_col) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) 
  clustmap
  
  # run GAM on PC1 and PC2, plot
  plot_pca1 <- function(clus){
    dat <- pca_points_clust %>% 
      filter(cluster == clus) # filter for cluster
    # run gam
    mod <- gam(PC1 ~ s(age) + s(site.name, bs="re"), data = dat, method = "REML")
    
    # check gam
    png(file = paste0("Diagnostics/GAM_check_PC1_", clus,".png"))
    par(mfrow = c(2,2))
    gam.check(mod)
    dev.off()
    
    # extract values
    values <- plot(mod, residuals = TRUE, select = 1, seWithMean = TRUE, shift = coef(mod)[1]) %>% pluck(1)
    fit <- data.frame(age = values$x, pc = values$fit,
                      se = values$se)
    resid <- data.frame(age = values$raw, pc = values$p.resid, PresenceAgri = dat$PresenceAgri)
    
    # plot
    p <- ggplot() +
      geom_point(data = resid, aes(x = age, y = pc), colour = clus_col[clus], size = 0.3) +
      ggtitle("PCA1") +
      scale_y_continuous("", limits = c(-5,9)) +
      scale_x_reverse("Time (BP)", breaks = c(10000, 0), limits = c(10000, 0)) +
      geom_ribbon(data = fit, alpha = 0.5, fill = "black",
                  aes(x = age, y = pc, ymin = pc - se, 
                      ymax = pc + se)) +
      geom_line(data = fit, aes(x = age, y = pc), size = 0.5, 
                colour = "black") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 8),
            panel.grid = element_blank(),
            plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
    p
  }
  # run for all clusters
  pca1 <- map(clus_lab, ~plot_pca1(.))
  names(pca1) <- clus_lab
  
  # PC2
  plot_pca2 <- function(clus){
    # run gam
    dat <- pca_points_clust %>% 
      filter(cluster == clus)
    mod <- gam(PC2 ~ s(age) + s(site.name, bs="re"), data = dat, method = "REML")
    
    # check gam
    png(file = paste0("Diagnostics/GAM_check_PC2_", clus,".png"))
    par(mfrow = c(2,2))
    gam.check(mod)
    dev.off()
    
    # extract values
    values <- plot(mod, residuals = TRUE, select = 1, seWithMean = TRUE, shift = coef(mod)[1]) %>% pluck(1)
    fit <- data.frame(age = values$x, pc = values$fit,
                      se = values$se)
    resid <- data.frame(age = values$raw, pc = values$p.resid, 
                        PresenceAgri = dat$PresenceAgri)
    
    # plot
    p <- ggplot() +
      geom_point(data = resid, aes(x = age, y = pc), colour = clus_col[clus], size = 0.3) +
      ggtitle("PCA2") +
      scale_y_continuous("", limits = c(-9,9)) +
      scale_x_reverse("Time (BP)",breaks = c(10000, 0), limits = c(10000, 0)) +
      geom_ribbon(data = fit, alpha = 0.5, fill = "black",
                  aes(x = age, y = pc, ymin = pc - se, 
                      ymax = pc + se)) +
      geom_line(data = fit, aes(x = age, y = pc), size = 0.5, 
                colour = "black") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 8),
            panel.grid = element_blank(),
            plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
    p
  }
  # Run for all clusters
  pca2 <- map(clus_lab, ~plot_pca2(.))
  names(pca2) <- clus_lab
  
  # Combine plots
  windows(width = 7,  height = 7)
  pcamap <- grid.arrange(pca1$Clust1, pca1$Clust2, pca1$Clust3, pca1$Clust4,
                         pca1$Clust5, pca1$Clust6, clustmap, pca2$Clust1,
                         pca2$Clust2, pca2$Clust3, pca2$Clust4, pca2$Clust5,
                         pca2$Clust6,
                         layout_matrix = rbind(rep(NA, 11),
                                               c(NA, 4, NA, 11, NA, NA, NA, NA, NA, NA, NA),
                                               rep(NA, 11),
                                               c(NA, 2, NA, 9, 7, 7,7, 3, NA, 10, NA),
                                               c(NA, NA, NA, NA, 7, 7, 7, NA, NA, NA, NA),
                                               c(NA, NA, NA, NA, 7, 7, 7, 6, NA, 13, NA),
                                               rep(NA, 11),
                                               c(NA, 1, NA, 8, NA, 5, NA, 12, NA, NA, NA),
                                               rep(NA, 11)),
                         widths = c(0.25, 1, 0.05, 1, 1, 1, 0.05, 1, 0.05, 1, 0.25),
                         heights = c(0.25, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.25))
  # Save
  ggsave("Figures/Fig3-PCA-map-clust.pdf", pcamap,  
         width = 200,  height = 175, units = "mm", dpi = 600)
}
graphics.off()