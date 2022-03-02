# Load libraries
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggmap)
library(gridExtra)

## Load and prepare data -----
dfCWM <- readRDS("RDS_files/05_multivariate_CWM_100.rds")

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

pca_hull <- 
  pca_points %>% 
  group_by(PresenceAgri) %>% 
  slice(chull(PC1, PC2))

pca_load <- 
  as_tibble(pca$rotation, rownames = 'variable') %>% 
  mutate(variable = dplyr::recode(variable,
                                  "SLA" = "Specific leaf area", 
                                  "PlantHeight" ="Plant height",
                                  "LDMC" = "Leaf dry matter content",
                                  "LA" = "Leaf area", 
                                  "LeafC" = "Leaf carbon", 
                                  "LeafN" =  "Leaf nitrogen",  
                                  "LeafP" =  "Leaf phosphorus", 
                                  "Seed.count" = "Seed count", 
                                  "Seed.length" = "Seed length",
                                  "Seed.mass" = "Seed mass"))

## PC1 & PC2 plot----
l <- ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = PresenceAgri), alpha = 0.5) +
  # geom_polygon(data = pca_hull, aes(fill = PresenceAgri, colour = PresenceAgri),
  #              alpha = 0.1,
  #              show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*9,
                   yend = PC2*9),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = c(-6.4,  1.90, 4, -4.3, 5, 4, -6, -5.5, 3.4,  0.91),
                   y = c(-0.40, -4.5, 2.15, -0.8, 1, 4.5, 1.26, 1.60, -3.8, 4.9),
           label = pca_load$variable,
           size = 4) +
  scale_x_continuous(name = paste0("PC1 (", round(100 * pca_sum$importance[2,"PC1"],2), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(100 * pca_sum$importance[2,"PC2"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.margin=grid::unit(c(1,1,1,1), "cm")) +
  scale_color_manual(values = c("#999999","#E69F00"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 

windows(6.84,6.84)
l
ggsave("Figures/Fig2-PCA-multivariate-v2.pdf", l, dpi = 600, width = 174,
       height = 174, units = "mm")
ggsave("Figures/Fig2-PCA-multivariate.png", l, dpi = 600, width = 174,
       height = 174, units = "mm")

## PC2 & PC3 plot----
pca_hull23 <- 
  pca_points %>% 
  group_by(PresenceAgri) %>% 
  slice(chull(PC2, PC3))

g <- ggplot(pca_points, aes(x = PC2, y = PC3)) +
  geom_point(aes(colour = PresenceAgri)) +
  geom_polygon(data = pca_hull23, aes(fill = PresenceAgri, colour = PresenceAgri),
               alpha = 0.1,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC2*10,
                   yend = PC3*10),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = c(-0.60,-4.6,2.15,-0.69,2.75,4.82,2.5,2.75,-5.1,4.8), 
           y = c(1.61,4,5.14,4.76,4.03,2.08,2.3,3.34,2.68,-0.6),
           label = pca_load$variable,
           size = 4) +
  scale_x_continuous(name = paste0("PC2 (", round(100 * pca_sum$importance[2,"PC2"],2), "%)")) +
  scale_y_continuous(name = paste0("PC3 (", round(100 * pca_sum$importance[2,"PC3"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("#999999","#E69F00"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 
g
ggsave("Figures/SI8-PCA23_agri-multivariate.png", g, dpi = 600, width = 174, 
       height = 174, units = "mm")

## PC1 & PC3 plot----
pca_hull13 <- 
  pca_points %>% 
  group_by(PresenceAgri) %>% 
  slice(chull(PC1, PC3))

p <- ggplot(pca_points, aes(x = PC1, y = PC3)) +
  geom_point(aes(colour = PresenceAgri)) +
  geom_polygon(data = pca_hull13, aes(fill = PresenceAgri, colour = PresenceAgri),
               alpha = 0.1,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*10,
                   yend = PC3*10),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x =  c(-6.5, 1.90,  2.60, -3.55, 5.2, 3.9, -5.9, -5.5,  4.4,  0.91) , 
           y = c(1.61, 4.1, 5.14, 4.76, 4.5, 2.08, 2.21, 3.34, 2.68, -0.60),
           label = pca_load$variable,
           size = 4) +
  scale_x_continuous(name = paste0("PC1 (", round(100 * pca_sum$importance[2,"PC1"],2), "%)")) +
  scale_y_continuous(name = paste0("PC3 (", round(100 * pca_sum$importance[2,"PC3"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("#999999","#E69F00"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 
p

ggsave("Figures/SI8-PCA13_agri-multivariate.png", p, dpi = 600, width = 174, 
       height = 174, units = "mm")

## PCA map ----
sites <- readRDS("RDS_files/03_PollenWEU-Harmonised.rds") %>% 
  purrr::map(., ~dplyr::select(., site.name, lat, long) %>% 
               unique) %>% 
  bind_rows()

timecat <- c(-200, seq(from = 1000, to = 8000, by = 1000),10000, Inf)
timecat.lab <- c("0-1000 BP", "1000-2000 BP", "2000-3000 BP", "3000-4000 BP", "4000-5000 BP",
                 "5000-6000 BP", "6000-7000 BP", "7000-8000 BP", "8000-10000 BP")


df <- left_join(pca_points, sites, by = c("site.name" = "site.name")) %>% 
  mutate(., time.bin = cut(as.numeric(age), breaks = timecat,right = TRUE,
                           labels = c(timecat.lab,NA))) %>% 
  select(PC1, site.name:time.bin)

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
  coord_map("sinusoidal")

(pcamap <- eumap + 
    geom_point(data = df[df$PresenceAgri == "After agriculture",],
               aes(x = long, y = lat), pch = 21, fill = NA, colour = "black", size = 2) +
    geom_point(data = df, aes(x = long, y = lat, color = PC1), size = 1.5) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    ggtitle("") +
    scale_color_viridis_c() +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    facet_wrap(~time.bin, ncol = 3))

ggsave("Figures/Fig3-PCA-map-multivariate.pdf", pcamap,  
       width = 174, 
       height = 174, units = "mm", dpi = 600)

## PCA scores table ----
compscores <- 
  as_tibble(pca$rotation, rownames = "variable") %>% 
  mutate(variable = dplyr::recode(variable,
                                  "SLA" = "Specific leaf area", 
                                  "PlantHeight" ="Plant height",
                                  "LDMC" = "Leaf dry matter content",
                                  "LA" = "Leaf area", 
                                  "LeafC" = "Leaf carbon", 
                                  "LeafN" =  "Leaf nitrogen",  
                                  "LeafP" =  "Leaf phosphorus", 
                                  "Seed.count" = "Seed count", 
                                  "Seed.lenght" = "Seed length",
                                  "Seed.mass" = "Seed mass")) %>% 
  select(variable, PC1,PC2,PC3) %>% 
  arrange(PC1,PC2,PC3) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 3)))

# add proportion explained variation
sum_pca <- summary(pca) %>% 
  pluck(.,"importance") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "variable") %>% 
  filter(variable == "Proportion of Variance") %>% 
  dplyr::select(variable:PC3) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 3)*100))

compscores  <- compscores %>% 
  bind_rows(sum_pca) 

write.csv(compscores, "Output/Sample_scores_CWM_PCA-multivariate.csv")

