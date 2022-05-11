# This script produces the PCA plots and map fig. 2, fig.3 and 
# supplementary S6

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggmap)
library(gridExtra)

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

cbf <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## PC1 & PC2 plot----
l <- ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = PresenceAgri)) +
  geom_polygon(data = pca_hull, aes(fill = PresenceAgri, colour = PresenceAgri),
               alpha = 0.1,
               show.legend = FALSE) +
  geom_segment(data = pca_load,
               aes(x = 0, y = 0, 
                   xend = PC1*9,
                   yend = PC2*9),
               size = .5,
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = c(-0.88,-4,-3.5,-4.2,-1.5,
                         2,-3.5,-2,-4,2.5),
                   y = c(5,-2,-3,0.5,-4.3,
                         -4.5,2.8,3.3,-0.9,-0.4),
           label = pca_load$variable,
           fontface = 2,
           size = 4) +
  scale_x_continuous(name = paste0("PC1 (", round(100 * pca_sum$importance[2,"PC1"],2), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(100 * pca_sum$importance[2,"PC2"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.margin=grid::unit(c(1,1,1,1), "cm")) +
  scale_color_manual(values = c("#009E73","#F0E442"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 

windows(6.84,6.84)
l
ggsave("Figures/Fig2-PCA-multivariate.pdf", l, dpi = 600, width = 174,
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
                   xend = PC2*9,
                   yend = PC3*9),
               size = .5,
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = c(6.3, -3, -2.8, 0.15, -4, 
                         -5.8, 4, 4, -3, -0.4), 
           y = c(1.3, 0.5, -1.2, 0.5, 3.1,
                 -0.8, 2, 0.2, 0, 8.5),
           label = pca_load$variable,
           fontface = 2,
           size = 4) +
  scale_x_continuous(name = paste0("PC2 (", round(100 * pca_sum$importance[2,"PC2"],2), "%)")) +
  scale_y_continuous(name = paste0("PC3 (", round(100 * pca_sum$importance[2,"PC3"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("#009E73","#F0E442"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 
windows(6.84,6.84)
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
               size = .5,
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x =  c(0.2, -4.3, -4, -4.4, -1.16,
                          2.5, -3.8, -3.5, -4, 1.39) , 
                   y = c(1.4, 0.75, -1.25, 0.25, 3.5,
                         -1.2, 2, 1.1, -0.3, 9.5),
           fontface = 2,
           label = pca_load$variable,
           size = 4) +
  scale_x_continuous(name = paste0("PC1 (", round(100 * pca_sum$importance[2,"PC1"],2), "%)")) +
  scale_y_continuous(name = paste0("PC3 (", round(100 * pca_sum$importance[2,"PC3"],2), "%)")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_color_manual(values = c("#009E73","#F0E442"), 
                     aesthetics = c("colour", "fill"),
                     name = "") 
windows(6.84,6.84)
p

ggsave("Figures/SI8-PCA13_agri-multivariate.png", p, dpi = 600, width = 174, 
       height = 174, units = "mm")

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

