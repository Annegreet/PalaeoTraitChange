library(tidyverse)

# read in trait data, gapfiled 2020 version
gapfilled <- read.csv("Data/mean_gap_filled_back_transformed.csv")
spec <- read.csv("Data/hierarchy.info.csv")
traitnames <- read.delim("Data/TryAccTraits.txt")

# add species column
dfTRAITS <- bind_cols(spec, gapfilled)

# rename column names
dfTRAITS <- dfTRAITS %>% 
  rename(SSD = X4,
         Root.depth = X6,
         SLA = X11,
         LeafC = X13,
         LeafN = X14,
         LeafP = X15,
         X18 = X18, # traitID 18 undefined in try and readme... 
         Stem.diam = X21,
         Seed.mass = X26,
         Seed.lenght = X27,
         Leaf.thick = X46,
         LDMC = X47,
         LeafN.LA = X50,
         LDM = X55,
         deltaN15 = X78,
         Seed.germ = X95,
         Seed.count = X138,
         Wood.ves.length = X282,
         Wood.fib.length = X289,
         Root.length.mass = X1080,
         LA = X3112,
         LA.leaflet = X3113,
         LA.undif = X3114,
         Leaf.water = X3120,
         Leaf.length = X144,
         Leaf.width = X145,
         LeafCN = X146,
         Leaf.fresh.mass = X163,
         Stem.conduit.dens = X169,
         Chromosome.number = X223,
         cDNA.content = X224,
         Disp.unit.length = X237,
         Stem.conduit.diam = X281) %>% 
  select(-X,-X1, -X18)

head(dfTRAITS)
saveRDS(dfTRAITS, "RDS_files/Gap_filled_traits_2020.rds")


