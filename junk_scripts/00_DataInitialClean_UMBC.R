# Modifications to specific datasets, before initial processing
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"

# Setting/finding location of the script, so that the WD is automatically within GITHub local copy location.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

# Bringing in the packages:
library("tidyverse")

# Read in the commonly used functions:
source("scripts/00_Functions.R")

# Let's bring in the UMBC data
site <- "UMBC"

# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/", site))[str_which(list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/",site)),
                                                                                  pattern = "_AllStem_")]

# Create a file for storing all of the Stem files
UMBC_raw <- lapply(stem_files, FUN = reader1, loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))


#### Lets work with the data ####################################################################################################################################
# We will start trying to harmonize the data with the other ForestGEO sites. 
# First, thing though, we need to harmonize the two datasets with each other.

# Census 1 has some extra columns that are doing nothing
UMBC_raw[[1]] <- UMBC_raw[[1]][,-c(15:16)]

# Census 2 has quadrat data split into two columns that are redundant (column "quadrat" is just "SITE" + "-" + "num")
UMBC_raw[[2]] <- UMBC_raw[[2]][,-c(4:5)]

# The 'lx/ly' values are different between the two censuses. In census 1, they are the relative position,
# while they seem to represent the position of the quadrat center in census 2
UMBC_raw[[2]] <- UMBC_raw[[2]] %>% mutate(lx = gx - lx, ly = gy - ly)


# Alright, columns are now identical to one another across both censuses.
# Now we need to harmonize UMBC to the other forestGEO sites.




unique(unlist(lapply(UMBC_raw,function(z){z$codes})))




View(UMBC_raw[[1]])

UMBC





head(SERC1)

head(UMBC1)

unique(unlist(str_split((c(UMBC1$codes,UMBC2$codes)),";")))
