# HVDF_GenerateExclusion.R
# Script to generate the exclusion zone for the plot (wetland)
# and the outer bounds of the plot.
# packages
for(i in c("terra", "sf", "raster", "tidyverse", "magrittr")){
  library(i, character.only = T)
  }

# Locations
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_1 <- "/data/ForestGEO/processed/standardized/"

# bring in functions
source("./scripts/00_Functions.R")


# Generate areas to exclude for the swamp in HVDF
load(paste0(loc_Gdr, loc_1, paste0("SERC", "_AllStemsOLD.r")))
AllStems <- SERC_AllStems

# marking stems incorporated into the fence:
fence_stems.v <- rep(F, nrow(AllStems[[1]]))
for(i in AllStems){
  fence_stems.v <- fence_stems.v | (CodesLogic(i$Codes, c = "F", sep = ";"))
}

# stems with spatial data
stem_sptl <- AllStems[[1]] %>%
  mutate(Fence = fence_stems.v) %>%
  filter(!(is.na(PX) | is.na(PY))) %>%
  st_as_sf(coords = c("PX", "PY"))

# box of area to be excluded (deer exclosure) 
df <- stem_sptl %>%
  filter(Fence) %>%
  summarise(geometry = st_combine( geometry)) %>% # not sure what this does
  st_concave_hull(ratio = .075) %>% # concave hull, ratio chosen by eye
  st_geometry()

# Generate outer box to the plot.
x <- AllStems[[1]]
site_bndry <- as(raster::extent(floor(min(x$PX, na.rm = T)),
                                ceiling(max(x$PX, na.rm = T)),
                                floor(min(x$PY, na.rm = T)),
                                ceiling(max(x$PY, na.rm = T))),
                 "SpatialPolygons") %>% st_as_sf()
# 
### visual confirmation graph
#ggplot() + geom_sf(data = stem_sptl, aes(color = Mnemonic)) + 
#geom_sf(data = site_bndry, aes(alpha = 0)) +
# geom_sf(data = df,aes( alpha = 0))

# --- Write out ---------------------------------------------
# boundaries
st_write(site_bndry, paste0(loc_Gdr, "/data/ForestGEO/spatial/",
                            "SERC_boundaries.shp"),
         delete_layer = T)

# End