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


# Generate areas to exclude for the swamp in HVDF
load(paste0(loc_Gdr, loc_1, paste0("HVDF", "_AllStemsOLD.r")))
AllStems <- HVDF_AllStems

stem_sptl <- AllStems[[1]] %>%
  filter(!(is.na(PX) | is.na(PY))) %>%
  st_as_sf(coords = c("PX", "PY"))

# Stems not sampled during the second 
swmp_stems <- AllStems[[2]]$StemID[AllStems[[2]]$DFstatus == "stem_gone"]

# draw a concave shape around these stems, the ratio was 
# found through trial and error.
df <- stem_sptl %>%
  filter(StemID %in% swmp_stems) %>%
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

## visual confirmation graph
# ggplot() + geom_sf(data = stem_sptl, aes(color = Mnemonic)) + 
# geom_sf(data = site_bndry, aes(alpha = 0))

# --- Write out ---------------------------------------------
# exlcusion
st_write(df, paste0("G:/Shared drives/NASA_ForestScaling/",
                    "LifeHistoryStrats/data/ForestGEO/",
                    "spatial/HVDF_exclusion.shp"),
         delete_layer = T)
# boundaries
st_write(site_bndry, paste0(loc_Gdr, "/data/ForestGEO/spatial/",
                            "HVDF_boundaries.shp"),
         delete_layer = T)

# End