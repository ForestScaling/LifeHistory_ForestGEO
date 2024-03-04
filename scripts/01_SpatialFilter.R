# ### Setup ###################################################################
# Packages
for(i in c("terra","sf","raster", "tidyverse", "magrittr")){
  library(i, character.only = TRUE)
}

# File locations - will work to stdize loc with a sep script
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_1 <- "/data/ForestGEO/processed/standardized/"

# set working directory to the current location.
# I think(?) this is currently dependent on the RStudio IDE, fix this later
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")


# ### Begin individual sites ##################################################
site.v <- c("SCBI", "HVDF", "SERC")

for(site in site.v){

# read in the site stem data
load(paste0(loc_Gdr, loc_1, paste0(site, "_AllStems.r")))

# --- Reading in the boundaries -----------------------------------------------
sptl_file <- paste0(loc_Gdr, "/data/ForestGEO/spatial/",
                    site, "_boundaries")
# Check that all of the components of the shp file exist
if(all(file.exists(paste0(sptl_file, c(".shp", ".shx", ".dbf"))))){
  # Pass
  sptl_bounds <- read_sf(paste0(sptl_file, ".shp"))
}else{
  # fail
  stop(paste("For", site, "boundaries '.shp' filis incomplete."))
}

# --- Reading in the exclusion areas ------------------------------------------
sptl_file <- paste0(loc_Gdr, "/data/ForestGEO/spatial/",
                    site, "_exclusion")
# Check that all of the components of the shp file exist
excl.log <- all(file.exists(paste0(sptl_file, c(".shp", ".shx", ".dbf"))))

# read in the files.
if(excl.log){
  # Pass
  sptl_excl <- read_sf(paste0(sptl_file, ".shp"))
  cat("Successfully read in exclusion-area for: ", site, "\n", sep = "")
}else{
  # fail
  cat(paste(site, "has no exclusion-area files."))
}

# dataset of the stems 
stem_sptl <- AllStems[[1]] %>%
  filter(!is.na(PX) & !is.na(PY)) %>%
  st_as_sf(., coords = c("PX","PY")) %>% 
  st_set_crs(., st_crs(sptl_bounds))

# ### Assign stems to the grid-pattern based on ###############################
subp_dim <- 25 # size of the subplots m x m
cat("Using subplots of ", subp_dim, "m x ", subp_dim, "m, please ensure that ",
    "the coordinate system is in meters.\n", sep = "")


# generate a grid of the subplots.
grd <- as(raster::extent(sptl_bounds),
          "SpatialPolygons") %>%
  st_as_sf() %>%
  st_make_grid(., cellsize = c(subp_dim, subp_dim)) %>%
  st_set_crs(., st_crs(sptl_bounds)) # keeps the same CRS as the bounds file


# get rid of the edge subplots
grd <- grd[st_contains(sptl_bounds, grd, sparse = F)]


# remove suplots that are in the excluded zone
if(excl.log){
  grd <- grd[!st_intersects(grd, sptl_excl, sparse = F)]
}


# subset stems that are contained in valid grid squares.
stem_kp.v <- (st_touches(grd, stem_sptl, # touches edge of subplt
                     sparse = F) | 
            st_contains(grd, stem_sptl, # w/in a subplt
                        sparse = F)) %>%
  # select the subplt for the stem
  apply(., MARGIN = 2, function(x){
    if(any(x)){ # found in at least one subplt
      return(max(which(x))) # max is to resolve cases on the boundary,
    }else{
      return(NA) # ones outside get an NA
    }
  })


# --- Generate plot for viewing areas excluded ----------------------------
# add in the data
stem_sptl <- stem_sptl %>%
  mutate(subplt = stem_kp.v) %>%
  st_set_crs(NA)



# Plot
plt <- ggplot() +
  geom_sf(data = stem_sptl, aes(color = is.na(subplt))) + 
  scale_color_manual(values = c("#40798C", "#AF3E4D")) # set colors
# add in the grid subplots and the excluded areas
plt <- plt +
  geom_sf(data = st_set_crs(grd, NA), alpha = 0) # add grid for subplot
# add excluded If applicable
if(excl.log){
  plt <- plt +
    geom_sf(data = st_set_crs(sptl_excl, NA),
            alpha = 0.5)
}


# formatting
plt <- plt +
  guides(color = guide_legend(title = "Excluded")) + # legend
  ggtitle(label = paste("Stem plot:", site)) + #title
  theme_bw() # theme

print(plt) # show the plot

# --- Now we go through and add the sub-plot for each stem ----------------
AllStems <- lapply(AllStems, FUN = function(X){
  X %>% mutate(AlgoSubplot = stem_kp.v) %>% return()
})




# HERE!!




# write out the file
save(AllStems,
     file = paste0(loc_Gdr, loc_1, paste0(site, "_AllStems.r")))
}


# Detach the spatial packages, the 'select' function isn't playing nicely
for(i in c("sf", "terra", "raster")){
  detach(paste0("package:", i), character.only = T)
}

# End



# 
# 
# 
# 
# # Read in the dataset
# x <- Raw_Stems[[1]]
# 
# 
# # Shape file denoting the outer edge of the plot
# site_bndry <- as(raster::extent(min(x$PX, na.rm = T),
#                                 max(x$PX, na.rm = T),
#                                 min(x$PY, na.rm = T),
#                                 max(x$PY, na.rm = T)),
#                  "SpatialPolygons") %>% st_as_sf()
# 
# # grid of canopy layers
# grd <- as(raster::extent(min(x$PX, na.rm = T),
#                             max(x$PX, na.rm = T),
#                             min(x$PY, na.rm = T),
#                             max(x$PY, na.rm = T)),
#              "SpatialPolygons") %>%
#   st_as_sf() %>%
#   st_make_grid(., cellsize = c(l, l))
# 
# # convert stem dataframe into spatial, does not handle NA values in position
# stem_sptl <- Raw_Stems[[1]] %>%
#   filter(!is.na(PX) & !is.na(PY)) %>%
#   st_as_sf(., coords = c("PX","PY"))
# 
# # grid square subset. grid must be totally enclosed in plot.
# grd <- grd[st_contains(site_bndry, grd, sparse = F)]
# 
# # subset stems that are contained in valid grid squares.
# asgns <- (st_touches(grd, stem_sptl, # touches edge of subplt
#                      sparse = F) | 
#             st_contains(grd, stem_sptl, # w/in a subplt
#                         sparse = F)) %>%
#   # select the subplt for the stem
#   apply(., MARGIN = 2, function(x){
#     if(any(x)){ # found in at least one subplt
#       return(max(which(x))) # max is to resolve cases on the boundary,
#       }else{
#         return(NA) # ones outside get an NA
#         }
#     })
# 
# # ### Exclusion Zones ###############################################
# # Dealing with areas WITHIN the boundaries that should be excluded.
# # Read in the file
# 
# 
# 
# # we will need to have these filtered out!!!!!
# stem_sptl %>%
#   mutate(subplt = asgns) %>%
#   ggplot(aes(color = subplt)) + geom_sf()
# 
# 
# 
# 
# 
# out <- st_contains( grd, site_bndry)
# plot(out)
# 
# 
# pts <- st_sfc(st_point(c(.5,.5)), st_point(c(1.5, 1.5)), st_point(c(2.5, 2.5)))
# pol = st_polygon(list(rbind(c(0,0), c(2,0), c(2,2), c(0,2), c(0,0))))
# (lst = st_intersects(pts, pol))
# (mat = st_intersects(pts, pol, sparse = FALSE))
# 
