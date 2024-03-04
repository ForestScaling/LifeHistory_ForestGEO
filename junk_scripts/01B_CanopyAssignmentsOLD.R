# File 01B_CanopyAssignments.R
# Author: JMR
# Last Update: 12023-01-19
# Description: Takes in Standardized/harmonized raw-stem data and outputs the
# canopy layer assignments 
# ### Setup ###################################################################
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_1 <- "/data/ForestGEO/processed/standardized/"

#'
#+ echo = FALSE
# For compiling a report, I don't want to display the whole working
# directory nonsense.

# Bringing in the packages:
packages <- c("tidyverse", "magrittr", "lme4", "rlist", "rstudioapi")
for(p in packages){library(p,character.only = T)}

# Setting/finding location of the script, so that the WD is automatically
# within GITHub local copy location.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")



#setwd("C:/Users/juan.m.rodriguez/Documents/ProgramFiles/Github/NASA_ForestScaling/LifeHistory/forestGEO")



# Get the commonly used functions
source("scripts/00_Functions.R")

#Site
site.v <- c("HVDF","SCBI","SERC")

# --- Looping through each of the sites. --------------------------------------
for(site in site.v){
# Find all the Stem files from this site
load(paste0(loc_Gdr, loc_1, paste0(site, "_AllStems.r")))
AllStems <- get(paste0(site, "_AllStems"))

# Read in the taxa table for this site
load(paste0(loc_Gdr, "/data/ForestGEO/taxa/", site, "_TaxaTable.r"))
Spp_table <- get(paste0(site, "_TaxaTable"))

# ### Begin Manipulations #####################################################
# We are going to add in some columns, however to keep the df from getting too
# bulky I am going to remove some of them before sending out the data. This is
# to keep track of which ones are worth keeping from the original data.
#
cn2k <- colnames(AllStems[[1]]) 

# Add in the taxonomic information:
AllStems <- lapply(AllStems, function(x){
  
  x %<>% left_join(., select(Spp_table,
                            -c(SubSpecies, Authority, No., Species,
                               PriorNames, SpeciesID, Latin)),
                  by = "Mnemonic")
  
  return(x)
})

# Start off with dropping taxa that are not to be used for canopy level
# estimates (i.e. lianas)
dropme <- Spp_table$Mnemonic[which(Spp_table$DropBegin == 1)]
for(i in 1:length(AllStems)){
  AllStems[[i]] <- AllStems[[i]][!(AllStems[[i]]$Mnemonic %in% dropme), ]}


# --- Generate methods for assigning/estimating Crown Area --------------------
# Method 1, try using the TALLO data, using mixed effects models
# load(paste0(loc_Gdr, "/data/ForestGEO/TALLO/Tallo_MEModelFits.r"))
# fit <- TalloLM_sep_fits$
# CF <- Tallo_rslt[[2]]
# b <- coef(fit)$biome_div
# 
# # Create the function
# CRfunc_temp <- function(clade, DBH, CF){
# 
#   
#   # Values for the intercept depending on the clade
#   b0 <- mutate(tibble(clade), b0 = case_when(clade == "U" ~ b$`(Intercept)`[8],
#                                             clade == "G" ~ b$`(Intercept)`[9],
#                                             clade == "A" ~ b$`(Intercept)`[8],
#                                             TRUE ~ -111111)) %>% pull(b0)
#   # Values for the exponent dependent on the clade
#   b1 <- mutate(tibble(clade), b1 = case_when(clade == "U" ~ b$`log(stem_diameter_cm)`[8],
#                                             clade == "G" ~ b$`log(stem_diameter_cm)`[9],
#                                             clade == "A" ~ b$`log(stem_diameter_cm)`[8],
#                                             TRUE ~ 0)) %>% pull(b1)
#   return( CF*exp(b0+b1*log(DBH)) )
# }
# 
# 
# # Apply the function to the data
# AllStems <- lapply(AllStems, function(x){
#   x <- mutate(x,
#              CRad = CRfunc_temp(Clade, DBH = DBH, CF = CF),
#              CArea = pi*(CRad^2))
#})

# --- Apply method 2 with separate lm for each biome-division -----------------
load(paste0(loc_Gdr, "./data/ForestGEO/TALLO/Tallo_sepLMModelFits.r"))

if(site %in% c("HVDF")){
  b0_A <- TalloLM_sep_fits$TempCon_A$b0
  b1_A <- TalloLM_sep_fits$TempCon_A$b1
  CF_A <- TalloLM_sep_fits$TempCon_A$CF
  
  b0_G <- TalloLM_sep_fits$TempCon_G$b0
  b1_G <- TalloLM_sep_fits$TempCon_G$b1
  CF_G <- TalloLM_sep_fits$TempCon_G$CF
  }else if(site %in% c("SERC", "SCBI")){
    # Extract the parameters estimated earlier
    b0_A <- TalloLM_sep_fits$TempBL_A$b0
    b1_A <- TalloLM_sep_fits$TempBL_A$b1
    CF_A <- TalloLM_sep_fits$TempBL_A$CF
  
    b0_G <- TalloLM_sep_fits$TempBL_G$b0
    b1_G <- TalloLM_sep_fits$TempBL_G$b1
    CF_G <- TalloLM_sep_fits$TempBL_G$CF
    }

AllStems %<>% lapply(., function(x){
  x %<>% mutate(CRad = case_when(Clade %in% c("A", "U") ~
                                  CF_A*exp(b0_A + b1_A*log(.$DBH)),
                                Clade == "G" ~
                                  CF_G*exp(b0_G + b1_G*log(.$DBH)),
                                is.na(Clade) ~ NA_real_,
                                TRUE ~ -(10^9)), CArea = pi*(CRad^2)) %>%
    select(-c(CRad))
  return(x)
  
})


# --- Apply the algorithm to assign canopy level/layer ------------------------
# Go through and perform the algorithm

# THIS ASSUMES THAT THE PLOT IS A RECTANGLE WITH SIDES ALONG THE X and Y Axis!!


# How big is the subplot for the algorithm
subp_dim <- 25 # dimension of each side of the subplot in meters for the
                  # canopy assignment algorithm, currently using Ruger's
                  # 31.25m x 31.25m value

# now we need to make the marks delineating where to assign canopy layers
plot(x = AllStems[[1]]$PX, y = AllStems[[1]]$PY, main = site)

y_min0 <- min(AllStems[[1]]$PY, na.rm = T)
x_min0 <- min(AllStems[[1]]$PX, na.rm = T)
x_max0 <- max(AllStems[[1]]$PX, na.rm = T)
y_max0 <- max(AllStems[[1]]$PY, na.rm = T)

nsub_x <- floor(round(max(AllStems[[1]]$PX, na.rm = T))/subp_dim)
nsub_y <- floor(round(max(AllStems[[1]]$PY, na.rm = T))/subp_dim)

x_min <- floor(min(AllStems[[1]]$PX, na.rm = T)) # minimum x, rounded down 
y_min <- floor(min(AllStems[[1]]$PY, na.rm = T)) # minimum y, rounded down 

# upper bound, note that since we need to have canopy levels assigned in entire
# m x m meter subplots, we may have extra "outside" to trim. For example
# if subplots are 10m x 10m, and stems at positions of y from 0.1 to 92.1
# the stems at y > 90m will be removed since we don't have a complete subplot
# to fill in stems for them.
x_max <- max(nsub_x)*subp_dim
y_max <- max(nsub_y)*subp_dim

# generate a message talking about trimming the stems outside of neatly
# assigned subplots
print(paste0(site, ": x = [", x_min0,", ", x_max0, "], y = [",
             y_min0,", ", y_max0, "]. ","Trimming site to x = [",
             x_min,", ", x_max, "], y = [", y_min,", ", y_max,
             "] to accomadate subplots of: ", subp_dim, "m x ", subp_dim,"m"))

# We need to go through and pick out the stems that are outside the boundary,
dropme <- data.frame(StemID = AllStems[[1]]$StemID, drop = FALSE)
for(df in AllStems){
  if(!identical(df$StemID, dropme$StemID)){
    warning("Order of StemID changed when removing out-of bounds stems!")}
  
  dropme$drop <- dropme$drop |
    ((!is.na(df$PX)) & ((df$PX < x_min) | (df$PX > x_max))) |
    ((!is.na(df$PY)) & ((df$PY < y_min) | (df$PY > y_max)))

  
}

print(paste0(length(which(dropme$drop)),
             " stems removed at ",
             site,
             " for being out of bounds",
             " for subplots used in canopy-level assignments (",
             round(100*length(which(dropme$drop))/length(dropme$drop), digits = 2),
             "%)"))

 
dropme <- dropme$StemID[dropme$drop]

# trying to drop sites after the algorithm sorts through it all, for HVDF
# it should get rid of the swamp
find_dropmeLater <- NULL
print("made dropmelat")
if(site == "HVDF"){
  find_dropmeLater <- AllStems[[2]]$StemID[AllStems[[2]]$DFstatus == "stem_gone"]
  
}

dropmeLater <- NULL
print(paste("201", is.null(dropmeLater)))

AllStemsTemp <- AllStems
for(i2 in 1:length(AllStems)){
  
  x <- AllStems[[i2]]

  # the number of subplots for each dimension, note that some measured stems will be dropped since they are on the edge
  
  x %<>% filter(!(StemID %in% dropme))

  
  lyr_df <- data.frame(StemID = x$StemID,
                      AlgoSubplot = NA,
                      CanopyLvl = NA)

  # Assumes that the area is a perfect quadrilateral, and there are no gaps that should be ignored.
  for(j in 1:nsub_y){
    for(i in 1:nsub_x){
      
      
      
      
     
      # print(paste0(i,", ",j))
      # Find all alive individuals, aren't supposed to be dropped and are in the right place
      sub_plot <- x[which(x$PX >= (i - 1)*subp_dim &
                            x$PX <= i*subp_dim &
                            x$PY >= (j - 1)*subp_dim &
                            x$PY <= j*subp_dim), ]
      
      # HERE!
      if(any(sub_plot$StemID[sub_plot$StemID %in% find_dropmeLater])){
        if(is.null(dropmeLater)){print(paste("NULL at i=",i,", j =", j))}
        dropmeLater <- c(dropmeLater, sub_plot$StemID)
        
      }
      #if(is.null(dropmeLater)){print(paste("Line 237", is.null(dropmeLater)))}
      
      sub_plot %<>% filter(sub_plot$RametStatus == "A" &
                             !is.na(sub_plot$DBH) &
                             sub_plot$DropBegin == 0)

      
      
      #if(any(sub_plot$StemID %in% stopme)){print(paste(i, j))}
      
      if(!(nrow(sub_plot) > 0)){
        next
        }
      
      # Add columns for the Algorithm subplot and canopy level/layer
      sub_plot %<>% mutate(AlgoSubplot = paste(i, j, sep = "-"), CanopyLvl = NA)
      

      
      # Reorder this data by crown radius, largest at beginning
      sub_plot <- sub_plot[order(sub_plot$CArea, decreasing = T),]
      
      # Start with the "total area" being 0, and the layer being 1
      lyr_area <- 0
      l <- 1
      
      # Go through all the stems in the subplot
      for(s in 1:nrow(sub_plot)){
        # If there are no individuals in that layer OR the area of the crown, if greater than or equal to 50% of
        # the crown area goes into the layer, keep it in the same layer, assign the stem automatically to the layer
        if((length(which(sub_plot$CanopyLvl == l)) == 0) |
           ((subp_dim^2) - lyr_area) >= (sub_plot$CArea[s]/2)){
          # Assign layer and add to the "accounted area
          sub_plot$CanopyLvl[s] <- l
          lyr_area <- lyr_area + sub_plot$CArea[s]
          
        # If <50% of the next stem's crown area fits in the layer, then we move on to the next layer
        }else{
            l <- l + 1 # Go on to the next layer
            sub_plot$CanopyLvl[s] <- l
            lyr_area <- sub_plot$CArea[s] # the "layer area" is now reset, only that individual is in the layer now
        }
      } 
      # Here we go and include the layer assignments and write down the subplot it was placed in
      lyr_df[which(lyr_df$StemID %in% sub_plot$StemID), ] <-
        sub_plot[, which(colnames(sub_plot) %in% c("StemID", "CanopyLvl", "AlgoSubplot"))]

    }
  }

  
  # Add it into the actual data now
  x <- left_join(x, lyr_df, by = "StemID")
  AllStemsTemp[[i2]] <- x
}

AllStems <- AllStemsTemp
dropmeLater <- unique(dropmeLater)


# ### Prep and send out the data to the next script #############################################################################
# Ok, now lets trim out some of the columns that we no longer need
AllStems <- lapply(AllStems, function(x){
  x <- x[, which(colnames(x) %in% c(cn2k, "DropParamtrz", "AlgoSubplot", "CanopyLvl"))]
  
  
  if(length(dropmeLater) > 0){
  # Add the flagging for those stems in weird spots (e.g. HVDF swamp)
  x$AlgoSubplot[x$StemID %in% dropmeLater] <- "REMOVE"

  # dropping Subplots in problematic areas
  x <- x[-1*c(which(x$StemID %in% dropmeLater)),]}

  return(x)
})

# print message about dropping problematic areas
if(!is.null(dropmeLater)){
  cat(paste0("Dropping ", length(dropmeLater)," stems (~",
            round( 100*length(dropmeLater)/nrow(AllStemsTemp[[1]]) ,1),
             "%) due to being in areas",
             " for exclusion (e.g. deer exclosure, swamp)."))}


# plot
plot.df <- AllStemsTemp[[1]] %>% mutate(Flag = case_when(!(StemID %in% dropmeLater) ~ "Ignored",
                                              StemID %in% find_dropmeLater ~ "Excluded",
                                              StemID %in% dropmeLater ~ "Trimmed",
                                              TRUE ~ "ERROR!"))

myplot <- ggplot(plot.df, aes(x = PX, y = PY, colour = Flag)) +
  geom_point() + ggtitle(site)
print(myplot)

for(i in 1:length(AllStems)){
  print(dim(AllStems[[i]]))
}

# Location of output for the canopy level assigned datasets
loc_2 <- "/data/ForestGEO/processed/canopy_assign/"
# Not the best solution to use assign, but I wanted to make sure I could individually retrieve data for each site.
assign(paste0(site,"_AllStemsCL"), AllStems)
save(list = paste0(site,"_AllStemsCL"),
     file = paste0(loc_Gdr, loc_2, site, "_AllStemsCL.r"))

}

#'
#' This is just the initial summary for the first three sites.
#'
#'





















