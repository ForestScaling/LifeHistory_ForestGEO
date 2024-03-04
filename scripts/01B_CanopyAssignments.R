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
packages <- c( "lme4", "rlist", "rstudioapi", "tidyverse", "magrittr")
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

# Read in the taxa table for this site
load(paste0(loc_Gdr, "/data/ForestGEO/taxa/",
            site, "_TaxaTable.r"))

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
subp_dim <- 25
cat("Going through algorithm with ", subp_dim, "m x ", subp_dim,
    "m dimensions. Please ensure that this matches the grid cells.\n",
    sep = "")

print(paste0("For: ",site, " dim of: ", dim(AllStems[[1]])))
AllStems <- lapply(AllStems, FUN = function(X){
  
  # Add column for the CanopyLvl
  X %<>% mutate(CanopyLvl = NA)
  
  # Begin looping through the subplpots at census X
  for(i in unique(X$AlgoSubplot)){
    # Skip the NA values
    if(is.na(i)){
      next
      }
    
    # Get stems only for that subplot
    sub_plot <- X %>%
      filter(AlgoSubplot == i & !is.na(AlgoSubplot)) %>% # must be in subplot
      filter(RametStatus == "A" & !is.na(DBH)) %>% #must be alive
      arrange(desc(CArea))
    
    # Skip forward if no live stems present
    if(!(nrow(sub_plot) > 0)){next}
    
    # set area initially to 0 and start with layer 1
    lyr_area <- 0
    l <- 1
    
    # Go through all the stems in the subplot
    for(s in 1:nrow(sub_plot)){
      # If there are no individuals in that layer OR the area of the crown,
      # if greater than or equal to 50% of the crown area goes into the layer,
      # keep it in the same layer, assign the stem automatically to the layer
      
      if((length(which(sub_plot$CanopyLvl == l)) == 0) |
         ((subp_dim^2) - lyr_area) >= (sub_plot$CArea[s]/2)){
        # Assign layer and add to the "accounted area
        sub_plot$CanopyLvl[s] <- l
        lyr_area <- lyr_area + sub_plot$CArea[s]
        
        # If <50% of the next stem's crown area fits in the layer,
        # then we move on to the next layer
      }else{
        l <- l + 1 # Go on to the next layer
        sub_plot$CanopyLvl[s] <- l
        lyr_area <- sub_plot$CArea[s]} # the "layer area" is now reset
      
      
    } # end loop through stems of the subplot
    
    # Add in the values of the subplot back to the df
    X[which(X$StemID %in% sub_plot$StemID),] <- sub_plot
    
  } # end loop through subplots, begin next census
  
  return(X)
}) # end lapply


# ### Generate values #########################################################
# --- Canopy Level ------------------------------------------------------------
# number of canopy levels is set automatically
df <- AllStems[[1]] %>%
  group_by(CanopyLvl) %>%
  summarize(n = n()) %>% 
  filter(!is.na(CanopyLvl)) %>% 
  arrange(CanopyLvl)

print(site)
print(df)

# proportional cutoff
k <- 0.1 # must have at least k times as many stems in the lowest layer



if(nrow(df) != 1){ # verify that there are at least two canopy levels
  while(df$n[nrow(df)] < df$n[nrow(df) - 1]*k){ 
    # check that the lowest CL has at least k times as many stems as the
    # next CL up
    
    # add the stems from the lowest CL to the next one up
    df$n[nrow(df) - 1] <- df$n[nrow(df) - 1] + df$n[nrow(df)]
    
    df <- df[-nrow(df),] # remove the lowest layer
    
    # Check that there are at least two layers remaining,if not, stop grouping
    if(nrow(df) < 2){
      break
    }
  }
}

# Message notifying user if number of canopy layers was reduced.
if(nrow(df) != max(AllStems[[1]]$CanopyLvl, na.rm = T)){
  # notify if reduced
  cat(paste0("number of canopy levels reduced from ",
             max(AllStems[[1]]$CanopyLvl, na.rm = T),
             " to ", nrow(df),
             " due to insufficient stems in lower levels.\n")) 
}else{
  # if not reduced
  cat(paste0(nrow(df), " canopy Levels used.\n"))
}

# set the number of canopy layers
cat("Manually setting number of canopy Levels to 2.\n")
n_Lvls <- 2
#n_Lvls <- nrow(df)

# group together all the rows that have insufficient stems
AllStems <- lapply(AllStems, function(df){
  df %<>% mutate(CanopyLvl = replace(CanopyLvl,
                                     CanopyLvl > n_Lvls,
                                     n_Lvls))
  return(df)
})




# trim the unnecessary variables
AllStems <- AllStems %>% lapply(., FUN = function(X){
  X <- X %>% select(all_of(c(cn2k, "CanopyLvl")))
  return(X)
})

# visual check
print(site)
print(head(AllStems[[1]]))

# ### Prep and send out the data to the next script ###########################
# Location of output for the canopy level assigned datasets
loc_2 <- "/data/ForestGEO/processed/canopy_assign/"
save(list = paste0("AllStems"),
     file = paste0(loc_Gdr, loc_2, site, "_AllStemsCL.r"))



}
# End