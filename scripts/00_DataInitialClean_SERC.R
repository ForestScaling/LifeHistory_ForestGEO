# Script: 00_DataInitialClean_SERC.R
# Author: JMR 
# Last Update: 12023-01-11
# Description:
# File in-taking the data for SCBI and spitting out a standardized version so that it is harmonized
# with the other ForestGEO sites. Outputs two files, a "taxa" file, containing the information to
# correlate stems to actual species of interest and (if you want to drop particular problematic 
# lianas for instance) as well as the AllStem file, which contains all the stems measured across 
# the censuses. Also runs QC checks for the data.
# Future:
# I need to check for what the different codes mean for this particular site, still waiting to hear
# back from SERC.
# ### Setup #######################################################################################
# Modifications to specific datasets, before initial processing
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"

# Setting/finding location of the script, so the WD is automatically within GITHub local location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

#' Commet
# Bringing in the packages:
library("tidyverse")
library("magrittr")

# Get the commonly used functions
source("scripts/00_Functions.R")

# Site Name:
site <- "SERC"

# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/", site))[
  str_which(list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/", site)), pattern = "_AllMeas_")]

# Create a file for storing all of the Stem files
Raw_Stems <- lapply(stem_files, FUN = reader1, loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))

# Read in the species table
Spp_table <- read.delim(paste0(loc_Gdr, "/data/ForestGEO/raw/", site,"/",
                               site,"_TaxaGenusSp_Report.txt"), sep =  "\t", header = T)

# ### Beginning Manipulations of Species table ####################################################
# First verifying that all species in the data actually have a corresponding row in the spp table
if(length(which(sapply(Raw_Stems, FUN = function(z){
  length(which(!(z$Mnemonic %in% Spp_table$Mnemonic)))}) == 0)) != length(Raw_Stems)){
  warning("At least one specieces identification code in the Stem data 
          is NOT represented in the table of taxa!")}
# No discrepancies of note

# The Taxanomic table seems to be in pretty good shape, the only exception seems to be that IDLevel
# is treating 'multiple' and 'genus' synonymously under 'genus'. I am going to take the conserva-
# -tive approach and assign "multiple" to all genus sp. that have congeners already present in the
# species table.

Spp_table <- as_tibble(Spp_table)

# Rename the column to be in accordance
colnames(Spp_table)[which(colnames(Spp_table) == "Subspecies")] <- "SubSpecies" 

# Sambucus nigra has a ssp. but they write it as spp. (plural species), we are also going to update
# the latin/scientific name
# Another reason I think ssp. are overrated!
Spp_table$SubSpecies[which(Spp_table$SubSpecies == "spp. canadensis")] <- "ssp. canadensis" 

# Change blank subspecies to NA
Spp_table$SubSpecies[which(Spp_table$SubSpecies == "")] <- NA

# Add in the Latin names
Spp_table %<>% mutate(Latin =
                        case_when(!is.na(.$SubSpecies) ~ {paste(.$Genus, .$Species, .$SubSpecies)},
                      TRUE ~ {paste(Genus, Species)}))

# Problematic assignments, do not include these at all, may include lianas which canopy cannot be
# really modeled effectively for crown class assignment due to them not supporting their own
# weight. Others may be included for crown canopy levels, but not have their parameters estimated.
# Also add whether taxon is gymnosperm 'G' or angiosperm 'A', or unknown 'U'
Spp_table %<>% mutate(
  DropBegin = case_when(Family  == "Vitaceae" | Genus == "Toxicodendron" ~ 1,
                                            TRUE ~ 0),
  DropParamtrz = case_when(IDlevel %in% c("species", "subspecies") ~ 0,
                           TRUE ~ 1),
  Clade = case_when(Family == "Unknown" ~ "U",
                    Family %in% c("Pinaceae", "Taxaceae", "Cupressaceae") ~ "G",
                    TRUE ~ "A"))

# ### Formatting the Stem data to conform with other sites ########################################
Raw_Stems <- lapply(Raw_Stems, function(y){
  y$SubSpecies[which(y$SubSpecies == "")] <- NA
  return(y)})

# Let's go back and sync up the Latin names again
Raw_Stems <- lapply(Raw_Stems, function(x){
  x <- x[, -2]
  x <- left_join(x, Spp_table[,c(11, 6)], by = "Mnemonic")
  x <- x[, c(1, 18, 2:17)]
  return(x)})

# Beginning to add new columns in for recruitment and status  (genet and ramet) -------------------
# We need to check what are all the current statuses and codes used in the SERC data:
codes <- unique(unlist(lapply(Raw_Stems, function(x){
  return(str_split(x$Codes, pattern = ";"))
  })))
statuses <- unique(unlist(lapply(Raw_Stems, function(x){return(x$Status)})))

# Show the unique codes and statues.
cat("List of unique codes:\n")
print(codes)
cat("List of unique statuses:\n")
print(statuses)

# Including the recruit status.
# Values are either:
# RR - ramet is recruit, the genet was recorded at t-1, but new ramet at t
# GR - genet recruited, the entire genet has no measurements prior to t, note that all RR are
#   GR but not all GR are RR
# NA - first census, or not found
# NO - not a recruit (repeat measure)

# Careful, it assumes that the censuses are in order
Raw_Stems <- lapply(Raw_Stems, function(x){return(mutate(x, RecruitStatus = "NO"))})
Raw_Stems[[1]]$RecruitStatus <- NA
Raw_Stems %<>%
  any1b4(out_col = "RecruitStatus",
         t_name = "Census",
         val_r = "RR",
         check_col = "StemID") %>% 
  any1b4(out_col = "RecruitStatus",
         t_name = "Census",
         val_r = "GR",
         check_col = "TreeID")

# Adding in a column for ramet Status
# The following are the possible values it can take.
# A - ramet alive
# D - ramet dead
# M - missing, no data for the individual, but was present last time.
# P - Prior, individual was not present in this survey (should not be applicable at this stage since SERC didn't store their data that way.)

# Right now whether a ramet is alive or dead is going to be set solely based on the status column, with broken stems considered dead for simplicity.
# Broken stems will be heandled with more precision later once I have a response from the SERC people. Stems with the code SL (only in census 2) will
# be given 'M' since these lack DBHs and are going to assumed to mean Stem-lost (we will verify this with the metadata, but this will be good to just
# provide a test case)

# Adding in the column itself
Raw_Stems %<>% lapply(., function(x){
  x %<>% mutate(RametStatus = case_when(as.vector(CodesLogic(col = .$Codes, "SL", sep = ";")) ~ "M",
                                        .$Status == "alive" ~ "A",
                                        .$Status %in% c("stem dead", "dead", "broken below") ~ "D",
                                        TRUE ~ NA_character_))
  # Stems with an "SL" code are considered lost (M)
  # Stems with status of stem dead, dead, of broken below are "D"
  # stems status alive is "A"
  # NA otherwise
  return(x)
  })

# Now we need to add data about the genet
# A - genet alive, some part of the tree/shrub is alive
# D - the entire tree/shrub is dead
# M - missing, no data for the individual, but was present in a previous survey
# P - Prior, individual was not present in this survey but is present in future survey(s)
# F - category which is a subset of A where the stem(s) alive are <10mm DBH, IS NOT BEING INCLUDED
#     RIGHT NOW!

Raw_Stems %<>% lapply(., genetAssign, evid_rcrt = T)


# Add in the "missing" stems, ramets with status "M" and "P".
Raw_Stems %<>% MP_adder(.)

# QC v2, all in one function
QC <- QC_Check(Raw_Stems,
               do.spp_chck = T,
               do.chng_HOM = T,
               do.zomb_ram = T)

Raw_Stems <- QC$QCdata

# --- Fill in values for the Stem (Main vs secondary) ---------------------------------------------
# Use a consistent algorithm.
Raw_Stems %<>% stemRankAsgn(., col_dbh = "DBH", col_out = "Stem", consider_if = "A",
                            col_con = "RametStatus", col_grp = "TreeID", col_ind = "StemID")

# ### Export the data back out for use in the actual analysis #######################################################################################
# I frequently am using chunks of this code in other data harmonization scripts for this project, so this is meant to change the
# name into a more identifiable one for later if need be.
SERC_AllStems <- Raw_Stems 
SERC_TaxaTable <- Spp_table

# Save SERC Stemdata
save(SERC_AllStems, file = paste0(loc_Gdr,"/data/ForestGEO/processed/standardized/SERC_AllStemsOLD.r"))

# Save the Taxatable
save(Spp_table, file = paste0(loc_Gdr,"/data/ForestGEO/taxa/SERC_TaxaTable.r"))


# 
# 
# 
# 
# 
# 
# 
# 
# # --- QC Zombies, Alternate POM & Changing species ID -----------------------------------------------------------------------------------------------
# # Create the dataframe to hold the individual StemID across time
# tb <- ListWiden(List = Raw_Stems, by = "StemID", kcols = c("GenetStatus", "RametStatus", "DBH", "Latin", "Codes", "HOM", "TreeID"))
# tb$Change_1 <- 0
# tb$Change_2 <- 0
# 
# # check that species did not change ID over time
# if(length(which(tb$Latin_1 != tb$Latin_2)) != 0){
#   warning(paste("There are",length(which(tb$Latin_1 != tb$Latin_2)), "StemIDs that change species ID over time!" ))}
# 
# # Check that those with a POM shift have code "A"
# if(nrow(tb[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";"))),]) > 0){
#   warning(paste("There are",
#                 nrow(tb[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";"))),]),
#                 "StemIDs with differing POM but lacking code 'A', adding the additional code now."))
#   
#   # Adding a marker column to denote which ones were adjusted
#   tb$Change_2[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))] <- 1
#   
#   # Adding the code
#   tb$Codes_2[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))] <-
#     paste(tb$Codes_2[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))], "A", sep = ";")
#   
# 
# }
# 
# # Check for zombie stems
# if(nrow(tb[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" ),]) > 0){
#   warning(paste("There are",
#                 nrow(tb[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" ),]),
#                 "StemIDs whose RametStatus goes from dead to alive. Changing their RametStatus to be 'A'."))
#   
#   # Adding a marker so that we know they were changed
#   tb$Change_1[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" )] <- 1
#   # Changing the ramet status to alive so that they aren't zombies
#   tb$RametStatus_1[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" )] <- "A"
#   # Change the genet status, since if the ramet is alive so is the genet
#   tb$GenetStatus_1[which(tb$RametStatus_1 == "A")] <- "A"
# 
# }
# 
# # Updating Genet status, we are just going to rerun the genet status algorithm from earlier, and record how many changed
# # Create a placeholder variable storing which TreeID at t1 have more than 1 
# x <- tb %>% group_by(TreeID_1) %>% summarize(uniq = length(unique(GenetStatus_1)))
# x <- x$TreeID_1[which(x$uniq != 1)] # Find all those TreeID with more than 1 GenetStatus (one was changed earlier)
# 
# # change the genet status to 'A' for all of these, and note that their data was edited.
# tb$GenetStatus_1[which(tb$TreeID_1 %in% x)] <- "A"
# tb$Change_1[which(tb$TreeID_1 %in% x)] <- 1
# 
# # Now is time to reincorporate the Edits for QC
# # We will also add a code of "Edit" for all those that were modified so that we can keep track of them.
# Raw_Stems <- lapply(Raw_Stems, FUN = function(x){
#   i <- unique(x$Census)
#   
#   # Separate only the columns of interest
#   tb_use <- tb[,c(1, str_which(string = colnames(tb), pattern = paste0("_",i,"$") ))]
#   
#   # Edit the column names so they are generic again
#   colnames(tb_use) <- sub(paste0("_", i, "$"), replacement = "", x = colnames(tb_use))
#   
#   # Drop rows that have "NA" as a species (they weren't represented in the original df)
#   tb_use <- tb_use[which(!is.na(tb_use$Latin)),]
#   
#   # Add the "Edit" code to those that were changed
#   tb_use$Codes[which(tb_use$Change == 1)] <- paste(tb_use$Codes[which(tb_use$Change == 1)], "Edit", sep = ";")
#   
#   # Let's re-add these data back to Raw_Stems
#   x2 <- colnames(x)
#   x <- x[, which((colnames(x) == "StemID" )| !(colnames(x) %in% colnames(tb_use)))]
#   x <- left_join(x, tb_use, by = "StemID")
#   # Rearrange the columns to match again, also drops the Change column
#   
#   x <- x[,match(x2, colnames(x))]
#   
#   return(x)
# })
# 
# # Checking that any genets with ramet recruits at t+1 aren't considered a dead genet at t
# z <- unique(Raw_Stems[[2]]$TreeID[which(Raw_Stems[[2]]$RecruitStatus == "RR")]) # TreeID with RR at t=2
# z <- (Raw_Stems[[1]]$TreeID %in% z) & # TreeID matches one with recruit ramets at t=2
#                                  Raw_Stems[[1]]$GenetStatus == "D"
# if(nrow(Raw_Stems[[1]][which(z),]) > 0){ # Genet was "dead" at t = 1
#   warning(paste("There are", nrow(Raw_Stems[[1]][which(z),]), "StemIDs that have a genet which was 'dead' at t = 1 but still had recruits at t = 2.
# Considering the genets to be alive at t = 1, you can always override this if you think they were just
# missed."))
#   
#   # Going through and reassigning the values of genet status
#   Raw_Stems[[1]]$GenetStatus[which(z)] <- "A"
#   Raw_Stems[[1]]$Codes[which(z & !(CodesLogic(Raw_Stems[[1]]$Codes, c = "Edit", sep = ";")))] <- 
#     paste(Raw_Stems[[1]]$Codes[which(z & !(CodesLogic(Raw_Stems[[1]]$Codes, c = "Edit", sep = ";")))], "Edit", sep = ";")
#   
# }
# 
# 
# # Because some places had no codes earlier, there are some codes like ";A:B" or maybe some like "A;B;", lets get rid of any
# # leading or trailing semi-colons, (this was an issue at HVDF, but I don't think it was present here)
# Raw_Stems <- lapply(Raw_Stems, function(x){
#   
#   # For the ";B;C"
#   if(length(str_which(x$Codes, pattern = "^;")) > 0){
#     x$Codes <- sub(x = x$Codes, pattern = "^;", replacement = "")
#   }
#   
#   # For the "A;B;"
#   if(length(str_which(x$Codes, pattern = "^;")) > 0){
#     x$Codes <- sub(x = x$Codes, pattern = "^;", replacement = "")
#   }
#   
#   return(x)
#   
# })
# 
# # --- Fill in values for the Stem (Main vs secondary) ---------------------------------------------
# # Use a consistent algorithm.
# Raw_Stems %<>% stemRankAsgn(col_dbh = "DBH", col_out = "Stem", consider_if = "A",
#                             col_con = "RametStatus", col_grp = "TreeID", col_ind = "StemID")
# 
# 
# 


