# Script: 00_DataInitialClean_ZOFN.R
# Author: JMR 
# Last Update: 12023-01-09
# Description:
# File in-taking the data for ZOFN and spitting out a standardized version so that it is harmonized
# with the other ForestGEO sites. Outputs two files, a "taxa" file, containing the information to
# correlate stems to actual species of interest and (if you want to drop particular problematic
# lianas for instance) as well as the AllStem file, which contains all the stems measured across
# the censuses.

# ### Setup #######################################################################################
# Modifications to specific datasets, before initial processing
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"

# Setting/finding location of the script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

# Bringing in the packages:
library("tidyverse")
#library("magrittr")

# Get the commonly used functions
source("scripts/00_Functions.R")

# Site Name:
site <- "ZOFN"

# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/", site))[
  str_which(list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/", site)), pattern = "_AllMeas_")]

# Create a file for storing all of the Stem files
Raw_Stems <- lapply(stem_files,
                    FUN = reader1,
                    loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))

# Read in the species table
Spp_table <- read.delim(paste0(loc_Gdr, "/data/ForestGEO/raw/", site,"/",
                               site,"_TaxaGenusSp_Report.txt"), sep =  "\t", header = T)

# ### Taxa Table Harmonization ####################################################################
# Check for any missing /unaccounted for species
# Holy cow, that's only a handful of spp.!
if(length(which(sapply(Raw_Stems, FUN = function(x){
  length(which(!(x$Mnemonic %in% Spp_table$Mnemonic)))}) == 0)) != length(Raw_Stems)){
  warning("At least one specieces identification code in the Stem data is NOT represented in the 
         table of taxa!")} # No discrepancies of note

# Let's change the colname Subspecies to SubSpecies, for consistency with the stem data
colnames(Spp_table)[which(colnames(Spp_table) == "Subspecies")] <- "SubSpecies"

# Adding data for potentially problematic taxa, ZOFN has so few, and doesn't  have Lianas that 
# this probably be moot for this site
Spp_table$DropBegin <- 0 
Spp_table$DropParamtrz <- 0

# Adding a categorization of conifers and Angiosperms
Spp_table$Clade <- "A"
Spp_table$Clade[which(Spp_table$Family %in% c("Pinaceae", "Taxaceae", "Cupressaceae"))] <- "G"


# ### Stem Harmonization ##########################################################################
# ZOFN lacks "Tag" numbers, I am just going to fill in this with TreeID and a "TAG"
# It also lacks information for StemID, but for this I will use the TreeID + tag number with
# filling zeros, e.g. TreeID 12345 Tag 12 will be 12345012
# Also changing subspecies to NA for "blanks"
Raw_Stems <- lapply(Raw_Stems, function(x){
  # Here is Tag "filled in"
  x$Tag <- paste0("TAG", x$TreeID)
  
  # going through and creating the stem ID
  x$StemID <- paste0(x$TreeID, str_pad(x$StemTag, width = 3, side = "left", "0"))
  
  # Subspecies
  x$SubSpecies[which(x$SubSpecies == 0)] <- NA
  
  # Changing the DBH units to cm to match others
  x$DBH <- x$DBH/10
  
  return(x)
})

# Let's take a look at what statuses and codes are present
codes <- unique(unlist(lapply(Raw_Stems, function(x){
  return(str_split(x$Codes, pattern = ";"))
})))
statuses <- unique(unlist(lapply(Raw_Stems, function(x){return(x$Status)})))

# Show the unique codes and statues.
cat("List of unique codes:\n")
print(codes)
cat("List of unique statuses:\n")
print(statuses)
# Not very helpful...

# --- Adding in the Ramet status information, genet status will be added post QC check ------------
# Adding in a column for ramet Status
# The following are the possible values it can take.
# A - ramet alive
# D - ramet dead
# M - missing, no data for the individual, but was present last time.
# P - Prior, individual was not present in this survey (should not be applicable at this stage
# since ZOFN didn't store their data that way.)

# Adding in the column itself
Raw_Stems <- lapply(Raw_Stems, function(x){
  x <- mutate(x, RametStatus = NA)
  x$RametStatus[which(x$Status == "alive")] <- "A" 
  x$RametStatus[which(x$Status %in% c("stem dead", "dead"))] <- "D"
  #x$RametStatus[which(CodesLogic(x$Codes, "L"))] <- "M" # Stems with an 'L' code are left as miss.
  return(x)
})

# Genet status:
# A - genet alive, some part of the tree/shrub is alive
# D - the entire tree is dead
# M - missing, no data for the individual, but was present in a previous survey
# P - Prior, individual was not present in this survey but is present in future survey(s)
# * Normally there is an "F" category... but not used here.
Raw_Stems <- lapply(Raw_Stems, function(x){
  x %<>%
    group_by(TreeID, ) %>%
    summarize(sts = paste(unique(RametStatus), collapse = ";"),
              GenetStatus = case_when(
                CodesLogic(sts, "A") ~ "A",
                CodesLogic(sts, "D") ~ "D",
                CodesLogic(sts, "M") ~ "M",
                sts == "P" ~ "P",
                TRUE ~ "Z-ERROR")) %>%
    select(!(sts)) %>%
    right_join(x) %>% 
    {select(., match(c(colnames(x), "GenetStatus"), colnames(.)))}
  
  return(x)
})

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


# --- QC Checks -----------------------------------------------------------------------------------
# Probably garbage come back to later
print("stop 12023-01-11, disregard QC")
# Create the widened list for easy comparison
tb <- ListWiden(List = Raw_Stems, by = "StemID", kcols = c("GenetStatus", "RametStatus", "DBH",
                                                           "Latin", "Codes", "HOM", "TreeID"))
tb$Change_1 <- 0
tb$Change_2 <- 0

# check that species did not change ID over time
if(length(which(tb$Latin_1 != tb$Latin_2)) != 0){
  warning(paste("There are",length(which(tb$Latin_1 != tb$Latin_2)),
                "StemIDs that change species ID over time!" ))}

# Check that those with a POM shift have code "A"
if(nrow(tb[which((tb$HOM_1 != tb$HOM_2) &
                 !(CodesLogic(tb$Codes_2, c = "A", sep = ";"))),]) > 0){
  warning(paste("There are",
                nrow(tb[which((tb$HOM_1 != tb$HOM_2) &
                                !(CodesLogic(tb$Codes_2, c = "A", sep = ";"))),]),
                " StemIDs (... so all of them) with differing POM but lacking code 'A', ",
                "adding the additional code now."))
  
  # Adding a marker column to denote which ones were adjusted
  tb$Change_2[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))] <- 1
  
  # Adding the code
  tb$Codes_2[which((tb$HOM_1 != tb$HOM_2) & !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))] <-
    paste(tb$Codes_2[which((tb$HOM_1 != tb$HOM_2) &
                             !(CodesLogic(tb$Codes_2, c = "A", sep = ";")))], "A", sep = ";")
}
# Check for zombie stems
if(nrow(tb[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" ),]) > 0){
  warning(paste("There are",
                nrow(tb[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" ),]),
                "StemIDs whose RametStatus goes from dead to alive. Changing their RametStatus to be 'A'."))
  
  # Adding a marker so that we know they were changed
  tb$Change_1[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" )] <- 1
  # Changing the ramet status to alive so that they aren't zombies
  tb$RametStatus_1[which(paste(sep = "_", tb$RametStatus_1, tb$RametStatus_2) == "D_A" )] <- "A"
  # Change the genet status, since if the ramet is alive so is the genet
  tb$GenetStatus_1[which(tb$RametStatus_1 == "A")] <- "A"
}








# Post QC
# Assign Stem (Main vs Secondary)
# Note that it is overriding their previous assignments 
warning("Overriding site's own stem assignments if present. Only live stems considered, largest 
  stem assigned 'Main' others assigned 'Secondary' (NA-s are considered last). This will NOT be 
  denoted with an 'Edit' code. Note that a stem's assignemt can change over time if a secondary 
  surpases the main at a later time.")

Raw_Stems %<>% stemRankAsgn(col_dbh = "DBH", col_out = "Stem", consider_if = "A",
                            col_con = "RametStatus", col_grp = "TreeID", col_ind = "StemID")


















