# Script: 00_DataInitialClean_HVDF.R
# Author: JMR
# Last Update: 12023-01-17
# Description:
# File in-taking the data provided by D. Orwig and spitting out a standardized version
# so that it is harmonized with the other ForesTGEO sites.
# ### Setup #######################################################################################
if(!file.exists("./scripts/AA_LocationManager.R")){setwd("..")}
source("./scripts/AA_LocationManager.R")

# Bringing in the packages:
for(i in c("tidyverse", "stringr", "magrittr")){
  library(i, character.only = T)
}

# Get the commonly used functions
source("scripts/00_Functions.R")

# Bringing in the HF Data from Dr. D. Orwig
load(paste0(loc_Gdr, "/data/ForestGEO/raw/HVDF/harvardforest.stem2.rdata"))
load(paste0(loc_Gdr, "/data/ForestGEO/raw/HVDF/harvardforest.stem1.rdata"))
load(paste0(loc_Gdr, "/data/ForestGEO/raw/HVDF/harvardforest.spptable.rdata"))

# Put the relevant site measurements in a single list for easier of manipulation
Raw_Stems <- list(harvardforest.stem1, harvardforest.stem2)

# Generating a species table
Spp_table <- harvardforest.spptable

# Making the simple edit to tree # 33363, supposedly went from dbh of 29.5 -> 308. cm
# within ~ 5yr, probably just a typo meant to be 30.8 not 308.
Raw_Stems[[2]]$dbh[Raw_Stems[[2]]$stemID == 33363] <- 30.8

# Fixing a name typo, Viburnum acerfolium to V. acerifolium
Spp_table$Species[Spp_table$sp == "vibuac"] <- "acerifolium"

# ### Harmonize the taxanomic table ###############################################################
# First verifying that all species in the stem data have a corresponding row in the spp table
if(length(which(sapply(Raw_Stems, FUN = function(z){
  length(which(!(z$sp %in% Spp_table$sp)))
  }) == 0)) != length(Raw_Stems)){
  warning(paste0("At least one specieces identification code in the Stem ",
                 "data is NOT represented in the table of taxa!"))}


# Adding No. column, rearranging columns and renaming them so they are in agreement w/ other sites
Spp_table %<>%
  mutate(No. = 1:nrow(.)) %>%
  select(No., Family, Genus, Species, subsp, sp, IDLevel, Authority, syn, SpeciesID, Latin) %>%
  rename(SubSpecies = subsp, IDlevel = IDLevel, PriorNames = syn, Mnemonic = sp)
  
# Adding Drop @ beginning, Drop-parametrize columns and column for over-arching clade
# Toxicodendron is being flagged since it is a liana
Spp_table %<>%
  mutate(DropBegin = case_when(Genus == "Toxicodendron" ~ 1,
                               TRUE ~ 0),
         DropParamtrz = case_when(IDlevel == "multiple" ~ 1,
                                  IDlevel %in% c("species", "genus") ~ 0,
                                  TRUE ~ NA_real_),
         Clade = case_when(Mnemonic == "deadhw" ~ "A",
                           Mnemonic == "deadsw" ~ "G",
                           Family == "Unknown" ~ "U",
                           Family %in% c("Pinaceae", "Taxaceae", "Cupressaceae") ~ "G",
                           TRUE ~ "A"))


# --- Integrating the species data into the stem-data ---------------------------------------------
# Adding columns for Stem (main vs. secondary) and formatting positions and column names
Raw_Stems <- lapply(Raw_Stems, function(x){
  x %<>%
    rename(Mnemonic = sp) %>%
    # Add in the spp data from the taxa table
    left_join(., y = select(Spp_table, SubSpecies, Mnemonic, IDlevel, SpeciesID, Latin),
              by = "Mnemonic") %>%
    # Add Stem and No. columns
    mutate(No. = 1:nrow(.),
           Stem = NA) %>%
    # Column repositioning
    select(No., Latin, Mnemonic, SubSpecies, quadrat, gx, gy, treeID, tag, stemID, StemTag,
           CensusID, dbh, hom,  ExactDate, codes, Stem, status, DFstatus) %>%
    # Column renames
    rename(Quadrat = quadrat, PX = gx, PY = gy, TreeID = treeID, Tag = tag, StemID = stemID,
           Census = CensusID, DBH = dbh, HOM = hom, Date = ExactDate, Codes = codes,
           Status = status)
  return(x)
  })

# Change some of the NA to "" and visa versa.
Raw_Stems <- lapply(Raw_Stems, function(x){
  x$Codes <- replace_na(x$Codes, "")
  x$Stem[which(x$Stem == "")] <- NA
  return(x)
  })

# Let's look at all the statuses and the codes
codes <- unique(unlist(sapply(Raw_Stems, function(z){
  return(unique(str_split(z$Codes, ";")))
  })))
statuses <- unique(unlist(lapply(X = Raw_Stems, FUN =  function(z){
  return(unique(z$Status))
  })))

# Print out the codes to see:
cat("List of unique codes:\n")
print(codes)
cat("List of unique statuses:\n")
print(statuses)

# --- Adding in the status columns (genet and ramet) & the recruitment ----------------------------
# We need to add in some other standardized columns
# we are going to add a standardized status column 'RametStatus' values are either:
# A - ramet alive
# D - ramet dead
# M - missing, no data for the individual, but was present last time.
# P - Prior, individual was not present in this survey

# Add the column, fill in NA for now, we can check which cases haven't been addressed by getting NA
Raw_Stems <- lapply(Raw_Stems, function(x){
  x %<>% mutate(RametStatus = case_when(Status == "P" ~ "P",
                                        CodesLogic(col = .$Codes, "DT") |
                                          CodesLogic(col = .$Codes, "D") ~ "D",
                                       (is.na(Date) & is.na(DBH) & is.na(Census)) |
                                         (CodesLogic(col = .$Codes, "Miss")) ~ "M",
                                       TRUE ~ "A"))
  return(x)
  })



# We also need to have a 'GenetStatus' column. This has information on the entire tree 
# A - genet alive, some part of the tree/shrub is alive
# D - the entire tree is dead
# M - missing, no data for the individual, but was present in a previous survey
# P - Prior, individual was not present in this survey but is present in future survey(s)
# * Normally there is an "F" category which is a subset of A where the stem(s) alive are 
# <10mm DBH, however at HVDF all stems are greater than 10mm DBH

# Add the column, fill in NA for now, we can check which cases haven't been addressed by getting NA
Raw_Stems <- lapply(Raw_Stems, function(x){
  x <- genetAssign(x, evid_rcrt = T, resprt_code = "S")
  return(x)
})


# Set a column up for recruitment
# Values are either:

# RR - ramet is recruit, the genet was recorded at t-1, but new ramet at t
# GR - genet recruit, the entire genet has no prior data, all RR are GR but not all GR are RR
# NA - first census, or not found
# NO - not a recruit (repeat measure)

# Add in the row
Raw_Stems <- lapply(Raw_Stems, function(y){
  y = y %>% mutate(RecruitStatus = NA)
  
})

# Assign the NO s, they should have measurements earlier, and not be missing later
a <- Raw_Stems[[1]]$StemID[Raw_Stems[[1]]$RametStatus != "P"] # stems that were "present" at t1
Raw_Stems[[2]]$RecruitStatus[
  which((Raw_Stems[[2]]$StemID %in% a) &
          (Raw_Stems[[2]]$RametStatus != "M"))] <- "NO" # assign these to not be considered recrts

# Assign the RR, should have had Ramet not recorded at t1 and have measurements at t2
a <- Raw_Stems[[1]]$StemID[Raw_Stems[[1]]$RametStatus == "P"] 
Raw_Stems[[2]]$RecruitStatus[which((Raw_Stems[[2]]$StemID %in% a))] <- "RR"

# Assign GR, note that all GR are RR, but not all RR are GR.
# Main difference is that the entire genet had to not be present previously
a <- Raw_Stems[[1]]$StemID[Raw_Stems[[1]]$GenetStatus == "P"] 
Raw_Stems[[2]]$RecruitStatus[which((Raw_Stems[[2]]$StemID %in% a))] <- "GR"

# --- QC V2.0 -------------------------------------------------------------------------------------
# I need to go through and ensure that all stems in a single df have only one census number.
# Right now some of them have NA in there.
Raw_Stems %<>% lapply(., function(x){
  x %<>% mutate(Census = ignore_NA( unique(.$Census)))
  return(x)
})

# Run the QC check
QC_res <- QC_Check(Raw_Stems, do.spp_chck = T, do.chng_HOM = T, do.zomb_ram = T)
Raw_Stems <- QC_res$QCdata

# Assign whether a stem is primary or secondary 
Raw_Stems %<>% 
  stemRankAsgn(data_list = ., col_dbh = "DBH", col_out = "Stem",
               consider_if = "A", col_con = "RametStatus", col_grp = "TreeID",
               col_ind = "StemID")

# add a column denoting if the stem was recorded (all were in HVDF, this is for consistency)
Raw_Stems %<>% lapply(., function(x){
   x %<>% mutate(., Recorded = 1)
  return(x)
})

# ### Export the data back out for use in the actual analysis #####################################
AllStems <- Raw_Stems
HVDF_AllStems <- AllStems

# Save HVDF Stemdata
save(HVDF_AllStems,
     file = paste0(loc_Gdr, "/data/ForestGEO/processed/standardized/HVDF_AllStemsOLD.r"))

# Save the Taxatable
save(Spp_table,
     file = paste0(loc_Gdr, "/data/ForestGEO/taxa/HVDF_TaxaTable.r"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # --- Adding a section for QC, changing species ID, Zombie trees and change in POM ----------------
# # We will start with creating a "wide" table to more easily compare stems over time
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
#   i <- i[-which(is.na(i))]
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
#   Raw_Stems[[1]]$GenetStatus == "D"
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
# # leading or trailing semi-colons
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


