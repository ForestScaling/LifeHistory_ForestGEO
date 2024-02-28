# Script: 00_DataInitialClean_SCBI.R
# Author: JMR 
# Last Update: 12022-01-17
# Description:
# File in-taking the data for SCBI and spitting out a standardized version so that it is harmonized
# with the other ForesGEO sites. Outputs two files, a "taxa" file, containing the information to 
# correlate stems to actual species of interest and (if you want to drop particular problematic 
# lianas for instance) as well as the AllStem file, which contains all the stems measured across 
# the censuses.
# ### Setup #######################################################################################
# Modifications to specific datasets, before initial processing
# Verify Locations
if(!file.exists("./scripts/AA_LocationManager.R")){setwd("..")}
source("./scripts/AA_LocationManager.R")

# Setting/finding location of the script so the WD is automatically within GITHub local copy loc.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

# Bringing in the packages:
library("tidyverse")
library("magrittr")

# Get the commonly used functions
source("scripts/00_Functions.R")

# Site Name:
site <- "SCBI"

# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/",
                                site))[
  str_which(list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/",
                              site)), pattern = "_AllMeas_")]

# Create a file for storing all of the Stem files
Raw_Stems <- lapply(stem_files, FUN = reader1, loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))

# Read in the species table
Spp_table <- read.delim(paste0(loc_Gdr, "/data/ForestGEO/raw/", site,"/",
                               site,"_TaxaGenusSp_Report.txt"), sep =  "\t", header = T)

# ### Harmonize Taxonomic table ###################################################################
# First verifying that all species in the stem data actually have a taxa entry
if(length(which(sapply(Raw_Stems, FUN = function(z){
  length(which(!(z$Mnemonic %in% Spp_table$Mnemonic)))}) == 0)) != length(Raw_Stems)
)warning(paste0("At least one specieces identification code in ",
                "the Stem data is NOT represented in the table of taxa!"))
# look's like there is a discrepancy

#  Fix column names and add the column for the latin name
Spp_table %<>% rename(SubSpecies = Subspecies) %>%
  mutate(Latin = paste(Genus, Species))

# We need to address this!
NoMnem <- unique(unlist(sapply(Raw_Stems, function(z){
  return(z$Latin[which(!(z$Mnemonic %in% Spp_table$Mnemonic))])
})))

# Message indicating the species not represented in the table:
cat(paste0("The following Latin names in the stems file aren't in the taxonomic table:\n"))
print(NoMnem)

# Now we need to add into the spp table those taxa that are not represented
NewMnem <- distinct(do.call(rbind, lapply(Raw_Stems, function(x){
  x[which(x$Latin %in% NoMnem),c(2:4)]
}))) #Get a set of not listed ones

# Create a new df with the taxa to be added
Spp_table <- NewMnem %>%
  mutate(Genus = str_split_fixed(.$Latin, " ", 2)[, 1],
         Species = str_split_fixed(.$Latin, " ", 2)[, 2],
         IDlevel = case_when(Species == "sp" ~ "multiple",
                             TRUE ~ "species"),
         Family = case_when(Genus == "Prunus" ~ "Rosaceae",
                            Genus == "Euonymus" ~ "Celastraceae",
                            Genus == "Viburnum" ~ "Adoxaceae",
                            Genus == "Rhododendron" ~ "Ericaceae",
                            TRUE ~ "ERROR_PLEASE_VERIFY!")) %>% bind_rows(Spp_table, .)

# A few other adjustments
Spp_table %<>%
  mutate(Species = case_when(Species == "sp" ~ "sp.", TRUE ~ Species),
         Latin = paste(Genus, Species),
         DropBegin = case_when(Family  == "Vitaceae" | Genus == "Toxicodendron" ~ 1,
                               TRUE ~ 0),
         DropParamtrz = case_when(IDlevel %in% c("species", "subspecies") ~ 0,
                                  TRUE ~ 1),
         Clade = case_when(Family == "Unknown" ~ "U",
                           Family %in% c("Pinaceae", "Taxaceae", "Cupressaceae") ~ "G",
                           TRUE ~ "A"))

# Rechecking that everyone is now represented in the spp table
if(length(which(sapply(Raw_Stems,FUN = function(z){
  length(which(!(z$Mnemonic %in%
                 Spp_table$Mnemonic)))}) == 0)) == length(Raw_Stems)){
  warning(paste("Assignments for missing taxa in species table were made for:",
                paste(NoMnem,collapse =  ", ")))
}else{
  warning(paste0("At least one specieces identification code in the Stem data is STILL ",
                 "NOT represented in the table of taxa!"))}

# Update taxa names in the Stem Data
Raw_Stems %<>% lapply(., function(x){
  x %<>% select(-c("Latin", "SubSpecies")) %>%
    left_join(select(Spp_table, Mnemonic, Latin, SubSpecies), by = "Mnemonic") %>%
    select(colnames(x))
  return(x)
})

# ### Harmonizing Stem Data #######################################################################
# Apparently, non-multistemmed individuals still have "main" and "secondary" classification schema,
# Main and Secondary should only apply to multistemmed indiv.
Raw_Stems <- lapply(Raw_Stems, function(x){
  df = x %>% group_by(TreeID) %>% summarise(n = n())
  df = df$TreeID[which(df$n == 1)]
  
  x$Stem[which(x$TreeID %in% df)] = ""
  
  return(x)
})

# Convert DBH to cm (for consistency), and switch it so that DBH of 0 is NA and HOM of 0 is NA
Raw_Stems %<>% lapply(., function(y){
  y %<>% mutate(DBH = case_when(is.na(DBH) | DBH == 0 ~ NA_real_,
                                DBH > 0 ~ DBH/10,
                                TRUE ~ -99999999999),
                HOM = case_when(HOM == 0 ~ NA_real_,
                                TRUE ~ HOM))
  
  return(y)
  
})


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

# Set a column up for recruitment
# Values are either:

# RR - ramet is recruit, the genet was recorded at t-1, but new ramet at t
# GR - genet recruited, the entire genet has no measurements prior to t,
# *note that all RR are GR but not all GR are RR
# NA - first census
# NO - not a recruit (repeat measure)
# * Note that if the individual was found but no data recorded (Missing) it still must of been
# present to be "seen" so assumed to have been present at time. so can still be cosidered a recruit
# (just not know if alive or not)
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
# P - Prior, individual was not present in this survey, applicable later

# Adding in the column itself
Raw_Stems %<>% lapply(., function(x){
  x %<>% mutate(RametStatus = case_when(.$Status == "missing" ~ "M",
                                        .$Status %in% c("stem dead", "broken below", "dead") ~ "D",
                                        .$Status == "alive" ~ "A",
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


# ### Export the data back out for use in the actual analysis #####################################
SCBI_AllStems <- Raw_Stems
SCBI_TaxaTable <- Spp_table

# Save SCBI Stemdata
save(SCBI_AllStems,
     file = paste0(loc_Gdr,"/data/ForestGEO/processed/standardized/SCBI_AllStemsOLD.r"))

# Save the Taxatable
save(Spp_table,
     file = paste0(loc_Gdr,"/data/ForestGEO/taxa/SCBI_TaxaTable.r"))


# My graveyard of old stuff, I guess intelligently designed systems also have vestigial characters...

# 
# # Set a column up for recruitment
# # Values are either:
# 
# # RR - ramet is recruit, the genet was recorded at t-1, but new ramet at t
# # GR - genet recruited, the entire genet has no measurements prior to t, note that all RR are GR but not all GR are RR
# # NA - first census, or not found
# # NO - not a recruit (repeat measure)
# Raw_Stems <- lapply(Raw_Stems, function(y){
#   y = y %>% mutate(RecruitStatus = NA)
#   return(y)
# })
# 
# # Add in values, lets start off with the NO-s 
# # for t3, check t2 and t1 if there already is one
# a1 <- Raw_Stems[[1]]$StemID[which(Raw_Stems[[1]]$Status != "missing")]
# a2 <- Raw_Stems[[2]]$StemID[which(Raw_Stems[[2]]$Status != "missing")]
# a3 <- Raw_Stems[[3]]$StemID[which(Raw_Stems[[3]]$Status != "missing")]
# Raw_Stems[[3]]$RecruitStatus[which((Raw_Stems[[3]]$StemID %in% c(a1, a2)))] <- "NO"
# # for t2 check t1
# Raw_Stems[[2]]$RecruitStatus[which((Raw_Stems[[2]]$StemID %in% c(a1)))] <- "NO"
#   
# # Add values for the RR, if they weren't in any previous censuses
# Raw_Stems[[3]]$RecruitStatus[which(!(Raw_Stems[[3]]$StemID %in% c(a1, a2)))] <- "RR"
# Raw_Stems[[2]]$RecruitStatus[which(!(Raw_Stems[[2]]$StemID %in% c(a1)))] <- "RR"
# 
# # Now we need to override the RR that had complete genet recruitment
# a1 <- unique(Raw_Stems[[1]]$TreeID[which(Raw_Stems[[1]]$Status != "missing")])
# a2 <- unique(Raw_Stems[[2]]$TreeID[which(Raw_Stems[[2]]$Status != "missing")])
# Raw_Stems[[3]]$RecruitStatus[which(!(Raw_Stems[[3]]$TreeID %in% c(a1, a2)))] <- "GR"
# Raw_Stems[[2]]$RecruitStatus[which(!(Raw_Stems[[2]]$TreeID %in% c(a1)))] <- "GR"
# 
# 
# # we are going to add a standardized status column 'RametStatus' values are either:
# # A - ramet alive
# # D - ramet dead
# # M - missing, no data for the individual, but was present last time.
# # P - Prior, individual was not present in this survey
# 
# # Add in the column and fill in the values
# Raw_Stems <- lapply(Raw_Stems, function(z){
#   z = mutate(z, RametStatus = Status)
#   z$RametStatus[which(z$RametStatus == "missing")] = "M"
#   z$RametStatus[which(z$RametStatus == "alive")] = "A"
#   z$RametStatus[which(z$RametStatus %in% c("stem dead", "dead", "broken below"))] = "D"
#   return(z)
# })
# 
# df <- lapply(Raw_Stems,function(z){
#   return(z[,c(8,10,13,15,16,17,18,20)])
# })
# 
# df2 <- full_join(df[[1]],df[[2]], by = "StemID") %>% full_join(df[[3]],by = "StemID")
# 
# df2$Order <- paste(df2$RametStatus.x,df2$RametStatus.y,df2$RametStatus,sep ="_")
# 
# View(df2[which(df2$Order %in% c("D_A_A", "A_D_A","D_A_D","D_D_A","D_M_A","NA_D_A")),])
# 
# 
# #lapply(Raw_Stems,function(y){unique(y$Status)})
# 
# 
# 
# 
# # We need to check what sorts of codes are present here:
# codes <- unique(unlist(sapply(Raw_Stems,function(z){unique(str_split(z$Codes,";"))})))
# statuses <- unique(unlist(sapply(X = Raw_Stems,FUN =  function(z){unique(z$Status)})))
# 
# View(Raw_Stems[[1]][which(Raw_Stems[[1]]$Status == statuses[4]),])
# 
# 
# 
# 
# 
# 
# lapply(Raw_Stems,function(z){unique(z$Status[which(CodesLogic(z$Codes,"X"))])})
# 
# 
# # 
# # 
# # Thoughts:
# #   Main vs Secondary Stem:
# #   I am not sure how much I trust assignments in this field (HF has some stems who have the same TreeID but aren’t coded with M or S). If we want to do only main individuals changing over time, we can just look at the biggest diameter individual from t1 and follow it (still follow it ‘till death even if it ends up being smaller than another ramet in the genet), and if a new recruit comes in who has multiple stems, just continue onwards focusing on the biggest one at the time of recruitment. 
# # This also brings into question what other methods could be used:
# #   1.	Mains only (see above)
# # 2.	Ramet at random from genet’s first recruitment (e.g. if a I. verticillata recruits at t2, at random choose one of the stems and have it be the individual, follow it through time even if it no longer is the “main stem”).
# # 3.	Ramet at random but sample as a multinomial with each stem weighed by its diameter (or basal area. E.g. an Ilex verticillate recruiting at t2 with three stems (A 10mm, B 12mm, C 17mm), random number choose A, B or C, but the probability of A is 10/(10+12+17) ~0.26, prob of B is ~0.31, and probability of C is ~0.44.
# #                                                                                        
# 
# 
# 
# 
# 
# 
# 
# x <- lapply(Raw_Stems, function(z){
#   return(z[,c(8,10,13,15,17,18)])})
# 
# df <- left_join(x[[1]],x[[2]],by = "StemID") %>% left_join(x[[3]],by = "StemID")
# lapply()
# 
# # Find all codes
# unique(unlist(lapply(Raw_Stems,function(z){
#   str_split(z$Codes, pattern = ";")
#   
# })))
# 
#wide <- ListWiden(SCBI_AllStems, kcols =c("DBH","Codes", "RametStatus"))


#SCBI_AllStems %>% bind_rows() %>% mutate(ha = paste0(RametStatus, RecruitStatus)) %>% group_by(StemID) %>%
#  summarize(h = paste0(RametStatus, collapse = "-")) %>% group_by(h) %>% summarize(n1 = n())%>%
#  mutate( perc = 100*n1/sum(.$n1))


#SERC_AllStems %>% bind_rows() %>% mutate(ha = paste0(RametStatus, RecruitStatus)) %>% group_by(StemID) %>%
#  summarize(h = paste0(RametStatus, collapse = "-")) %>% group_by(h) %>% summarize(n1 = n())%>%
#  mutate( perc = 100*n1/sum(.$n1))


#HVDF_AllStems %>% bind_rows() %>% mutate(ha = paste0(RametStatus, RecruitStatus)) %>% group_by(StemID) %>%
#  summarize(h = paste0(RametStatus, collapse = "-")) %>% group_by(h) %>% summarize(n1 = n())%>%
#  mutate( perc = 100*n1/sum(.$n1))


