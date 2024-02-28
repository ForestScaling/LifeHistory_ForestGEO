# 00_ReformatStdz.R
# Author: JMR
# convert from the prev. type of data format to the new one (SCBI, HVDF, SERC)
# ### Setup ###################################################################
# File locations - will work to stdize loc with a sep script
if(!file.exists("./scripts/AA_LocationManager.R")){setwd("..")}
source("./scripts/AA_LocationManager.R")
loc_1 <- "/data/ForestGEO/processed/standardized/"

# General options
min_DBH <- 10 # minimum size in cm
filter_DBH <- T # Whether or not to filter trees by a minimum DBH

# Load in the packages
for(i in c("tidyverse", "magrittr")){
  library(i, character.only = T)
}

# Load common functions
source("./scripts/00_Functions.R")

# ### Reformat ################################################################
# Load in the file (note that this will overwrite the prev. standardized one)
site.v <- c("HVDF", "SCBI", "SERC") # sites to reformat

# loop through the sites
for(site in site.v){

# Load the data
load(paste0(loc_Gdr, loc_1, paste0(site, "_AllStemsOLD.r")))
AllStems <- get(paste0(site, "_AllStems"))
# going forward the loaded data won't store the object with the site code,
# this is cumbersome.




# --- Adding artificial cutoffs into the data ---------------------------------
# rerunning quality control for stems so that artificial cutoff in DBH can be 
# added
Raw_Stems <- AllStems

if(filter_DBH){
# Trimming the trees smaller than minimum DBH
  cat("Filtering out observations where tree DBH was less than ",
      min_DBH,
      "cm ...",
      sep = "")
  
  Raw_Stems <- lapply(Raw_Stems, FUN = function(X){
    X <- X %>%
      filter(., DBH >= min_DBH)
  })
}

#Raw_Stems <- lapply(Raw_Stems, function(x){return(mutate(x, RecruitStatus = "NO"))})
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


# Now we need to add data about the genet
# A - genet alive, some part of the tree/shrub is alive
# D - the entire tree/shrub is dead
# M - missing, no data for the individual, but was present in a previous survey
# P - Prior, individual was not present in this survey but is present in future survey(s)
# F - category which is a subset of A where the stem(s) alive are <10mm DBH, IS NOT BEING INCLUDED
#     RIGHT NOW!
Raw_Stems %<>% lapply(., FUN = function(X){
  X <- X %>% mutate(Codes = "")
  return(X)
})
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


AllStems <- Raw_Stems

# --- Add commonality in coordinates (also check for changes) -----------------
# Need to make it so locations are set to the stem at all observations

# get the unique positions for X and Y, it should fail if there are 
# multiple values for a single stem (not including NA)
PX_uniq.v <- AllStems %>%
  # Extract the PX column
  sapply(., FUN = function(X){
    return(dplyr::pull(X, PX))}) %>% 
  # take the unique value (ignore NA)
  apply(MARGIN = 1, FUN = function(Y){
    if(all(is.na(Y))){
      out <- NA
    }else{
      out <- Y[which(!is.na(Y))]
    }
    return(unique(out))
  })


# repeat for Y
PY_uniq.v <- AllStems %>%
  # Extract the PY column
  sapply(., FUN = function(X){
    return(dplyr::pull(X, PY))}) %>% 
  # take the unique value (ignore NA)
  apply(MARGIN = 1, FUN = function(Y){
    if(all(is.na(Y))){
      out <- NA
    }else{
      out <- Y[which(!is.na(Y))]
    }
    return(unique(out))
  }) 


# Verify that only one PX and PY value per stem (exclude NA)
if(!(is.atomic(PY_uniq.v) && is.atomic(PX_uniq.v))){
  stop(paste("Multiple coordinate (PX and/or PY, for a single",
             "stem over the censuses (ignoring NA)! Please fix."))
} 


# Now fill in the values
AllStems <- lapply(AllStems, FUN = function(X){
  X %>% mutate(PX = PX_uniq.v, PY = PY_uniq.v) %>%
    return(.)
})

# --- Add column to easily flag change in HOM -----------------------
AllStems <- AllStems %>% lapply(X = ., FUN = function(X){
  X %<>% mutate(Chng_HOM = as.integer(CodesLogic(Codes,
                                                 sep = ";",
                                                 c = "A")))
  return(X)
})



# --- Filter out the random columns -------------------------------------------
# columns that are useful.
cols <- c("StemID", "TreeID", "Mnemonic", "PX", "PY", 
          "Census", "DBH","Date", "Stem", "HOM", "Chng_HOM", "Codes",
          "RametStatus", "GenetStatus", "RecruitStatus")

# Only keep necessary columns
AllStems <- lapply(AllStems, FUN = function(X){
  X <- X %>% select(all_of(cols))
  return(X)
})


# --- Write outputs -----------------------------------------------------------
save(AllStems, file = paste0(loc_Gdr, loc_1, paste0(site, "_AllStems.r")))
}

# End of script
