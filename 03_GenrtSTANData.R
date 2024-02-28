# 03_GenrtSTANData
# Generate STANData_lists

# ### Setup ###################################################################
if(!file.exists("./scripts/AA_LocationManager.R")){setwd("..")}
source("./scripts/AA_LocationManager.R")
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")

# Bringing in the packages:
for(x in c("magrittr", "rstudioapi", "tidyverse")){
  library(x, character.only = T)
}

site <- "HVDF"
for(site in c("HVDF", "SCBI", "SERC")){


load(paste0(loc_data, site, "_STANdata.r"))

# ### Parameters ##############################################################
# Empty List of STAN objects
All_STAN_data.L <- list()

# Empty dataframe of parameters and their corresponding n_obs, TempID 
All_params_SpeciesID.df <- data.frame(NULL)

# --- Continuous growth -------------------------------------------------------
# Growth data for estimation
growth_all.df <- STANdata.L$growth_data[[1]] # only look at first 2 censuses
n_Lvls <- unique(growth_all.df$CanopyLvl) %>% length # number of canopy Lvls 

# Loop through for growth
# create a table for storing which Species numbr correspond to which mnemonic
# for after the MCMC
Growth_Species_ID.df <- data.frame()

for(Lvl in 1:n_Lvls){
  # Only get values for this canopy level
  growth_CL.df <- growth_all.df %>%
    filter(CanopyLvl == Lvl)
  
  # temporary dataframe that lists the Mnemonic with the ID number,(STAN)
  # only deals with integers, this will help us get these later
  temp_species.df <- data.frame(Mnemonic = unique(growth_CL.df$Mnemonic),
                                tempID = 1:length(unique(
                                  growth_CL.df$Mnemonic)))
  
  # Include the number of observations for each 
  temp_species.df <- growth_CL.df %>%
    group_by(Mnemonic) %>%
    summarize(n_obs = n()) %>%
    left_join(temp_species.df, ., by = "Mnemonic")
  
  # Save the temporary numeric ID for each OTU so it can be extracted later
  Growth_Species_ID.df <- temp_species.df %>%
    mutate(CanopyLvl = Lvl) %>%
    bind_rows(Growth_Species_ID.df, .)
  
  growth_CL.df %<>%
    left_join(., temp_species.df, by = "Mnemonic")
  
  # Lsit holding the data for the STAN script
  STAN_growthData.L <- list(n_stem = nrow(growth_CL.df),
                            n_spp = nrow(temp_species.df),
                            species = growth_CL.df$tempID,
                            loggrowth_vals = log(growth_CL.df$growth2))

  # Add it to the list of the STAN data lists, along with modifying the name
  All_STAN_data.L <- append(All_STAN_data.L,
                            list(STAN_growthData.L))
  # rename the list item so it is can be retrieved later
  names(All_STAN_data.L)[
    length(All_STAN_data.L)] <- paste0("growthCont_CL", Lvl)
  
  
}

# update the chart showing how Mnemonic relates to TempID
All_params_SpeciesID.df <- Growth_Species_ID.df %>%
  mutate(param = paste0("growthCont_CL",CanopyLvl)) %>%
  bind_rows(All_params_SpeciesID.df, .)

# --- Growth Binary -----------------------------------------------------------
# Binary growth data extracted
print(paste0("Currently only getting the FIRST set of",
             " censuses for binary-growth parameters"))
growthBinary_All.df <- STANdata.L$growth_bin[[1]] 

# Set the number of canopy levels to deal with
n_Lvls <- length(unique(growthBinary_All.df$CanopyLvl))

# Create a dataframe to hold the speciesID for going into STAN
GrowthBinary_Species_ID.df <- data.frame()

# loop though the levels
for(Lvl in 1:n_Lvls){
  # Get the data for just the current CL
  growthBinary_CL.df <- growthBinary_All.df %>%
    filter(CanopyLvl == Lvl)
  
  # Create the temporary space to manipulate the speceis ID
  temp_species.df <- data.frame(Mnemonic = unique(growthBinary_CL.df$Mnemonic),
                                tempID = 1:length(unique(
                                  growthBinary_CL.df$Mnemonic)))
  
  # Include the number of observations for each 
  temp_species.df <- growthBinary_CL.df %>%
    group_by(Mnemonic) %>%
    summarize(n_obs = n()) %>%
    left_join(temp_species.df, ., by = "Mnemonic")
  
  # save the tempID so we can retrieve individula taxa later
  GrowthBinary_Species_ID.df <- temp_species.df %>%
    mutate(CanopyLvl = Lvl) %>%
    bind_rows(GrowthBinary_Species_ID.df, .)
  
  growthBinary_CL.df %<>%
    left_join(., temp_species.df, by = "Mnemonic") %>%
    mutate(grow_Bin = as.integer(grow_Bin))
  
  # Create the list holding the data
  STANModeldata.L <- list(n_stem = nrow(growthBinary_CL.df),
                          n_spp = max(growthBinary_CL.df$tempID),
                          species = growthBinary_CL.df$tempID,
                          bindat = growthBinary_CL.df$grow_Bin,
                          t_interv = growthBinary_CL.df$Date_dur)
  

  # Save the data used
  All_STAN_data.L <- append(All_STAN_data.L, list(STANModeldata.L))
  names(All_STAN_data.L)[
    length(All_STAN_data.L)] <- paste0("growthBinary_CL", Lvl)

} # End of loop through CL

# update the chart showing how Mnemonic relates to TempID
All_params_SpeciesID.df <- GrowthBinary_Species_ID.df %>%
  mutate(param = paste0("growthBinary_CL", CanopyLvl)) %>%
  bind_rows(All_params_SpeciesID.df, .)

# --- survival binary ---------------------------------------------------------

cat("Currently only getting the FIRST set of",
    "censuses for survival parameters.\n")

survival_All.df <- STANdata.L$survival[[1]]

# Create a data frame to hold the speciesID for going into STAN
Bsurvival_Species_ID.df <- data.frame()

# loop though the levels
for(Lvl in 1:n_Lvls){
  # Get the data for just the current CL
  survival_CL.df <- survival_All.df %>%
    filter(CanopyLvl == Lvl) %>%
    filter(is.na(Stem) | Stem == "Main")
  
  # Create the temporary space to manipulate the speceis ID
  temp_species.df <- data.frame(Mnemonic = unique(survival_CL.df$Mnemonic),
                                tempID = 1:length(unique(
                                  survival_CL.df$Mnemonic)))
  
  # Include the number of observations for each 
  temp_species.df <- survival_CL.df %>%
    group_by(Mnemonic) %>%
    summarize(n_obs = n()) %>%
    left_join(temp_species.df, ., by = "Mnemonic")
  
  # save the tempID so we can retrieve individula taxa later
  Bsurvival_Species_ID.df <- temp_species.df %>%
    mutate(CanopyLvl = Lvl) %>%
    bind_rows(Bsurvival_Species_ID.df, .)
  
  survival_CL.df %<>%
    left_join(., temp_species.df, by = "Mnemonic") %>%
    mutate(survived = as.integer(ram_survived)) # change to ram_surv?
  
  # Create the list holding the data
  STANModeldata.L <- list(n_stem = nrow(survival_CL.df),
                          n_spp = max(survival_CL.df$tempID),
                          species = survival_CL.df$tempID,
                          bindat = survival_CL.df$survived,
                          t_interv = survival_CL.df$Date_dur
  )
  
  # Save the data used
  All_STAN_data.L <- append(All_STAN_data.L, list(STANModeldata.L))
  names(All_STAN_data.L)[
    length(All_STAN_data.L)] <- paste0("BinarySurv_CL", Lvl)
  
}# End of loop through CLs

# update the chart showing how Mnemonic relates to TempID
All_params_SpeciesID.df <- Bsurvival_Species_ID.df %>%
  mutate(param = paste0("BinarySurv_CL",CanopyLvl)) %>%
  bind_rows(All_params_SpeciesID.df, .)

# ### Output the files ########################################################
save(list = c("All_STAN_data.L",
              "All_params_SpeciesID.df"),
     file = paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/", site,
                   "_STANData_MCMC.r"))
}

# End of script