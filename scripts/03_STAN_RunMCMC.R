# File: 03_STAN_RunMCMC.R
# Author: JMR
# Last Update: 12023-07-02
# ### Setup ###################################################################
t_ini <- Sys.time() # Timer for the entire script
set.seed(123)

# location of the data
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")

# get directory of the script itself (WARNING MAY NOT WORK FOR BATCH JOBS)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr",
           "rstudioapi", "tidyverse",
           "cowplot")){
  library(x, character.only = T)
}

# STAN options/defaults for the script
options(mc.cores = parallel::detectCores() - 2) # cores JMR mchn w/ 12
rstan_options(auto_write = TRUE,
              javascript = FALSE) # Don't recompile the code if unmodified

# Loading in diagnostics written by M. Betancourt
source('scripts/STAN/stan_utility_rstan.R')

# Place to store the data used for fitting the models, (useful for PPC)
All_Stan_data.L <- list()

# store all the stan models in one place.
All_STAN_fits.L <- list()


# --- Site specific setup -----------------------------------------------------
# Load in the site's data
site <- "HVDF"
load(paste0(loc_data, site, "_STANdata.r"))
#STANdata.L <- get(paste0(site, "_STANModeldata.L"))

# --- Growth ------------------------------------------------------------------
# File of the model
GrowthSTAN_file <- "scripts/STAN/growth_lognormal_estimate.stan"

# Check for the STAN model's file
if(!stanc(GrowthSTAN_file)$status){
  warning("Growth file for STAN script does not exist!")
}

# Check how many intervals are available for the site (n census - 1)
if(STANdata.L$n_intervals !=1){
  warning(paste0("For ", site, " there are ",
                 STANdata.L$n_intervals,
                 " intervals. Currently only estimating values ",
                 "from the first interval"))
}

# Growth data for estimation
growth_all.df <- STANdata.L$growth_data[[1]] # only look at first 2 censuses
n_Lvls <- unique(growth_all.df$CanopyLvl) %>% length # number of canopy Lvls 

# repeat the analysis for each Canopy Level
# create a list output for the growth fits
growth_fitsAllLvl.L <- list(NULL) 

# Parameters for running the model
warm_its <- 1000  # number of warmups
smpl_its <- 1500  # number of USED iterations, (MCMC will run for smpl + warm)
n_coreuse <- 5    # number of cores to use
n_chains <- 5     # number of chains to run

# Loop through for growth
# create a table for storing which Species numbers correspond to which mnemonic
# for after the MCMC
Growth_Species_ID.df <- data.frame()

for(Lvl in 1:n_Lvls){
  # Only get values for this canopy level
  growth_CL.df <- growth_all.df %>%
    filter(CanopyLvl == Lvl)
  
  # temporary dataframe that lists the Mnemonic along with the ID number,(STAN)
  # only deals with integers, this will help us retrieve these values later
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
  
  growth_CL.df %<>% left_join(., temp_species.df, by = "Mnemonic")
  
  # Lsit holding the data for the STAN script
  STAN_growthData.L <- list(n_stem = nrow(growth_CL.df),
                            n_spp = nrow(temp_species.df),
                            species = growth_CL.df$tempID,
                            loggrowth_vals = log(growth_CL.df$growth2))
  
  
  
  t1 <- Sys.time() 
  print(t1)
  # Run the Stan code
  fit1 <- stan(file = GrowthSTAN_file,
               data = STAN_growthData.L,
               warmup = warm_its, 
               iter = smpl_its + warm_its,
               chains = n_chains,
               cores = n_coreuse,
               thin = 1,
               seed = 123,
               control = list(adapt_delta = 0.925))
  
  
  print(round(Sys.time() - t1, 2))
  
  # Run the diagnostics for the chain
  cat("Running dianostics:\n")
  check_all_expectand_diagnostics(fit1)
  check_all_hmc_diagnostics(fit1)
  check_energy(fit1)
  
  # Transfer the results to the list
  growth_fitsAllLvl.L[[Lvl]] <- fit1
  
  All_Stan_data.L <- append(All_Stan_data.L, list(STAN_growthData.L))
  names(All_Stan_data.L)[length(All_Stan_data.L)] <- paste0("growthCont_CL", Lvl)
  
  # Adding to the fits
  All_STAN_fits.L <- append(All_STAN_fits.L, list(fit1))
  names(All_STAN_fits.L)[length(All_STAN_fits.L)] <- paste0("growthCont_CL",Lvl)
  
}

# --- Growth Binary Variable --------------------------------------------------
STAN_file <- "scripts/STAN/growthBinary.stan"
# Check for the STAN model's file
if(!stanc(STAN_file)$status){
  warning("Growth file for STAN script does not exist!")
}

# Binary growth data extracted
print(paste0("Currently only getting the FIRST set of",
      " censuses for binary-growth parameters"))
growthBinary_All.df <- STANdata.L$growth_bin[[1]] 

# Set the number of canopy levels to deal with
n_Lvls <- length(unique(growthBinary_All.df$CanopyLvl))

# Create a dataframe to hold the speciesID for going into STAN
GrowthBinary_Species_ID.df <- data.frame()

growthBin_fitsAllLvl.L <- list(NULL) 

# set the parameters for running the MCMC
warm_its <- 1000  # number of warmups
smpl_its <- 1500  # number of USED iterations, (MCMC will run for smpl + warm)
n_coreuse <- 5    # number of cores to use
n_chains <- 5     # number of chains to run


# loop though the levels
for(Lvl in 1:n_Lvls){
  cat(paste0("Beginning Binary-Growth parameter MCMC on CL: ", Lvl, "\n"))
  
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
  
  growthBinary_CL.df %<>% left_join(., temp_species.df, by = "Mnemonic") %>%
    mutate(grow_Bin = as.integer(grow_Bin))
    
  # Create the list holding the data to send to STAN
  STANModeldata.L <- list(n_stem = nrow(growthBinary_CL.df),
                     n_spp = max(growthBinary_CL.df$tempID),
                     species = growthBinary_CL.df$tempID,
                     grow_bindat = growthBinary_CL.df$grow_Bin,
                     t_interv = growthBinary_CL.df$Date_dur)
  
  # Keeping track of time for the runs
  t1 <- Sys.time() 
  cat(paste0("1", as.character(t1),
             " - Beginning STAN MCMC call for Binary Growth Variable CL-",
             Lvl, "\n"))
  
  # fit the model
  fit1 <- stan(file = STAN_file,
               data = STANModeldata.L,
               warmup = warm_its, 
               iter = smpl_its + warm_its,
               chains = n_chains,
               cores = n_coreuse,
               thin = 1,
               seed = 123,
               control = list(adapt_delta = 0.925))
  
  # End of MCMC
  t2 <- Sys.time() 
  
  cat(paste0("1", as.character(t2),
             " - Finished call.\n"))
  
  print(round(t2 - t1, 2))
  
  # Run the diagnostics for the chain
  cat(paste0("1", as.character(Sys.time()), "- Begin running dianostics:\n"))
  check_all_expectand_diagnostics(fit1)
  check_all_hmc_diagnostics(fit1)
  check_energy(fit1)
  
  # Transfer the results to the list
  growthBin_fitsAllLvl.L[[Lvl]] <- fit1
  
  # Save the data used
  All_Stan_data.L <- append(All_Stan_data.L, list(STANModeldata.L))
  names(All_Stan_data.L)[length(All_Stan_data.L)] <- paste0("growthBinary_CL", Lvl)
  
  # Adding to the fits
  All_STAN_fits.L <- append(All_STAN_fits.L, list(fit1))
  names(All_STAN_fits.L)[length(All_STAN_fi=ts.L)] <- paste0("growthBinary_CL",Lvl)
}


# --- Survival setup ----------------------------------------------------------
# survival_All.df <- STANdata.L$survival[[1]]
# 
# # ### FOR TREATING RAMET PRIMARY STEMS ONLY
# # Process the requirements for the Survival:
# # Survival at the level of the ramet, but only considering the main stem
# # cat("For survival looking at survival of the ramet, only main/single stems.")
# # survival_All.df %<>%
# #  filter(Stem != "Secondary" | is.na(Stem)) %>%
# #  filter(!is.na(ram_survived)) # HERE
# 
# # Process the requirements for the Survival:
# # Survival at the level of the ramet, but only considering the main stem
# cat(paste0("For survival looking at survival of the whole 'genet' ",
#            "Using the maximum Canopy Level assigned to the stem. Other",
#            " options to explore are using median weighted by basal area.\n"))
# survival_All.df %<>%
#   filter(!is.na(gen_survived)) %>%
#   mutate(BA_bef = pi*(DBH_bef/2)^2,
#     w = case_when(is.na(BA_bef) ~ 0,
#                   !is.na(BA_bef) ~ BA_bef),
#     gen_survived = as.integer(gen_survived),
#     ram_survived = as.integer(ram_survived))
# 
# 
# survival_All.df %<>%
#   group_by(TreeID) %>%
#   summarize(Mnemonic = unique(Mnemonic), # Name
#             CanopyLvl = case_when( # take the maximum CL (weight median try?)
#               all(is.na(CanopyLvl)) ~ NA_real_,
#               !all(is.na(CanopyLvl)) ~
#                 suppressWarnings(max(CanopyLvl,na.rm = T)),
#               TRUE ~ -999999999999999),
#             gen_survived = unique(gen_survived), # binary of whether gen surv
#             prop_surv_w = sum(w*ram_survived, # proportion of surv basal Area
#                               na.rm = T)/
#               sum(w, na.rm = T),
#             n_stem = n(), # number of stems
#             Census = unique(Census), # census
#             Date_dur = weighted.mean(Date_dur, # take date by weighing by BA
#                                      w = w,
#                                      na.rm = T)) %>% 
#   filter(!is.na(CanopyLvl))
# 
# # Weight proportional survival by the average duration in measurement
# #
# k <- 6 # round everything to be within 10^k
# survival_All.df %<>% # rescale proportion loss to 1 year
#   mutate(prop_surv_w = round(prop_surv_w^(1/Date_dur), k),
#          prop_surv_w = case_when(prop_surv_w == 1 ~ 1 - 10^(-k),
#                                  prop_surv_w == 0 ~ 0 + 10^(-k),
#                                  TRUE ~ prop_surv_w))
# 
# 
# # getting -Inf values for alive genet ID values in census 2 but had no CL
# # assigned to them, I am going to remove them.
# 
# 
# 
# 
# # Set the number of canopy levels to deal with
# n_Lvls <- length(unique(survival_All.df$CanopyLvl))

# ---- binary Survival --------------------------------------------------------
STAN_file <- "scripts/STAN/Bsurvival.stan"
# Check for the STAN model's file
if(!stanc(STAN_file)$status){
  warning(paste0('"', STAN_file,'"', " does not exist!"))
}

print(paste0("Currently only getting the FIRST set of",
             " censuses for survival parameters"))

# Create a data frame to hold the speciesID for going into STAN
Bsurvival_Species_ID.df <- data.frame()

Bsurvival_fitsAllLvl.L <- list(NULL) 

# set the parameters for running the MCMC
warm_its <- 1000  # number of warmups
smpl_its <- 1500  # number of USED iterations, (MCMC will run for smpl + warm)
n_coreuse <- 5    # number of cores to use
n_chains <- 5     # number of chains to run

# loop though the levels
for(Lvl in 1:n_Lvls){
  cat(paste0("Beginning survival parameter MCMC on CL: ", Lvl, "\n"))
  
  # Get the data for just the current CL
  survival_CL.df <- survival_All.df %>%
    filter(CanopyLvl == Lvl)
  
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
  
  survival_CL.df %<>% left_join(., temp_species.df, by = "Mnemonic") %>%
    mutate(survived = as.integer(gen_survived)) # change to ram_surv. if desird
  
  # Create the list holding the data
  STANModeldata.L <- list(n_stem = nrow(survival_CL.df),
                     n_spp = max(survival_CL.df$tempID),
                     species = survival_CL.df$tempID,
                     surv_bindat = survival_CL.df$survived,
                     t_interv = survival_CL.df$Date_dur
                     )
  
  # Keeping track of time for the runs
  t1 <- Sys.time() 
  cat(paste0("1", as.character(t1),
             " - Beginning STAN MCMC call for Survival CL-",
             Lvl, "\n"))
  
  # fit the model
  fit1 <- stan(file = STAN_file,
               data = STANModeldata.L,
               warmup = warm_its, 
               iter = smpl_its + warm_its,
               chains = n_chains,
               cores = n_coreuse,
               thin = 1,
               seed = 123,
               control = list(adapt_delta = 0.90))
  
  # End of MCMC
  t2 <- Sys.time() 
  
  cat(paste0("1", as.character(t2),
             " - Finished call.\n"))
  
  print(round(t2 - t1, 2))
  
  # Run the diagnostics for the chain
  cat(paste0("1", as.character(Sys.time()), "- Begin running dianostics:\n"))
  check_all_expectand_diagnostics(fit1)
  check_all_hmc_diagnostics(fit1)
  check_energy(fit1)
  
  # Transfer the results to the list
  Bsurvival_fitsAllLvl.L[[Lvl]] <- fit1
  
  # Save the data used
  All_Stan_data.L <- append(All_Stan_data.L, list(STANModeldata.L))
  names(All_Stan_data.L)[length(All_Stan_data.L)] <- paste0("BinarySurv_CL", Lvl)
  
  # Adding to the fits
  All_STAN_fits.L <- append(All_STAN_fits.L, list(fit1))
  names(All_STAN_fits.L)[length(All_STAN_fits.L)] <- paste0("BinarySurv_CL",Lvl)
  
}

# ---- proportional Survival --------------------------------------------------
# STAN_file <- "scripts/STAN/Psurvival.stan"
# # Check for the STAN model's file
# if(!stanc(STAN_file)$status){
#   warning(paste0('"', STAN_file,'"', " does not exist!"))
# }
# 
# print(paste0("Currently only getting the FIRST set of",
#              " censuses for survival parameters"))
# # Create a data frame to hold the speciesID for going into STAN
# Psurvival_Species_ID.df <- data.frame()
# 
# Psurvival_fitsAllLvl.L <- list(NULL)
# 
# # set the parameters for running the MCMC
# warm_its <- 1000  # number of warmups
# smpl_its <- 1500  # number of USED iterations, (MCMC will run for smpl + warm)
# n_coreuse <- 5    # number of cores to use
# n_chains <- 5     # number of chains to run
# 
# # loop though the levels
# for(Lvl in 1:n_Lvls){
#   cat(paste0("Beginning P-survival parameter MCMC on CL: ", Lvl, "\n"))
# 
#   # Get the data for just the current CL
#   survival_CL.df <- survival_All.df %>%
#     filter(CanopyLvl == Lvl)
# 
#   # Create the temporary space to manipulate the speceis ID
#   temp_species.df <- data.frame(Mnemonic = unique(survival_CL.df$Mnemonic),
#                                 tempID = 1:length(unique(
#                                   survival_CL.df$Mnemonic)))
# 
#   # Include the number of observations for each
#   temp_species.df <- survival_CL.df %>%
#     group_by(Mnemonic) %>%
#     summarize(n_obs = n()) %>%
#     left_join(temp_species.df, ., by = "Mnemonic")
# 
#   # save the tempID so we can retrieve individula taxa later
#   Psurvival_Species_ID.df <- temp_species.df %>%
#     mutate(CanopyLvl = Lvl) %>%
#     bind_rows(Psurvival_Species_ID.df, .)
# 
#   survival_CL.df %<>% left_join(., temp_species.df, by = "Mnemonic")
# 
#   # Create the list holding the data
#   STANModeldata.L <- list(n_obs = nrow(survival_CL.df),
#                      n_OTU = length(unique(survival_CL.df$tempID)),
#                      OTU = survival_CL.df$tempID,
#                      propn = survival_CL.df$prop_surv_w)
# 
#   # Keeping track of time for the runs
#   t1 <- Sys.time()
#   cat(paste0("1", as.character(t1),
#              " - Beginning STAN MCMC call for P-Survival CL-",
#              Lvl, "\n"))
# 
#   # fit the model
#   fit1 <- stan(file = STAN_file,
#                data = STANModeldata.L,
#                warmup = warm_its,
#                iter = smpl_its + warm_its,
#                chains = 1,
#                cores = n_coreuse,
#                thin = 1,
#                seed = 123,
#                control = list(adapt_delta = 0.90))
# 
#   # End of MCMC
#   t2 <- Sys.time()
# 
#   cat(paste0("1", as.character(t2),
#              " - Finished call.\n"))
# 
#   print(round(t2 - t1, 2))
# 
#   # Run the diagnostics for the chain
#   cat(paste0("1", as.character(Sys.time()), "- Begin running dianostics:\n"))
#   check_all_expectand_diagnostics(fit1)
#   check_all_hmc_diagnostics(fit1)
#   check_energy(fit1)
# 
#   # Transfer the results to the list
#   Psurvival_fitsAllLvl.L[[Lvl]] <- fit1
# 
# }
# 

# ### Output Results ##########################################################
# Group together the SpeciesIDs
All_params_SpeciesID.df <- Growth_Species_ID.df %>%
  mutate(param = paste0("growthCont_CL",CanopyLvl))
# add the others
All_params_SpeciesID.df <- GrowthBinary_Species_ID.df %>%
  mutate(param = paste0("growthBinary_CL",CanopyLvl)) %>%
  bind_rows(All_params_SpeciesID.df, .)

All_params_SpeciesID.df <- Bsurvival_Species_ID.df %>%
  mutate(param = paste0("BinarySurv_CL",CanopyLvl)) %>%
  bind_rows(All_params_SpeciesID.df, .)

# Save the output MCMC so we don't have to keep re-running it
save(list = c("All_STAN_fits.L",
              "All_Stan_data.L",
              "All_params_SpeciesID.df"),
     file = paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/", site,
                   "_STAN_MCMC.r"))

# Timer for the entire script
t_fnl <- Sys.time()

cat(paste0("All runs completed.\n", 
           t_ini, " - Began\n",
           t_fnl, " - Ended\n"))

t_fnl - t_ini
