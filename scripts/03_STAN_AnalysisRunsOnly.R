# 03_STAN_AnalysisRunsOnly.R
t_ini <- Sys.time() # Timer for the entire script

options(warn = 1)
# set seed for the runs
set.seed(123)


# location of the data
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/")

# get directory of the script itself (WARNING MAY NOT WORK FOR BATCH JOBS)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr",
           "rstudioapi", "tidyverse",
           "parallel",
           "cowplot")){
  library(x, character.only = T)
}

# STAN options/defaults for the script
options(mc.cores = parallel::detectCores() - 2) # cores JMR mchn w/ 12
rstan_options(auto_write = TRUE,
              javascript = FALSE) # Don't recompile the code if unmodified

# Loading in diagnostics script written by M. Betancourt
cat("Loading MCMC NUTS diagnostics script written by M. Betancourt...\n")
source('scripts/STAN/stan_utility_rstan.R')

# Function to check all relevant diagnostics:
JMR_smry_diagnostics <- function(STANfit, fitname = "Model"){
  cat("Running dianostics for ", fitname,":\n", sep ="")
  check_all_expectand_diagnostics(STANfit)
  check_all_hmc_diagnostics(STANfit)
  check_energy(STANfit)
}

# Default vals for runs (override individually if desired.)
warm_its <- 2000  # number of warmups
smpl_its <- 5000 # number of USED iter, (MCMC will run for smpl + warm)
n_coreuse <- 5    # number of cores to use
n_chains <- 5     # number of chains to run
n_thin <- 1       # thinning, 1 = none
NUTS_adapt_delta <- 0.85 # adapt delta, default is 0.801, larger means
# ~smaller "steps" so less divergence, but less
# efficient

STAN_file.L <- list(growthCont = "growth_lognormal_reparam.stan", # alt: growth_lognormal_reparam.stan
                    growthBinary = "binary_logitparam.stan",
                    BinarySurv = "binary_logitparam.stan"
)


# ### STAN MCMC calls #########################################################
# Begin running through the sites
site <- "HVDF"
for(site in c("HVDF", "SCBI", "SERC")){
load(paste0(loc_data, site, "_STANData_MCMC.r")) # load the data, ready to go.

fit_names <- All_STAN_data.L %>% names(.) # names of all the stand Data sets
#fit_names <- "growthBinary_CL1"
# store all the stan models in one place.
All_STAN_fits.L <- vector("list", length = length(All_STAN_data.L))
names(All_STAN_fits.L) <- fit_names

# --- Begin calls ------------------------------------------------------------
cat("Beginning calls for ", site, " at: ",
    as.character(Sys.time()),
    "\n\n", sep = "")


# loop through the different datasets
for(f in names(All_STAN_data.L)){
  cat("\nStarting: ", f, "\n", sep = "")
  
  
  # extract the type of model (ignore CL for now)
  f_type <- gsub(x = f, pattern = "_CL\\d", replacement = "")
  #if(f_type == "growthCont"){next}
  # switch for discerning which script
  STAN_file <- do.call(switch, args = c(EXPR = f_type, STAN_file.L))
  
  # handle when the dataset doesn't match one of the types and that the
  # stan file works
  if(is.null(STAN_file)){
    cat(f, "doesn't match one of the types for the STAN files.",
        "Continuing to the next dataset.\n\n")
    next
    
  # Verify that the STAN file exists
  }else{
    # add the folder navigation
    STAN_file <- paste0("./scripts/STAN/", STAN_file) 
    cat("Using script at: ", STAN_file,"\n", sep = "") # message to user
    
    if(!stanc(STAN_file)$status){
      warning("Growth file for STAN script does not exist!\n")
      next
      }
  }
  
  # check that the file is
  
  # Actually run the MCMC
  t1 <- Sys.time() # start time
  cat("Start MCMC call at: ", as.character(t1), "\n\n", sep = "")
  
  fit1 <- stan(file = STAN_file,
               data = All_STAN_data.L[[f]],
               warmup = warm_its, 
               iter = smpl_its + warm_its,
               chains = n_chains,
               cores = n_coreuse,
               thin = n_thin,
               control = list(adapt_delta = NUTS_adapt_delta))
  
  t1.1 <- Sys.time()
  # run diagnostics
  #JMR_smry_diagnostics(fit1, fitname = f)
  
  All_STAN_fits.L[[f]] <- fit1
  
  t2 <- Sys.time() # end time
  
  # print message for duration
  cat("Ended MCMC call at:",
      as.character(t1.1),"\n",
      sep = "")
  cat("Time needed for ", f,"\n",
      "MCMC:  ",
      round(difftime(t1.1, t1, units = "mins"), 2),"mins\n",
      "Validation: ",
      round(difftime(t2, t1.1, units = "mins"), 2), "mins\n",
      "Total: ",
      round(difftime(t2, t1, units = "mins"), 2),
      "mins.\n\n")
}


### Storing results ###########################################################
save(list = c("All_STAN_fits.L", "All_params_SpeciesID.df"),
     file = paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/",site,
                   "_STAN_ResultsMCMC.R"))
}

# Timer for the entire script
t_fnl <- Sys.time()

cat("All runs completed.\n", 
    as.character(t_ini), " - Began\n",
    as.character(t_fnl), " - Ended\n",
    "Total duration: ",
    round(difftime(t2, t1, units = "mins"), 2),
    "min\n", sep = "")



Sys.sleep(3)
# End