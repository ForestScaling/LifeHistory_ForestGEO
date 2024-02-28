# location of the data
loc_data <- paste0("G:/Shared drives/NASA_ForestScaling/",
                   "LifeHistoryStrats/",
                   "data/ForestGEO/processed/")

# load in the data
load(paste0(loc_data,"exampleData.R"))

# Packages
invisible(lapply(list("tidyverse",
                      "rstan",
                      "parallel"),
                 FUN = require,
                 character.only = T))

# some Rstan options I like to use
rstan_options(auto_write = TRUE,  # whether to write the compiled C++ script
                                  #   prevents recompiling each time (?)
              javascript = FALSE) # prevents using the github C++ compiler
                                  #   on TRUE, will use the internet to access
                                  #   to the github server which has caused me
                                  #   issues before

# Default vals for runs (override individually if desired.)
warm_its <- 2000  # number of warmups
smpl_its <- 5000  # number of USED iter, (MCMC will run for smpl + warm)
n_coreuse <- 5    # number of cores to use
n_chains <- 5     # number of chains to run
n_thin <- 1       # thinning, 1 = none
NUTS_adapt_delta <- 0.85 # adapt delta, default is 0.801, larger means
                         # ~smaller "steps" so less divergence, but less
                         # efficient, 

# Prep Data for stan
# stan likes data in lists, where the name of each item in the list 
# corresponds to however it is named in the 'data' block of the script.

# Growth data
exmpl_grow.L <- list(
  # number of stems, an integer
  n_stem = nrow(exmpl_grow.df),
  # number of species, an integer
  n_spp = length(unique(exmpl_grow.df$species)),
  # vector of the species integers
  species = exmpl_grow.df$species,
  # vector of growth values (log transformed)
  loggrowth_vals = exmpl_grow.df$grow_val) 

# Survival/binary data
exmpl_surv.L <- list(
  # number of stems, an integer
  n_stem = nrow(exmpl_surv.df),
  # number of species, an integer
  n_spp = length(unique(exmpl_surv.df$species)),
  # vector of the species integers
  species = exmpl_surv.df$species,
  # vector of binary survival values
  bindat = exmpl_surv.df$survival,
  # Length of time between observations
  t_interv = exmpl_surv.df$t)

# --- Run the models ----------------------------------------------------------
### Growth Model ###
t_i <- Sys.time()

# Run the fit in STAN
fit_grow <- stan(file = "./STAN/growth_lognormal_estimate.stan", # stan script
                 data = exmpl_grow.L,        # list of data
                 warmup = warm_its,          # burn-in
                 iter = smpl_its + warm_its, # number of iterations (total)
                 chains = n_chains,          # number of chains
                 cores = n_coreuse,          # number of cores to use
                 thin = n_thin,              # thinning
                 control = list(adapt_delta = NUTS_adapt_delta))
fit_grow@model_name <- "exmpl_Growth" # update the name to the model 

# print duration
cat("Growth run time: ",
    round(difftime(Sys.time(), t_i, units = "mins"), 2),
    "min\n", sep = "")

### Survival/Binary Data Model ###
t_i <- Sys.time()

# Run the fit in STAN
fit_surv <- stan(file = "./STAN/binary_logitparam.stan", # stan script
                 data = exmpl_surv.L, # data used for the fit
                 warmup = warm_its,          # burn-in
                 iter = smpl_its + warm_its, # number of iterations (total)
                 chains = n_chains, # number of chains
                 cores = n_coreuse, # number of cores run
                 thin = n_thin, # thinning chains
                 control = list(adapt_delta = NUTS_adapt_delta))
fit_surv@model_name <- "exmpl_Surv" # fix model name

# print duration
cat("Survival run time: ",
    round(difftime(Sys.time(), t_i, units = "mins"), 2),
    "min\n", sep = "")


# Save the outputs
save(list = c("fit_surv", "fit_grow"),
     file = paste0("G:/Shared drives/NASA_ForestScaling/",
                   "LifeHistoryStrats/",
                   "data/ForestGEO/STAN_outputs/example_STANOutput.R"))

