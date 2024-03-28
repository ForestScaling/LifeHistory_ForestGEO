# File: 05_AggrParams
# Author: JMR 
# Date: 12023/07/16
# ### Setup ###################################################################
set.seed(123) # Set seed if necessary

# location of the data and models
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/taxa/")
loc_stanout <- paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr",
           "rstudioapi", "tidyverse",
           "cowplot",
           "stringr")){
  library(x, character.only = T)}

# --- Site specific setup -----------------------------------------------------
# Load in the site's data
site <- "HVDF"

for(site in c("HVDF", "SCBI", "SERC")){
load(paste0(loc_data, site, "_STANdata.r")) # load the data used (for PPC)
#AllData.L <- STANdata.L # assign to a common name

#load(paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/", site,
#            "_STAN_MCMC.r"))

load(paste0(loc_stanout, site, "_STANData_MCMC.r")) # load the data used
load(paste0(loc_stanout, site, "_STAN_ResultsMCMC.R")) # load the fits

# load the taxa table
load(paste0(loc_data2, site, "_TaxaTable.r"))
TaxaTable <- Spp_table
use_taxa.v <- STANdata.L$stature$Mnemonic[!is.na(STANdata.L$stature$number)]

# --- Extract summaries of the posterior --------------------------------------

# table for holding the values
param_ests.df <- TaxaTable %>%
  # only have taxa we considiered for parametrization
  filter(!as.logical(DropParamtrz) & !as.logical(DropBegin)) %>%
  select(., Mnemonic, Latin, IDlevel) %>%
  # Make sure that there are >=1 stem for the species to be included,
  # the texanomic table will include stuff that would be excluded due to
  # size cutoffs
  filter(Mnemonic %in% use_taxa.v)


# Logical check that the names of the parameters match up with the fits
if(!identical(sort(names(All_STAN_fits.L)),
              sort(names(All_STAN_data.L))) |
   !identical(sort(unique(All_params_SpeciesID.df$param)),
              sort(names(All_STAN_data.L)))){
  stop(paste0("Parameter names in the MCMC data, parameter-species",
              " assignment, and/or MCMC fits don't all match up!"))
}


# Common strings of text for parameters to disregard (log likelihood and
# simulation from the posterior) will use text matching to remove these
param_2_disregard <- c("g_rep", "log_lik", "lp__")

sumrz_MCMC_ALL <- data.frame() # Make an empty data frame to hold everything in it

# loop through the MCMCs to extract the summarized posteriors for each param
for(k in 1:length(All_STAN_fits.L)){

  # Name the parameter being used
  data_type_k <- names(All_STAN_fits.L)[k]
  cat(paste0("Summarizing - ", data_type_k,"\n"))

  # get the fit of the particular parameter
  param_fit <- All_STAN_fits.L[[k]]

  # Names of the parameters of interest
  raw_param_names <- names(param_fit)[
    # only keep those parameters of interests (e.g. don't keep log_lik)
    !grepl(pattern = paste(param_2_disregard, collapse = "|"),
           x = names(param_fit))]

  # Generate a table summarizing for each of the parameters of interes
  sumrz_MCMC_spp <- param_fit %>%
    as.matrix(pars = raw_param_names) %>%
    as_tibble %>%
    pivot_longer(everything(),
                 names_to = "param",
                 values_to = "val") %>%
    group_by(param) %>%
    summarize(., ecdf_p0 = ecdf(val)(0))
  
  sumrz_MCMC_spp <- summary(param_fit,
                        probs = c(0.025, 0.5, 0.975), # Using a 95% Credible Interval
                        pars = raw_param_names)$summary %>%
    as.data.frame() %>%
    mutate(., param = rownames(.)) %>%
    rename_with(~c("post_qLwr", "post_qMed", "post_qUpr"),
                .col = contains("%")) %>%
    mutate(wCrI =  post_qUpr - post_qLwr) %>%
    left_join(., sumrz_MCMC_spp, by = "param") %>%
    relocate(., param, .before = everything()) %>% # include with ecdf_0 values
    as_tibble()

  # separate these out into the hyperparameters and the OTU specific
  sumrz_MCMC_spp <- sumrz_MCMC_spp %>%
    mutate(is_sppParam = grepl("[[:digit:]]", param),
           tempID = case_when(!is_sppParam ~ NA_character_,
             is_sppParam ~ str_extract(string = param,
                                       pattern = "(?<=\\[).*(?=\\])")),
           param_label = str_replace(param,
                                     pattern = "\\[.*\\]",
                                     replacement = ""))




  # separate the hyperparameters (may be used for the "mean")
  sumrz_MCMC_hyp <- sumrz_MCMC_spp %>%
    filter(!is_sppParam)
  sumrz_MCMC_spp <- sumrz_MCMC_spp %>%
    filter(is_sppParam) 

  
  # link up with the temporaryID for each taxon
  sumrz_MCMC_spp <- All_params_SpeciesID.df %>%
    filter(param == data_type_k) %>%
    select(-c("param")) %>%
    mutate(tempID = as.character(tempID)) %>%
    left_join(., sumrz_MCMC_spp, by = "tempID") %>%
    mutate(param_label = paste0(param_label, CanopyLvl),
           param_est = TRUE)
  
  print(unique(sumrz_MCMC_spp$CanopyLvl)) # print canopy Level just as a check
  
  # add the canopyLVl data into th
  sumrz_MCMC_hyp %<>% mutate(CanopyLvl = unique(sumrz_MCMC_spp$CanopyLvl))
  
  # get the 
  sumrz_MCMC_hyp <- sumrz_MCMC_hyp %>% as.data.frame()
  rownames(sumrz_MCMC_hyp) <- sumrz_MCMC_hyp$param
  
  # Fill in the "missing" using estimates of the mean from the hyper parameter
  #I NEED TO ADD MISSING VALUES DO A SET OF IFS
  if(grepl("growthCont_CL", x = data_type_k)){
    # Get expected value of the group parameters for growth sd and average
    summrz_filler <- data.frame(param_label = paste0(c("growparam_sp_CL",
                                                       "growsd_sp_CL"),
                                       unique(sumrz_MCMC_spp$CanopyLvl)),
                        mean = c(sumrz_MCMC_hyp["mu_growparam_CL","mean"],
                                 exp(sumrz_MCMC_hyp["logmu_growsd_CL","mean"] +
                                       ((sumrz_MCMC_hyp["logsigma_growsd_CL",
                                                        "mean"])^2)/2)))

  }else if(grepl("growthBinary_CL", x = data_type_k)){
    
    # fix name issue for parameter
    sumrz_MCMC_spp %<>% mutate(param_label = gsub(pattern = "odds_rat",
                                                  replacement = "godds_rat",
                                                  x = param_label))
    
    # We don't need actual survival prob, just log odds
    sumrz_MCMC_spp %<>% filter(grepl(param_label, pattern = "odds_rat"))
    
    # extract mean and var
    a <- sumrz_MCMC_hyp["mean_odds","mean"]
    b <- sumrz_MCMC_hyp["sd_odds","mean"]
    
    # fill in expected value for p
    summrz_filler <- data.frame(
      param_label = paste0("godds_rat", unique(sumrz_MCMC_spp$CanopyLvl)),
      mean = a)
    
  }else if(grepl("BinarySurv_CL", x = data_type_k)){
    sumrz_MCMC_spp %<>% mutate(param_label = gsub(pattern = "odds_rat",
                                                  replacement = "sodds_rat",
                                                  x = param_label))
    
    # remove the non-odds ratio
    sumrz_MCMC_spp %<>% filter(grepl(param_label, pattern = "odds_rat"))
    
    
    # extract alpha and beta
    a <- sumrz_MCMC_hyp["mean_odds","mean"]
    b <- sumrz_MCMC_hyp["sd_odds ","mean"]
    
    # fill in expected value for p
    summrz_filler <- data.frame(
      param_label = paste0("sodds_rat", unique(sumrz_MCMC_spp$CanopyLvl)),
      mean = a)
    
    
  }else{
    stop("Parameter/Data type can't handle missing values!")
    }
  
  # What OTU are Missing
  miss_OTU <- param_ests.df$Mnemonic[
    !(param_ests.df$Mnemonic %in% sumrz_MCMC_spp$Mnemonic)]
  
  summrz_filler <- data.frame(Mnemonic = sort(rep(miss_OTU,
                                                  nrow(summrz_filler))),
                              summrz_filler) %>% mutate(param_est = F,
                                                        n_obs = 0,
                                                        CanopyLvl = unique(sumrz_MCMC_spp$CanopyLvl))
  
  
  
  
  
  # add "estimated"/filler values to the main ones
  sumrz_MCMC_spp <- bind_rows(sumrz_MCMC_spp, summrz_filler) %>%
    select(-c("tempID", "ecdf_p0", "is_sppParam"))
  
  # Add to the collection
  sumrz_MCMC_ALL <- bind_rows(sumrz_MCMC_ALL, sumrz_MCMC_spp)
  
}

# Add weights
W_min <- 10^(-6)
W_max <- 1

sumrz_MCMC_ALL %<>%
  group_by(param_label) %>%
  summarise(wCrI99 = quantile(wCrI, 
                              probs = 0.99,
                              na.rm = T)) %>% 
  left_join(sumrz_MCMC_ALL, ., by = "param_label") %>%
  mutate(W = 1 - (wCrI/wCrI99)) %>% # calculate weights
  mutate(W = case_when(W < 0 & !is.na(W) ~ W_max,
    is.na(W) | W < W_min ~ W_min,
                       TRUE ~ W))

# --- Add parameters for non-bayesian estimated parameters --------------------
data_All <- sumrz_MCMC_ALL %>%
  select(Mnemonic, param_label, mean, W, param_est) %>%
  rename(val = mean, param = param_label)

# Adding stature
data_All <- STANdata.L$stature %>%
  filter(Mnemonic %in% use_taxa.v) %>%
  filter(!is.na(number)) %>%
  as_tibble() %>%
  transmute(Mnemonic = Mnemonic,
            param = "log_stature",
            val = case_when(is.na(statureDBH) ~ mean(log(statureDBH), na.rm = T),
                            !is.na(statureDBH) ~ log(statureDBH)),
            W = case_when(is.na(statureDBH) ~ W_min,
                          !is.na(statureDBH) ~ 1),
            param_est = !is.na(statureDBH)) %>%
  bind_rows(data_All, .)

# adding two types of recruitment
cat("Getting recruitment data from first pair of censuses only.\n")
keep_rcrt <- c("recrt_p1_stemsperMBA", "recrt_stemsperMInd")
recrt_df <- STANdata.L$recruitment[[1]] %>%
  as_tibble %>%
  filter(Mnemonic %in% use_taxa.v)

colnames(recrt_df)[colnames(recrt_df) == "area_stem_mat"] <- "MBA_Weight"
colnames(recrt_df)[colnames(recrt_df) == "n_stem_mat"] <- "MSt_Weight"
#***
# Recruitment weights, check later
recrt_df %<>%
  mutate(MBA_Weight = case_when(
    is.na(MBA_Weight) ~ NA_real_,
    !is.na(MBA_Weight) ~ 1 - log(MBA_Weight+1)/quantile(log(MBA_Weight+1),
                                                        na.rm = T,
                                                        probs = 0.99)),
    MBA_Weight = case_when(MBA_Weight < 0 ~ 1,
                           MBA_Weight < W_min ~ W_min,
                           TRUE ~ MBA_Weight),
    MSt_Weight = case_when(
      is.na(MSt_Weight) ~ NA_real_,
      !is.na(MSt_Weight) ~ 1 - log(MSt_Weight)/quantile(log(MSt_Weight),
                                                          na.rm = T,
                                                          probs = 0.99)),
    MSt_Weight = case_when(MSt_Weight < 0 ~ 1,
                           MSt_Weight < W_min ~ W_min,
                           TRUE ~ MSt_Weight))
# Weights are whack right now
recrt_df %<>% pivot_longer(c(recrt_p1_stemsperMInd, recrt_p1_stemsperMBA),
                          names_to = "param", values_to = "val") %>%
  mutate(W = case_when(param == "recrt_p1_stemsperMBA" ~ MBA_Weight,
                       param == "recrt_p1_stemsperMInd" ~ MSt_Weight),
         param_est = !is.na(W),
         W = case_when(is.na(W) ~ W_min,
                       TRUE ~ W))
### HERE ########

# keep only the used columns
recrt_df %<>% select(Mnemonic, param, val, W, param_est)

recrt_df %<>%
  group_by(param) %>%
  summarize(m_val = mean(val, na.rm = T)) %>%
  left_join(recrt_df, ., by = "param") %>%
  mutate(val = case_when(is.na(val) ~ m_val,
                         !is.na(val) ~ val))

# Scraping the weights for now, just setting those that weren't
# estimated to the acerage and tanking their weight
recrt_df %<>% mutate(W = case_when(param_est ~ 1,
                                  !param_est ~ W_min)) 

# Getting riof the m_val for filling in un-estimated ones
recrt_df %<>% select(-c(m_val))

data_All %<>% bind_rows(., recrt_df)

# Ovwerview of how much data is there:
data_All %>%
  group_by(param) %>%
  summarize(prop_w_estimates = sum(as.integer(param_est))/n()) %>%
  print
data_All %>%
  group_by(Mnemonic) %>% summarize(prop_w_estimates = sum(as.integer(param_est))/n()) %>% View()


# Output the table
save(list = "data_All",
     file = paste0(loc_Gdr, "/data/ForestGEO/processed/param_matrix/",
                   site, "_vals_wPCA.r"))
}