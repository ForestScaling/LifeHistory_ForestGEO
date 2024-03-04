# File: 04_STAN_ModelChecks.R
# Author: JMR 
# Date: 12023/07/16
# ### Setup ###################################################################
set.seed(123) # Set seed if necessary

# location of the data and models
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr",
           "rstudioapi", "tidyverse",
           "cowplot")){
  library(x, character.only = T)}

# Loading in diagnostics written by M. Betancourt
source('scripts/STAN/stan_utility_rstan.R')

# --- Site specific setup -----------------------------------------------------
# Load in the site's data
site <- "HVDF"
load(paste0(loc_data, site, "_STANdata.r")) # load the data used (for PPC)
AllData.L <- get(paste0(site, "_STANdata.L")) # assign to a common name

load(paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/", site,
             "_STAN_MCMC.r"))



# --- Run checks for Growth ---------------------------------------------------
n_Lvls <- All_STAN_fits.L %>%
  names %>%
  str_extract(., pattern = "(?<=_CL)\\d+") %>%
  as.integer() %>% 
  max()

growth_all.df <- AllData.L$growth_data[[1]]

# Generate the plots
# parameters to viusally inspect
params <- c("growparam_sp_CL[3]", # names of the parameters to inspect
            "growparam_sp_CL[5]",
            "growsd_sp_CL[7]",
            "growsd_sp_CL[2]",
            "logmu_growsd_CL",
            "mu_growparam_CL")


for(Lvl in 1:n_Lvls){
  cat(paste0("Beginning check for Growth parameter in Canopy Lvl ",
             Lvl, "\n" ))
  
  fit1 <- growth_fitsAllLvl.L[[Lvl]]
  
  # get the original data that was used
  growth_CL.df <- growth_all.df %>% filter(CanopyLvl == Lvl)
  
  # Plot other diagnostics:
  posterior1 <- rstan::extract(fit1) # get the posterior of the fit
  posterior_array1 <- as.array(fit1)
  g_reps <- fit1 %>% as.matrix(., pars = "g_rep") # reps
  rep_rows <- sample(nrow(g_reps),
                     size = 100,
                     replace = F)
  
  # "Rank Plot" histogram of number of iterations per rank order amongst all
  # chains. Used to check that chains are mixing properly. If they are mixing
  # properly, we should expect to see that the histogram across all chains
  # should be approximately uniform.
  mcmc_rank_hist <- mcmc_rank_hist(fit1, pars = params) + 
    ggtitle(paste0(site,"_growCL",Lvl,"_rankHist"))
  
  # Let's look at the posterior predictive checks, and make sure our
  # model describes the data well.
  ppc_dens_obs <- ppc_dens_overlay(log(growth_CL.df$growth2),
                                   g_reps[rep_rows,]) +
    ggtitle(paste0(site,"_growCL", Lvl,"_PPC_allobs"))
  
  
  OTUS_check_table <- growth_CL.df %>%
    group_by(Mnemonic) %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(rare_groups = case_when(n >= 10 ~ Mnemonic,
                                   TRUE ~ "<10 'rare'")) %>% 
    left_join(growth_CL.df, ., by = "Mnemonic")
  
  # Posterior predictive check for observations
  ppc_obs <- ppc_dens_overlay_grouped(log(growth_CL.df$growth2),
                                      g_reps[rep_rows,],
                                      group = as.factor(OTUS_check_table$rare_groups)) +
    ggtitle(paste0(site,"_growCL",Lvl, "_PPC_Obs_bySpp"))
  
  # temporary dataframe to hold the values for each OTU
  df_mean <- as.data.frame(matrix(NA,
                                  ncol = length(unique(growth_CL.df$Mnemonic)),
                                  nrow = nrow(g_reps)))
  colnames(df_mean) <- unique(growth_CL.df$Mnemonic)
  df_sd <- df_median <- df_mean # duplicate for the other statistics
  
  
  for(OTU in unique(growth_CL.df$Mnemonic)){
    # Select all columns (obs), that correspond to the OTU
    temp.df <- g_reps[, which(growth_CL.df$Mnemonic == OTU)] %>%
      as.data.frame()
    
    # Take the statsitics 
    df_mean[, OTU] <- temp.df %>%
      apply(., MARGIN = 1, FUN = mean)
    df_median[, OTU] <- temp.df %>%
      apply(., MARGIN = 1, FUN = median)
    df_sd[, OTU] <- temp.df %>%
      apply(., MARGIN = 1, FUN = sd)
  }
  
  # dataframe holding onto the percentiles 
  df_PPC <- growth_CL.df %>%
    group_by(Mnemonic) %>%
    summarise(n = n(),
              CanopyLvl = Lvl,
              mean = mean(log(growth2)),
              median = median(log(growth2)),
              sd = sd(log(growth2))) %>%
    mutate(p_mean = NA, p_median = NA, p_sd = NA)
  
  for(i in 1:nrow(df_PPC)){
    OTU <- df_PPC$Mnemonic[i] # The OTU being used
    # mean
    myfun <- ecdf(df_mean[, OTU])
    df_PPC$p_mean[i] <- myfun(df_PPC$mean[i])
    # median
    myfun <- ecdf(df_median[, OTU])
    df_PPC$p_median[i] <- myfun(df_PPC$median[i])
    # sd
    if(is.na(df_PPC$sd[i])){
      next
    }else{
      myfun <- ecdf(df_sd[, OTU])
      df_PPC$p_sd[i] <- myfun(df_PPC$sd[i])}
  }
  
  
  # mean by sp.
  ppc1 <- ppc_stat_grouped(log(growth_CL.df$growth2),
                           g_reps, 
                           stat = "mean",
                           group = growth_CL.df$Mnemonic) +
    ggtitle(paste0(site, "_growCL", Lvl, "_PPC_avg_bySpp"))
  
  # median by sp.
  ppc2 <- ppc_stat_grouped(log(growth_CL.df$growth2),
                           g_reps, stat = "median",
                           group = growth_CL.df$Mnemonic) +
    labs(title = paste0(site,"_growCL",Lvl,"_PPC_med_bySpp"))
  
  # Standard deviation by sp.
  ppc3 <- ppc_stat_grouped(log(growth_CL.df$growth2),
                           g_reps,
                           stat = "sd",
                           group = growth_CL.df$Mnemonic) +
    ggtitle(paste0(site,"_growCL", Lvl, "_PPC_sd_bySpp"))
  
  # Auto correlation function
  plot_autocor <- mcmc_acf(posterior_array1, pars = params,
                           lags = 20) + ggtitle(paste0(site,"_growCL", Lvl, "_mcmc_autocor_lag"))
  
  # save location
  loc_sav <- "C:/Users/juan.m.rodriguez/Downloads/"
  
  # create a folder just for this set of plots
  loc_sav2 <- paste0(loc_sav,"diagn_growCL",Lvl)
  dir.create(loc_sav2)
  
  kableExtra::save_kable(x = knitr::kable(df_PPC, digits = 2),
                         file = paste0(loc_sav2, "/",site,"_growCL", Lvl,
                                       "_PPC_bySpp_table.html"))
  
  for(plot in list(ppc3, ppc1, ppc2,
                   plot_autocor, ppc_obs,
                   ppc_dens_obs, mcmc_rank_hist)){
    
    ggsave(filename = paste0(plot$labels$title, ".pdf"),
           plot = plot,
           device = "pdf",
           height = 11.32,
           width = 12.92,
           path = loc_sav2)
    
    
  }
  
  
}



