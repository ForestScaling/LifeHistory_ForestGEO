# File: 04_STAN_ModelChecksV3.R
# Author: JMR 
# Date: 12023/07/16

# ### Setup ###################################################################
set.seed(123) # Set seed if necessary
site.v <- c("HVDF", "SCBI", "SERC") # vector of the site codes

# location of the data and models
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_stanout <- paste0(loc_Gdr, "/data/ForestGEO/STAN_outputs/") # where stan
loc_out <- "C:/Users/juan.m.rodriguez/Downloads/" # Where to dump diagnostics
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr", "scales",
           "moments",
           "rstudioapi", "tidyverse",
           "cowplot", "zoo", "bayesplot")){
  library(x, character.only = T)
  }

# Loading in diagnostics script written by M. Betancourt
cat("Loading MCMC NUTS diagnostics script written by M. Betancourt...\n")
source('./scripts/STAN/stan_utility_rstan.R')

# JMR's diagnostics (General)
JMR_smry_diagnostics <- function(STANfit, fitname = "Model"){
  cat("Running dianostics for ", fitname,":\n", sep ="")
  
  check_all_expectand_diagnostics(STANfit)
  check_all_hmc_diagnostics(STANfit)
  
  cat("Check Energy:\n")
  k <- check_energy(STANfit)
  if(!is.null(k)){
    print(k)
  }else{
      cat("E-BFMI good.\n")}

}
# --- Check for locations of output files -------------------------------------
loc_out <- paste0(loc_out,
                  "diagnostics_",
                  Sys.Date())

# verifying the output location's exsitence
if(file.exists(loc_out)){
  cat("Output for directory diagnostics: '",
      loc_out,
      "' already exists, this may override plots!\n",
      sep = "")
}else{
    dir.create(loc_out)
  }

# verify the site folders
if(!all(file.exists(paste0(loc_out,"/",site.v)))){
  # create the files
  walk(paste0(loc_out,"/",site.v)[
    which(!file.exists(paste0(loc_out,"/",site.v)))],
    ~dir.create(.x))
}



# --- Site specific setup -----------------------------------------------------
site <- "HVDF"
for(site in site.v){
load(paste0(loc_stanout, site, "_STANData_MCMC.r")) # load the data used
load(paste0(loc_stanout, site, "_STAN_ResultsMCMC.R")) # load the fits






fit_names <- All_STAN_data.L %>%
  names(.) # names of all the stand Data sets

for(f in fit_names){  # Name of the model
  cat("Currently generating plots for ", site,"_",f,":\n\n", sep ="")
  
f_type <- gsub(x = f, pattern = "_CL\\d", replacement = "")

fit1 <- All_STAN_fits.L[[f]] # the fit
stan_data <- All_STAN_data.L[[f]] # the data for the model

# print(identical(length(unique(stan_data$species)), stan_data$n_spp))
# 
# }}


curr_params.df <- All_params_SpeciesID.df %>%
  filter(param == f) # TempID -Mnemonic data

# number of iterations 
n_iter_a <- sum(fit1@sim$n_save - fit1@sim$warmup2)
n_iter_c <- median(fit1@sim$n_save - fit1@sim$warmup2)

# choose which parameters to investigate most closely
# we will use the 2 most observed the most oberservations
curr_params.df <- curr_params.df %>%
  arrange(desc(n_obs))

g <- 3 # must be greater than 2

# check_all_hmc_diagnostics(fit1)
# check_all_expectand_diagnostics(fit1)

par_names <- names(fit1)[!(grepl(pattern = "log_lik|g_rep",
                                 x= names(fit1)))]
 



# Get the 'hyper parameters'
base_params <- switch(f_type,
                     growthBinary = c("mean_odds", "sd_odds"),
                     growthCont = c("mu_growparam_CL",
                                    "sigma_growparam_CL",
                                    "logmu_growsd_CL",
                                    "logsigma_growsd_CL"),
                     BinarySurv = c("mean_odds", "sd_odds"))


# chooses the g spaced-out positions (e.g. g = 5, the quantiles + min + max)
low_param_number <- c(1, round(quantile(1:nrow(curr_params.df),
                         probs = ((2:(g-1))/g))),nrow(curr_params.df))

# pick out the specific letters
low_param <- par_names[grepl(pattern = paste(paste0("\\[",
                                                    low_param_number,"\\]"),
                                     collapse = "|"), x = par_names)]

# separate the other duplicate parameters (e.g. we don't need both odds
# ratio and p for binary data since one is just a transformation of the other)
low_param <- low_param[str_which(low_param, pattern = switch(f_type,
                                     growthBinary = "odds_rat|odds",
                                     growthCont = "growsd_sp_CL|growparam_sp_CL",
                                     BinarySurv = "odds_rat"))]

# add them to the parameters to use
use_params <- c(base_params, low_param, "lp__")



raw_MCMC_iter.df <- as.matrix(fit1) %>%
  as_tibble %>%
  mutate(iter = row_number(.),
         divergent = get_divergent_iterations(fit1))

# Extract the reads for each parameter, stores the chain as well
raw_MCMC_iter_wChain.df <- bayesplot::mcmc_trace_data(fit1)

# --- Generate the outputs (All models) -----------------------------------------
plt.L <- list()

# Diagnostics (text)

k <- getOption("width")
options("width" = 250)
sink(file = paste0(loc_out, "/", # output text file
            site, "_", f, "_QuickDgnst.txt"))

cat("STAN Model ", site," ",f, ":\n", sep = "")

cat("Model Info (chains, iterations, thinning:\n")
print(fit1@sim[2:4]) # print n-chains, n-iter, and n-thin

cat("Quick diagnostics from STAN and M. Betancourt:\n")
JMR_smry_diagnostics(fit1, fitname = paste(site, f, sep = "_"))

cat("\nSummary:\n")
print(summary(fit1)$summary)
sink(NULL)
options("width" = k)

# Generate mcmc_trace_highlight (rank "hist", but easier to see)
plt.L <- list(mcmc_ranks.plt <-
    mcmc_rank_overlay(fit1, pars = use_params) +
                         theme_bw() +
                         scale_x_continuous(
                           labels =label_percent(scale = 100/n_iter_a),
                         n.breaks = 5,
                         name = "Rank (percentile of iterations)") +
      scale_y_continuous(labels = label_percent(scale = 100/n_iter_c),
                         n.breaks = 5,
                         limits = c(0,n_iter_c/10)) +
      ggtitle(paste0(site,"_", f, "_rankHist")),
  
      width = 11.75,
      height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# trace plot, not super helpful by itself, just make sure you don't have giant
# chunks stuck in one spot (sign of poor mixing)
plt.L <- list(
  mcmc_trace_data(x = fit1, pars = use_params) %>%
    ggplot(., aes( x = iteration, y = value, color = chain)) +
    geom_line() +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE),
                       n.breaks = 4) + 
    scale_color_viridis_d(alpha = 0.5, option = "G") +
    facet_wrap("parameter", scales = "free") +
    theme_bw() +
    ggtitle(paste0(site, "_", f, "_trace")),
  width = 11.75,
  height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# Autocorrelation vs Lag plot, make sure it rapidly approaches 0.
plt.L <- list(mcmc_acf.plt <- mcmc_acf(x = fit1,
                                       pars = use_params,
                                       lags = 25)$data %>%
                mutate(Chain = as.factor(Chain)) %>%
                ggplot(data = ., aes(x = Lag, y = AC, color = Chain)) +
                geom_line() +
                facet_wrap(facets = "Parameter") +
                scale_y_continuous(name = "Autocorrelation") +
                scale_color_viridis_d(alpha = 0.5, option = "G") +
                theme_bw() +
                ggtitle(paste0(site, "_", f, "_lagAcf")),
              width = 11.75,
              height = 7.5) %>%
                list() %>%
                c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# plot showing the running average across multiple chains
plt.L <- list(mcmc_CmlAvg.plt <- raw_MCMC_iter_wChain.df %>%
                arrange(iteration) %>%
                mutate(.by = c(parameter, chain),
                       w = 1,
                       roll_avg = cumsum(value)/cumsum(w),
                       roll_max = cummax(value),
                       move_avg = rollmean(x = value,
                                           k = 10,
                                           fill = NA,
                                           align = "right")) %>%
                select(-c(w)) %>%
                filter(parameter %in% use_params) %>%
                ggplot(data = .,
                       aes(x = iteration,
                           y = roll_avg,
                           color = chain)) +
                geom_line() + 
                scale_color_viridis_d(alpha = 0.5, option = "G") +
                facet_wrap("parameter", scales = "free") +
                scale_x_continuous(n.breaks = 4) +
                scale_y_continuous(name = "Cummulative average") +
                theme_bw() + 
                ggtitle(label = paste(site, f, "CmlAvg", sep = "_")),
              width = 11.75,
              height = 7.5) %>%
                list() %>%
                c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# ratio of Neff:N, should be >1:10
# mcmc_NRatio.plt <- mcmc_neff(neff_ratio(fit1, pars = par_names)) +
#   yaxis_text(hjust = 0) +
#   ggtitle(label = paste(site, f, "NRatio", sep = "_"))
# 
# # Rhat diagnostic
# mcmc_RHat.plt <- mcmc_rhat(rhat(fit1, pars = par_names)) +
#   yaxis_text(hjust = 0) +
#   scale_x_continuous(limits = c(0.95, 1.2)) +
#   ggtitle(label = paste(site, f, "Rhat", sep = "_"))

 
# Effective sample size
plt.L <- list(mcmc_NRatioRhat.plt <-
 data.frame(param = par_names,
           ess_bulk = apply(as.array(fit1, pars = par_names),
                            MARGIN = 3,
                            FUN = ess_bulk),
           ess_tail = apply(as.array(fit1, pars = par_names),
                            MARGIN = 3,
                            FUN = ess_tail)) %>%
   mutate(label = if_else(ess_bulk < 1000 | ess_tail < 1000,
                         param, "")) %>%
   ggplot(., aes(x = ess_bulk, y = ess_tail)) +
   geom_point() +
   # label any past thresholds
   geom_text(aes(label = label), hjust = 0, vjust = -.75) +
   # axes
   scale_x_continuous(limits = c(NA, NA)) +
   scale_y_continuous(limits = c(NA, NA)) +
   # lines for thresholds, dashed is a "concern", solid is bad
   geom_vline(xintercept = c(100,1000),
             linetype = c(1,2),
             color = "navy") +
   geom_hline(yintercept =  c(100,1000),
             linetype = c(1,2),
             color = "navy") +
   theme_bw() +
   ggtitle(label = paste(site, f, "ESS", sep = "_")),
 width = 7.5,
 height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

#rhat and Nrat
plt.L <- list(ESS.plt <-
                data.frame(param = par_names,
                           Nratio = neff_ratio(fit1, pars = par_names),
                           Rhat = rhat(fit1, pars = par_names)) %>%
                mutate(label = if_else(Rhat > 1.05 | Nratio <= 0.5,
                                       param, "")) %>%
                ggplot(., aes(x = Nratio, y = Rhat)) +
                geom_point() +
                # label any past thresholds
                geom_text(aes(label = label), hjust = 0, vjust = -.75) +
                # axes
                scale_x_continuous(limits = c(0, NA), name = "Neff/N") +
                scale_y_continuous(limits = c(NA, NA)) +
                # lines for thresholds, dashed is a "concern", solid is bad
                geom_vline(xintercept = c(.1,0.5),
                           linetype = c(1,2),
                           color = "navy") +
                geom_hline(yintercept =  c(1.05,1.1),
                           linetype = c(2,1),
                           color = "navy") +
                theme_bw() +
                ggtitle(label = paste(site, f, "NRat_v_Rhat", sep = "_")),
              width = 7.5,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)


# Nuts Energy
plt.L <- list(mcmc_EDistr.plt <- mcmc_nuts_energy(nuts_params(fit1),
                 bayesplot::log_posterior(fit1), binwidth = 1) +
                theme(plot.background = element_rect(fill = "white"))+
                ggtitle(label = paste(site, f, "EDistr", sep = "_")),
              width = 11.75,
              height = 7.5) %>%
                list() %>%
                c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# Pairs plot
plt.L <- list(mcmc_pairs.plt <-
                mcmc_pairs(fit1, pars = c(base_params, "lp__"),
                           off_diag_args = list(size = .5,
                                                alpha = .1),
                           lp = log_posterior(fit1),
                           np = nuts_params(fit1),
                           np_style = pairs_style_np()),
              width = 16,
              height = 12) %>%
  list %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# Parcoord


parc_dat <- mcmc_parcoord_data(
  x = as.array(fit1, pars = c(low_param, base_params)),
  pars = c(low_param, base_params),
  np = nuts_params(fit1))

div_draws <- summarize(parc_dat,
                       .by = "Draw",
                       Divergent = sum(Divergent)) %>%
  filter(Divergent != 0) %>%
  pull(Draw)


keep_draws <- sample(with(filter(parc_dat,
                                 !(Draw %in% div_draws)),
                          unique(Draw)),
                     size =  1000,
                     replace = F) %>%
  c(., div_draws)

# Add the levels
parc_dat %<>% mutate(Divergent = as.factor(Divergent))
levels(parc_dat$Divergent) <- c("0","1")


plt.L <- list(mcmc_parcoord.plt <- parc_dat %>%
  filter(Draw %in% keep_draws) %>%
  ggplot(., aes(x = Parameter,
                y = Value,
                group = Draw,
                color = Divergent)) +
  geom_line() +
  scale_color_manual(values = c("navy", "red")) +
  ggtitle(label = paste(site, f, "ParCoord", sep = "_")) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)),
height = 7.5, width = 7.5) %>%
  list %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)
# --- Posterior predictive check ----------------------------------------------
cat("Generating PostPC plots...\n\n")
n_ysim <- 100 # number of simulations for posterior pc density plots

# For the binary variables
if(f_type %in% c("growthBinary", "BinarySurv")){
post_p.m <- as.matrix(fit1, pars = "p")

time.v <- stan_data$t_interv
spp.v <- stan_data$species
n_spp <- stan_data$n_spp
n_stem <- stan_data$n_stem


# vector for the realized p
p_rlz.v <- vector(mode = "numeric", length = stan_data$n_stem)

# current random draw
sim_curr_res.df <- data.frame(spp = spp.v,
                              res = NA)

# data.frame of the results of p for each spp.
post_pc1.m <- matrix(data = NA,
                     nrow = n_ysim,
                     ncol = n_spp)

post_pc2.m <- matrix(data = NA,
                     ncol = n_stem,
                     nrow = n_ysim)

# vector of random iterations to use for the simulated poster draws
s.v <- sample(x = 1:nrow(post_p.m),
              size = n_ysim,
              replace = F)

# generate the y_sims
for(i in 1:n_ysim){
  
  s <- s.v[i] # which iteration of the MCMC to use
  
  # realized prob
  p_rlz.v <- post_p.m[s, c(spp.v)]^(time.v/5)
  
  # rng the draw
  post_pc2.m[i,] <- rbinom(n = n_stem,
                           size = 1,
                           prob = p_rlz.v)
  
  # summarize to get the p at each simulation
  post_pc1.m[i,] <- post_pc2.m[i,] %>%
    as_tibble %>%
    mutate(spp = spp.v) %>%
    rename(res = value) %>%
    summarize(., n = n(),
              p = sum(res)/n(),
              .by = "spp") %>%
    pull("p")
  
}
  dat <- data.frame(spp = stan_data$species,
                    bindat = stan_data$bindat) %>%
    summarise(p = sum(bindat)/n(), .by = spp)
  
  # function for retrieving the proportion
  ppn <- function(x){return(sum(x)/length(x))}
  
  plt.L <- list(ppc_all.plt <- ppc_dens_overlay(y = dat$p, yrep = post_pc1.m) +
    ggtitle(label = paste0(site, "_", f, "_PostPC_distrPpn")) +
    ylab("Density") +
    theme(plot.background = element_rect(fill = "white"))+
    xlab("p"), width =7.5, height = 7.5) %>%
    list() %>% c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  plt.L <- list(ppc_grouped.plt <- ppc_stat_grouped(yrep = post_pc2.m,
                   y = stan_data$bindat,
                   group = spp.v, stat = ppn, freq = F) +
    theme(plot.background = element_rect(fill = "white"))+
    ggtitle(label = paste0(site, "_", f, "_PostPC_ppnBySpp")) +
    scale_x_continuous(n.breaks = 3),height = 12, width = 16) %>%
    list() %>% c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
 
  
}else{
  # Statement for PostPC for lognormal
  dat <- data.frame(species = stan_data$species,
                    vals = stan_data$loggrowth_vals)
  
  # Plot other diagnostics:
  posterior1 <- rstan::extract(fit1) # get the posterior of the fit
  posterior_array1 <- as.array(fit1)
  g_reps <- fit1 %>%
    as.matrix(., pars = "g_rep") # reps
  rep_rows <- sample(nrow(g_reps),
                     size = 100,
                     replace = F)
  

  # All vals PostPC
plt.L <- list(PostPCAllVals.plt <- ppc_dens_overlay(dat$vals,
                   g_reps[rep_rows,]) +
    scale_y_continuous(limits = c(0, 1)) +
    ylab("Density") +
    xlab("Value") +
      theme(plot.background = element_rect(fill = "white"))+
    ggtitle(label = paste(site, f, "PostPCAllVals", sep = "_")),
    width = 7.5,
    height = 7.5) %>%
  list %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# skew PostPC  
  plt.L <- list(PostPCSkew.plt <- ppc_stat_grouped(dat$vals,
                   g_reps, 
                   stat = skewness,
                   group = dat$species) +
    scale_x_continuous(n.breaks = 3) +
      theme(plot.background = element_rect(fill = "white"))+
    ggtitle(paste(site,f,"PostPCSkew", sep ="_")),
  width = 16,
  height = 12) %>%
  list %>%
  c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  # Kurtosis PostPC
  plt.L <- list(PostPCKurt.plt <- ppc_stat_grouped(dat$vals,
                                     g_reps, 
                                     stat = kurtosis,
                                     group = dat$species) +
    scale_x_continuous(n.breaks = 3) +
      theme(plot.background = element_rect(fill = "white"))+
    ggtitle(paste(site,f,"PostPCKurt", sep ="_")),
    width = 16,
    height = 12) %>%
    list %>%
    c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
}





# Create directory for the model in R
if(!dir.exists(paste(loc_out, site, f, sep ="/"))){
  # create the folder for the specific model
  dir.create(paste(loc_out, site, f,sep = "/"))
}

# Loop through the plots and output them
cat("Outputing plots...")


for(plt_i.L in plt.L){
  # extract the plot
  plt <- plt_i.L[[1]]
  
  if(is.null(plt$labels$title)){
    ttl <- paste(site,f,"Pairs", sep = "_")
    }else{
      ttl <-  plt$labels$title
    }
  ggsave(plot = plt,
         filename = paste0(paste(loc_out,
                      site,f,
                      ttl, sep = "/"), ".png"),
         device ="png", width = plt_i.L$width,
         height = plt_i.L$height)
}

cat(" done.\n")


}


}