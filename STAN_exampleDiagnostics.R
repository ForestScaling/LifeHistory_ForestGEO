# --- Options -----------------------------------------------------------------
option_WriteOutputDiagn <- T
output_loc <- "C:/Users/juan.m.rodriguez/Downloads/test1"

# --- Setup -------------------------------------------------------------------
# Load the packages
invisible(lapply(c("bayesplot", "rstan", "gdata", "magrittr",
                   "scales", "moments", "rstudioapi", "tidyverse",
                   "cowplot", "zoo", "bayesplot"),
                 FUN = require,
                 character.only = T))

# Read the outputs from the model fits
load(file = paste0("G:/Shared drives/NASA_ForestScaling/",
                   "LifeHistoryStrats/",
                   "data/ForestGEO/STAN_outputs/example_STANOutput.R"))

# read in the original data
load(paste0("G:/Shared drives/NASA_ForestScaling/",
            "LifeHistoryStrats/",
            "data/ForestGEO/processed/exampleData.R"))


# read Betancourt's script
# Pulled from Michael Betancourt's github, betanalpha on 2023-06-02
# betanalpha/mcmc_diagnostics/rstan/stan_utility_rstan.R
source("./STAN/stan_utility_rstan.R")

# --- Functions ---------------------------------------------------------------

# My first pass, runs Betancourt's diagnostics, and the EBMFI diagnostic
JMR_smry_diagnostics <- function(STANfit, fitname = "Model"){
  cat("Running dianostics for ", fitname, ":\n",
      "Start: ", as.character(Sys.time()), "\n", sep = "")
  # Diagnostics (not mine)
  check_all_expectand_diagnostics(STANfit)
  check_all_hmc_diagnostics(STANfit)
  check_energy(STANfit)
  cat(capture.output(check_energy(fit_surv), type = "message"),"\n")
  
  cat("End: ", as.character(Sys.time()),"\n", sep = "")
}

# Output Writer:
STAN_diag_output.f <- function(fit,
                               modelName = "model",
                               output_file = "./Model_QuickDgnst.txt"){
  cat("Writing diagnostic data to: ",output_file, "\n", sep = "")
  cat("Writing diagnostics...")
  k <- getOption("width")
  options("width" = 250)
  
  sink(file = paste0(output_file))
  
  cat("STAN Model: ", modelName, ":\n", sep = "")
  
  cat("Model Info (chains, iterations, thinning:\n")
  print(fit@sim[2:4]) # print n-chains, n-iter, and n-thin
  
  cat("Quick diagnostics from STAN and M. Betancourt:\n")
  JMR_smry_diagnostics(fit, fitname = modelName)
  
  cat("\nSummary:\n")
  print(summary(fit)$summary)
  sink(NULL)
  
  options("width" = k)
  
  cat("done.\n")
  return(invisible(NULL))
  }

# ### Check Models ############################################################
# --- Survival ----------------------------------------------------------------
if(option_WriteOutputDiagn){
STAN_diag_output.f(fit = fit_surv,
                   modelName = "SampleSurvival",
                   output_file = paste0(output_loc,
                                        "ExmplSurv_QuickDiagn.txt"))
  }else{
    JMR_smry_diagnostics(STANfit = fit_surv, "SampleSurvival")
    }

# We are going to start with the survival model since the diagnostics
# are easier
# Number of iterations total (sum across all chains)
n_iter_a <- sum(fit_surv@sim$n_save - fit_surv@sim$warmup2)
# number of iterations for each chain
n_iter_c <- median(fit_surv@sim$n_save - fit_surv@sim$warmup2)

# Parameters to keep:
base_params <- c("mean_odds", "sd_odds")
use_params <- c(base_params,
                 "odds_rat[1]", "odds_rat[2]",
                 "odds_rat[10]","odds_rat[7]")


# ### Generate Plots ##########################################################
# plots are generated in a loop for each model fit supplied. There are two
# types of plots generated:
# 1. - simple diagnostics checking the validity of the MC. These are the same
# for all of the models regardless of their type and check things like
# divergences and sample sizes
# 2. - Posterior predictive checks (PostPC), these use the posterior as a model
# to simulate out data to compare it to the observed data. The point is to see
# how the posterior would predict the data, although MCMC results can be 
# 'valid' PostPC ensure that the model is appropriate for the data. The
# comparison is often done with summary statistics, chosen depending on the
# type of data, so different data types need different PostPC plots.

for(fit1 in list(fit_surv, fit_grow)){

# Data extracted
# In principle you only need to be able to extract the latter since]
# pivoting the data is possible, but this is just easier.
# Extract the draws, does not store information about the different chains
raw_MCMC_iter.df <- as.matrix(fit1) %>%
  as_tibble %>%
  mutate(iter = row_number(.),
         divergent = get_divergent_iterations(fit1))
# Same data but stored long with information from the chain
raw_MCMC_iter_wChain.df <- bayesplot::mcmc_trace_data(fit1)


# Which parameters are being used
# all names
par_names <- names(fit1)[!(grepl(pattern = "log_lik|g_rep",
                                 x = names(fit1)))]

  if(fit1@model_name == "exmple_Surv"){
    # survival will use the mean odds ratio and the variance in the odds ratio
    # and a few random species
  base_params <- c("mean_odds", "sd_odds")
  use_params <- c(base_params,
                "odds_rat[1]", "odds_rat[2]",
                "odds_rat[10]","odds_rat[7]")
  }else if(fit1@model_name == "exmpl_Growth"){
    # growth uses all the parameters in the upper level of the heirarchy
    base_params <- c("mu_growparam_CL",
                     "sigma_growparam_CL",
                     "logmu_growsd_CL",
                     "logsigma_growsd_CL")
    use_params <- c(base_params,
                    "growsd_sp_CL[1]", "growparam_sp_CL[2]",
                    "growsd_sp_CL[10]","growparam_sp_CL[7]")
  }
  
  # each plot is going in as a list, including the plot and the desired
  # output dimensions, at the end all of the plots will be generated and saved
  # into a single folder.
plt.L <- list()

# --- Common Diagnostic -------------------------------------------------------
# --- Rank order plot (chain) --- ####
# Pass: Chains are all exploring same location, all
# chains are ~uniform in distribution across the ranks.
plt.L <- list(
  mcmc_ranks.plt <-
    mcmc_rank_overlay(fit1, pars = use_params) +
    theme_bw() +
    # x axis
    scale_x_continuous(
      labels =label_percent(scale = 100/n_iter_a),
      n.breaks = 5,
      name = "Rank (percentile of iterations)") +
    # Y axis ( make labels easier to interpret)
    scale_y_continuous(labels = label_percent(scale = 100/n_iter_c),
                       n.breaks = 5,
                       limits = c(0,n_iter_c/10)) +
                ggtitle(paste0(fit1@model_name, "_rankHist")) +
                # description of the plot
                labs(caption = str_wrap(
                  paste("Histograms of rank order of parameters/quanitites",
                        "if mixing properly, rank order should be uniformly",
                        "distributed for all chains (flat). ",
                        "Look-out for sloped or uneven histograms."))) +
                theme(plot.caption = element_text(hjust = 0)),
              # dimensions
              width = 11.75,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)

# print the name to keep track
print(plt.L[[1]][[1]]$labels$title)

# --- Trace plot --- ####
plt.L <- list(mcmc_trace_data(x = fit1, pars = use_params) %>%
                ggplot(., aes(x = iteration, y = value, color = chain)) +
                geom_line() +
                # reformat x-axis
                scale_x_continuous(labels = function(x) format(x, scientific = TRUE),
                                   n.breaks = 4) + 
                scale_color_viridis_d(alpha = 0.5, option = "G") +
                facet_wrap("parameter", scales = "free") +
                theme_bw() +
                # title
                ggtitle(paste0(fit1@model_name, "_trace")) +
                # caption to aid interpretation
                labs(caption = str_wrap(
                  paste("Trace-plot of key quantities. We should be seeing",
                        "consistent and overlaping chains.",
                        "Look-out for sloped chains, (burn-in not complete/slow mixing),",
                        "clear tracing of individual chains (slow mixing),",
                        "or noticable periods of horizontal values within chain ",
                        "(proposals too dramatic).",
                        "This plot isn't terrible informative."))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 11.75,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# --- Autocorrelation Lag plot --- ####
# Autocorrelation vs Lag plot, make sure it rapidly approaches 0.
plt.L <- list(mcmc_acf.plt <- mcmc_acf(x = fit1,
                                       pars = use_params,
                                       lags = 25)$data %>%
                mutate(Chain = as.factor(Chain)) %>%
                ggplot(data = ., aes(x = Lag, y = AC, color = Chain)) +
                geom_line() +
                # facet wrap by the parameter
                facet_wrap(facets = "Parameter") +
                scale_y_continuous(name = "Autocorrelation") +
                scale_color_viridis_d(alpha = 0.5, option = "G") +
                theme_bw() +
                # add title
                ggtitle(paste0(fit1@model_name, "_lagAcf")) +
                # add caption to aid interpretation
                labs(caption = str_wrap(
                  paste("Plot showing non-independence among draws",
                        "for different quantities.",
                        "Seeing some autocorrelation is normal;",
                        "however, it should decay rapidly within",
                        "<10 draws."))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 11.75,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# --- Cummulative Averages --- ####
# plot showing the running average across multiple chains if everything
# is going well we shouldn't see a trend in the average value for the
# quantities 
rol_avg_n <- 50
plt.L <- list(mcmc_RolAvg.plt <-
                # Get the raw chains (needed for rolling average)
                raw_MCMC_iter_wChain.df %>%
                arrange(iteration) %>%
                # calculate the rolling average
                mutate(.by = c(parameter, chain),
                       w = 1,
                       roll_avg = cumsum(value)/cumsum(w),
                       roll_max = cummax(value),
                       move_avg = rollmean(x = value,
                                           k = rol_avg_n,
                                           fill = NA,
                                           align = "right")) %>%
                select(-c(w)) %>%
                filter(parameter %in% use_params) %>%
                # Start plotting
                ggplot(data = .,
                       aes(x = iteration,
                           y = move_avg,
                           color = chain)) +
                geom_line() + 
                scale_color_viridis_d(alpha = 0.5, option = "G") +
                # plot each parameter separately
                facet_wrap("parameter") +
                scale_x_continuous(n.breaks = 4) +
                scale_y_continuous(name = "Moving average") +
                theme_bw() + 
                ggtitle(label = paste(fit1@model_name, "RollAvg", sep = "_")) +
                # caption for adding information to aid interpretation
              labs(caption = str_wrap(
                paste("Plot rolling average (", rol_avg_n,
                      " vals) of different quantities ",
                      "We shouldn't be seeing any obvious trends. ",
                      "Note the difference in scales for different quantities",
                      sep = ""))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 11.75,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# --- Effective Sample sizes --- ####
# This plots the effective sample sizes for estimates of the bulk and tail
# of the posterior. Tails can often act up so hence why we need them to be
# looked at separately
flag_ESSTail <- flag_ESSBulk <- 1000 # recommended cutoffs
min_ESSTail <- min_ESSBulk <- 100 # absolute minimum cutoffs

plt.L <- list(mcmc_NRatioRhat.plt <-
                # setup the data 
                data.frame(param = par_names, # parameter names
                           # bulk effect. sample size for each param
                           ess_bulk = apply(as.array(fit1, pars = par_names),
                                            MARGIN = 3,
                                            FUN = ess_bulk),
                           # tail effective sample size for each param
                           ess_tail = apply(as.array(fit1, pars = par_names),
                                            MARGIN = 3,
                                            FUN = ess_tail)) %>%
                # flag quantities that are below the recommended
                mutate(label = if_else(ess_bulk < flag_ESSBulk |
                                         ess_tail < flag_ESSTail,
                                       param, "")) %>%
                ggplot(., aes(x = ess_bulk, y = ess_tail)) +
                geom_point() +
                # label any past thresholds
                geom_text(aes(label = label), hjust = 0, vjust = -.75) +
                # axes
                scale_x_continuous(limits = c(NA, NA)) +
                scale_y_continuous(limits = c(NA, NA)) +
                # lines for thresholds, dashed is a "concern", solid is bad
                # dashed is the recommended
                geom_vline(xintercept = c(min_ESSBulk, flag_ESSBulk),
                           linetype = c(1,2),
                           color = "navy") +
                geom_hline(yintercept =  c(min_ESSTail, flag_ESSTail),
                           linetype = c(1,2),
                           color = "navy") +
                theme_bw() +
                # plot title
                ggtitle(label = paste(fit1@model_name, "ESS", sep = "_")) +
                # caption for plot context
                labs(caption = str_wrap(
                  paste("Plot effective sample size in the main",
                        "'body' (bulk), and in the tails for all quantities",
                        "Solid lines are the absolute minimums, while dashed",
                        "are recommended. Plot automatically labels points",
                        "that are flagged for another look.",
                        sep = " "))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 7.5,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# --- Rhat vs N_ratio --- ####
# rhat and Nrat, plots the gelman-Rubin statistic R-hat against the 
# ratio of effective sample size to the number of kept draws
# Gelman-Rubin Rhat guidlines, above 1.1 is bad, keep an eye if above 1.05
max_Rhat <- 1.1
flag_Rhat <- 1.05
# Ratio of sample sizes, minimum is 1:10 Neff/N
# recommended is 1:2
min_Nrat <- 0.1
flag_Nrat <- 0.5


plt.L <- list(ESS.plt <-
                # setup the data
                data.frame(param = par_names, # parameter names to use
                           # ratio of Neff/N 
                           Nratio = neff_ratio(fit1, pars = par_names),
                           # R-hat statistic
                           Rhat = rhat(fit1, pars = par_names)) %>% 
                mutate(label = if_else(Rhat > flag_Rhat | Nratio <= flag_Nrat,
                                       param, "")) %>%
                ggplot(., aes(x = Nratio, y = Rhat)) +
                geom_point() +
                # label any past thresholds
                geom_text(aes(label = label), hjust = 0, vjust = -.75) +
                # axes
                scale_x_continuous(limits = c(0, NA), name = "Neff/N") +
                scale_y_continuous(limits = c(NA, NA)) +
                # lines for thresholds, dashed is a "concern", solid is bad
                geom_vline(xintercept = c(min_Nrat, flag_Nrat),
                           linetype = c(1,2),
                           color = "navy") +
                geom_hline(yintercept =  c(flag_Rhat,max_Rhat),
                           linetype = c(2,1),
                           color = "navy") +
                theme_bw() +
                # title
                ggtitle(label = paste(fit1@model_name,
                                      "NRat_v_Rhat", sep = "_"))+
                # caption for the plot to aid interpretation
                labs(caption = str_wrap(
                  paste("Plot ratio of effective sample size of MC (N_eff):",
                        "the number of draws against Gelman-Rubin R-hat.",
                        "Solid lines are the absolute minimums, while dashed",
                        "are recommended. Plot automatically labels points",
                        "that are flagged for concern. Low n-ratios show",
                        "that chain is too correlated, high R-hat means",
                        "that different chains are exploring different areas.",
                        sep = " "))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 7.5,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)


# --- Energy Distributions --- ####
# Hamiltonian specific diagnostic useful for seeing if the sample is
# behaving well. We should see that the energy transfer distribution pi_dE 
# is roughly matching the marginal energy distribution (pi_E).
plt.L <- list(mcmc_EDistr.plt <-
                # call plotting
                mcmc_nuts_energy(nuts_params(fit1),
                                 bayesplot::log_posterior(fit1),
                                 binwidth = 1) +
                theme(plot.background = element_rect(fill = "white"))+
                ggtitle(label = paste(fit1@model_name, "EDistr", sep = "_")) +
                # label the plot caption to aid interpretation
              labs(caption = str_wrap(
                paste("Energy transition density (delta_E) vs the marginal",
                      "energy distribution (E). They should be similar.",
                      "If the energy transition density (delta_E) is notably",
                      "'narrower' than the marginal energy distribution, then",
                      "the chain main not be properly sampling the tails of",
                      "the distribution",
                      sep = " "))) +
                theme(plot.caption = element_text(hjust = 0)),
              width = 11.75,
              height = 7.5) %>%
  list() %>%
  c(., plt.L)
print(plt.L[[1]][[1]]$labels$title)

# --- Pairs plot --- ####
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
print(plt.L[[1]][[1]]$labels$title <-
        paste0(fit1@model_name, "_pairs")) # can't get title onto plot directly

# --- Parallel Coordinates plot --- ####
# parallel coordinated plot 
n_other <- 1000 # number of draws to keep that AREN'T divergences

# prepare the data
parc_dat <- mcmc_parcoord_data(
  x = as.array(fit1, pars = c(use_params)),
  pars = c(use_params),
  np = nuts_params(fit1))
#parc_dat$Divergent[which(parc_dat$Draw == 1)] <- 1

# checks if there are any divergences
div_draws <- summarize(parc_dat,
                       .by = "Draw",
                       Divergent = sum(Divergent)) %>%
  filter(Divergent != 0) %>%
  pull(Draw)
# keeps all divergences plus a subsample of other draws
keep_draws <- sample(with(filter(parc_dat,
                                 !(Draw %in% div_draws)),
                          unique(Draw)),
                     size =  n_other,
                     replace = F) %>%
  c(., div_draws)

# Add the levels/recodes them to make sure they plot properly
parc_dat <- parc_dat %>%
  mutate(Divergent = factor(Divergent, ordered = T))
parc_dat$Divergent <- factor(parc_dat$Divergent, levels = c("1","0"))

# generate the plot
plt.L <- list(mcmc_parcoord.plt <-
                parc_dat %>%
                # subset to remove most of the normal draws.
                filter(Draw %in% keep_draws) %>%
                mutate(Draw = Draw*(-1)) %>% # ensures divergences plot last
                arrange(., desc(Divergent)) %>%
                ggplot(., aes(x = Parameter,
                              y = Value,
                              group = as.factor(Draw),
                              color = Divergent)) +
                geom_line(alpha = .75) +
                scale_color_manual(values = c("navy", "red")) +
                ggtitle(label = paste(fit1@model_name, "ParCoord", sep = "_")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              # add caption to improve interpretation
              labs(caption = str_wrap(
                paste("Parallel coordinates plot. Use this to 'follow'",
                      "divergences, to see if divergences are caused by",
                      "particular values/combinations of parameter(s).",
                      "Divergences are highlighted in red if present.",
                      "Plot only displays", n_other,"draws that are",
                      "not-divergences to prevent overcrowding in the",
                      "plotting window.",
                      sep = " "))) +
                theme(plot.caption = element_text(hjust = 0)),
              height = 7.5, width = 7.5) %>%
  list %>%
  c(., plt.L)

print(plt.L[[1]][[1]]$labels$title)
# --- Posterior predictive check ----------------------------------------------
cat("Generating PostPC plots...\n\n")
n_ysim <- 100 # number of simulations for posterior pc density plots

# --- For the binary variables ------------------------------------------------
if(fit1@model_name == "exmpl_Surv"){
  # extract the proportions of surving 
  post_p.m <- as.matrix(fit1, pars = "p")
  
  # Recover the initial values used to fit the model
  time.v <- exmpl_surv.df$t
  spp.v <- exmpl_surv.df$species
  n_spp <- length(unique(spp.v))
  n_stem <- nrow(exmpl_surv.df)
  
  
  # vector for the realized p
  p_rlz.v <- vector(mode = "numeric", length = n_stem)
  
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
  
  # generate the y_sims (simulate the observations)
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
  dat <- data.frame(spp = spp.v,
                    bindat = exmpl_surv.df$survival) %>%
    summarise(p = sum(bindat)/n(),.by = spp)
  
  # function for retrieving the proportion
  ppn <- function(x){return(sum(x)/length(x))}
  
  # --- All species proportion surviving --- ####
  plt.L <- list(ppc_all.plt <-
                  # generate plot
                  ppc_dens_overlay(y = dat$p, yrep = post_pc1.m) +
                  # formatting (title, labels, etc.)
                  ggtitle(label = paste0(fit1@model_name, "_PostPC_distrPpn")) +
                  ylab("Density") +
                  theme(plot.background = element_rect(fill = "white"))+
                  # caption to aid in interpretation
                  labs(caption = str_wrap(
                    paste("Posterior predictive check of the ",
                          "density of observed proportions of surviving stems ",
                          "across all species",
                          "(n_sims = ", n_ysim, "). ",
                          "We should see that our true values should lie within ",
                          "those predicted with the fitted model ",
                          "otherwise, the model is likely inappropriate.",
                          sep = ""))) +
                  theme(plot.caption = element_text(hjust = 0)) +
                  xlab("p"), width =7.5, height = 7.5) %>%
    list() %>% c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  # --- Proportion surviving by species --- ####
  plt.L <- list(ppc_grouped.plt <-
                  ppc_stat_grouped(yrep = post_pc2.m, # simulated data matrix
                                   y = exmpl_surv.df$survival, # true survival obs
                                   group = spp.v, # grouping (by species)
                                   stat = ppn, # calculate proportion surviving
                                   freq = F) +
                  theme(plot.background = element_rect(fill = "white"))+
                  ggtitle(label = paste0(fit1@model_name, "_PostPC_ppnBySpp")) +
                  scale_x_continuous(n.breaks = 3) +
                  # add caption to plot to help interpretation
                  labs(caption = str_wrap(
                    paste("Posterior predictive by check of the proprotion of",
                          "surviving stems by species.",
                          "Histograms showing simulated predicted",
                          "number of stems surviving",
                          sep = " "))) +
                  theme(plot.caption = element_text(hjust = 0)),
                height = 12, width = 16) %>%
    list() %>% c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  
}else{
  # --- Lognormal PostPC ------------------------------------------------------
  # Posterior predictive check for the growth data. Using skew and kurtosis
  # since these aren't directly used in the parametrization.
  dat <- data.frame(species = exmpl_grow.df$species,
                    vals = exmpl_grow.df$grow_val)
  
  # Extracting simulations
  # two different ways since different plotting functions use different formats
  posterior1 <- rstan::extract(fit1) # get the posterior of the fit
  posterior_array1 <- as.array(fit1)
  g_reps <- fit1 %>%
    as.matrix(., pars = "g_rep") # reps
  rep_rows <- sample(nrow(g_reps),
                     size = 100,
                     replace = F)
  
  
  # --- Distribution of obs --- ####
  # PostPC for all observations to see how they map out.
  plt.L <- list(PostPCAllVals.plt <-
                  # generating plot
                  ppc_dens_overlay(dat$vals, g_reps[rep_rows,]) +
                  # modify axes
                  scale_y_continuous(limits = c(0, 1)) +
                  ylab("Density") +
                  xlab("Value") +
                  theme(plot.background = element_rect(fill = "white"))+
                  # title
                  ggtitle(label = paste(fit1@model_name, "PostPCAllVals", sep = "_")) +
                  # adding caption
                  labs(caption = str_wrap(
                    paste("Posterior predictive of growth across all species",
                          "true values should be within simulated data, or",
                          "the model could be missing key features/misspecified.",
                          sep = " "))) +
                  theme(plot.caption = element_text(hjust = 0)),
                
                width = 7.5,
                height = 7.5) %>%
    list %>%
    c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  # --- skew PostPC --- ####
  # Post PC using the skew of the obs, broken up by species
  plt.L <- list(PostPCSkew.plt <- ppc_stat_grouped(dat$vals,
                                                   g_reps, 
                                                   stat = skewness,
                                                   group = dat$species) +
                  scale_x_continuous(n.breaks = 3) +
                  theme(plot.background = element_rect(fill = "white"))+
                  ggtitle(paste(fit1@model_name, "PostPCSkew", sep ="_")) +
                  # add caption
                  labs(caption = str_wrap(
                    paste("Posterior predictive of skew in growth across species",
                          "true values should be within simulated data, or",
                          "the model could be missing key features.",
                          sep = " "))) +
                  theme(plot.caption = element_text(hjust = 0)),
                width = 16,
                height = 12) %>%
    list %>%
    c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
  
  # --- Kurtosis PostPC --- ####
  # post pc using the kurtosis, broken up by species
  plt.L <- list(PostPCKurt.plt <- ppc_stat_grouped(dat$vals,
                                                   g_reps, 
                                                   stat = kurtosis,
                                                   group = dat$species) +
                  scale_x_continuous(n.breaks = 3) +
                  theme(plot.background = element_rect(fill = "white"))+
                  ggtitle(paste(fit1@model_name,"PostPCKurt", sep ="_")) +
                  # adding caption
                  labs(caption = str_wrap(
                    paste("Posterior predictive of kurtosis in growth across species",
                          "true values should be within simulated data, or",
                          "the model could be missing key features.",
                          sep = " "))) +
                  theme(plot.caption = element_text(hjust = 0)),
                width = 16,
                height = 12) %>%
    list %>%
    c(., plt.L)
  print(plt.L[[1]][[1]]$labels$title)
}
# --- Output plots ----------------------------------------------------------
cat("Outputing plots...")
# location for output
# do a separate folder for each fit
loc_out_model <- paste0(output_loc,"/",fit1@model_name,"/")

# loop through the plots, the list stores the plot, and its dimensions
for(plt_i.L in plt.L){
  plt <- plt_i.L[[1]] # extract the plot
  ttl <-  plt$labels$title # extract plot title
  # save the plot, use the title for filename
  ggsave(plot = plt,
         filename = paste0(loc_out_model,
                                 ttl, ".png"),
         device ="png", width = plt_i.L$width,
         height = plt_i.L$height)
} # end of plotting loop

} # end loop across fits