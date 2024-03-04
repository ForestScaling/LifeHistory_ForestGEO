# File 03_STANSetup.R
# Author: JMR
# Last Update: 12023-05-22

# ### Setup ###################################################################
# Setting the seed:
set.seed(123)

# location of the data
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")

# set the directory of the scripts
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Bringing in the packages:
for(x in c("bayesplot", "rstan",
           "gdata", "magrittr",
           "rstudioapi", "tidyverse")){
  library(x, character.only = T)
}

# STAN options/defaults for the script
options(mc.cores = min(1, parallel::detectCores() - 2)) # cores JMR mchn w/ 12
rstan_options(auto_write = TRUE,
              javascript = FALSE) # Don't recompile the code if unmodified

# Loading in diagnostics written by M. Betancourt
source('scripts/STAN/stan_utility_rstan.R')

# --- Site specific setup -----------------------------------------------------
# Load in the site's data
site <- "HVDF"
load(paste0(loc_data, site, "_AllStems_Vals.r"))
AllData.L <- get(paste0(site, "_AllStems_Vals.L"))

# Bring out the values
taxa.v <- AllData.L$taxa
stature.df <- AllData.L$stature
surv.L <- AllData.L$survival
grow.L <- AllData.L$growth
n_census <- AllData.L$n_census
recrt.L <- AllData.L$recruitment
n_stem <- 

common10 <- stature.df %>%
  arrange(desc(number)) %>%
  select(Mnemonic) %>%
  mutate(Mnemonic = as.character(Mnemonic))

common10 <- common10$Mnemonic[1:10]


# --- Growth Data -------------------------------------------------------------
# Growth 2
# removing all NA values
 grow.L %<>% lapply(., function(x){
   cnsus <- unique(x$Census) 
   n_stem_ini <- nrow(x)
   x %<>% filter(!is.na(growth2))
   
   cat(paste0("For: ", site, "-", cnsus,", ", nrow(x),
              " have growth values (",round(100*nrow(x)/n_stem_ini, 1),"%)."))
   return(x)
   
 })

# Starting with a single year:
grow.df <- grow.L[[1]] %>% mutate(., CanopyLvl = as.character(CanopyLvl))

# plot breaking down growth by the most common species
plot.df <- grow.df %>% filter(Mnemonic %in% common10)
ggplot(plot.df, aes(x = growth2, colour = Mnemonic)) +
  geom_density() + scale_x_continuous() +
  ggtitle(label = paste0(site, " Growth-II by Species")) +
  theme_bw()

# plot breaking down growth rate by canopy position
plot.df <- grow.df
ggplot(data = plot.df, aes(x = growth2, colour = CanopyLvl)) +
  geom_density() +
  ggtitle(label = paste0(site, " Growth-II by Canopy Level")) +
  theme_bw()


# prep the data for some simple coomparisons with the "best data"
file_STAN <- "scripts/STAN/growth2_lognormal_diag_byCL.stan"
stan.df <- grow.df %>% filter(!is.na(CanopyLvl) & !is.na(growth2))
#stan_data.L <- list(n_stem = nrow(stan.df),
#                  n_lvl = length(unique(stan.df$CanopyLvl)),
##                  levels = as.integer(stan.df$CanopyLvl),
#                  growth2 = stan.df$growth2)

# Check that the file is ready to go
stanc(file_STAN)$status

# Run the MCMC
stan_dur <- system.time(
  fit1 <- stan(file = file_STAN,
               data = stan_data.L,
               warmup = 1000, # 500
               iter = 2000,
               chains = 5,
               cores = 5,
               thin = 1,
               seed = 123,
               control = list(adapt_delta = 0.90))
)



fake_data <- stan.df %>%
  mutate(., CanopyLvl = case_when(CanopyLvl == "4" ~ "3",
                                 TRUE ~ CanopyLvl)) %>%
  filter(Mnemonic %in% c("tsugca","acerru","querru",
                         "kalmla","betule","pinust"))


STAN_Lognormal_res.L <- list(1, 2, 3)

STAN_Lognormal_taxa.L <- list(1,2,3)
for(l in 1:3){
  
  print(l)
  print(timestamp())
  use.df <- fake_data %>% filter(CanopyLvl == as.character(l))
  
  STAN_Lognormal_taxa.L[[l]] <- use.df %>%
    group_by(., Mnemonic) %>%
    summarize(., n_obs = n()) %>%
    mutate(., spp_numb = 1:nrow(.))
  
  use.df %<>% left_join(., STAN_Lognormal_taxa.L[[l]], by = "Mnemonic")
  
  stan_data.L <- list(n_stem = nrow(use.df),
                      n_spp = length(unique(use.df$spp_numb)),
                      species = use.df$spp_numb,
                      growth2 = use.df$growth2)
  
  file_STAN <- "scripts/STAN/growth2_lognormal_CLxSpp2.stan"
  
  # Run the MCMC
  STAN_Lognormal_res.L[[l]]  <- stan(file = file_STAN,
                 data = stan_data.L,
                 warmup = 1000, # 500
                 iter = 2000,
                 chains = 5,
                 cores = 5,
                 thin = 1,
                 seed = 123)
                 #control = list(adapt_delta = 0.90))
  
  
  
}


traceplot(STAN_Lognormal_res.L[[1]], pats)


CL1 <- rstan::extract(STAN_Lognormal_res.L[[1]])
plot(density(CL1$sp_growth[,1]), main = expression(beta[1]))

STAN_Lognormal_taxa.L[[1]]









