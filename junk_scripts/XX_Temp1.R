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
site <- "SERC"


param_2_disregard <- c("g_rep", "log_lik", "lp__")
growL1data.df <- data.frame()

for(site in c("HVDF", "SCBI", "SERC")){
  load(paste0(loc_data, site, "_STANdata.r")) # load the data used (for PPC)
  load(paste0(loc_stanout, site, "_STANData_MCMC.r")) # load the data used
  load(paste0(loc_stanout, site, "_STAN_ResultsMCMC.R")) # load the fits
  load(paste0(loc_data2, site, "_TaxaTable.r"))
  TaxaTable <- Spp_table
  k<-1
  param_fit <- All_STAN_fits.L[[k]]
  fit_name <- names(All_STAN_fits.L)[k]
  
  # get the temporary ids to actually get the proper Mnemonic ready
  tempID.df <- All_params_SpeciesID.df %>%
    filter(param == fit_name) %>%
    mutate(tempID = as.character(tempID)) %>%
    left_join(., select(TaxaTable, c(Mnemonic, Latin)),
              by = "Mnemonic")
  
  # get the fit of the particular parameter
  

 
  
  
  # Names of the parameters of interest
  raw_param_names <- names(param_fit)[
    # only keep those parameters of interests (e.g. don't keep log_lik)
    !grepl(pattern = paste(param_2_disregard, collapse = "|"),
           x = names(param_fit))]  
  
  # covert the draws into a data.frame
  draws.df <- param_fit %>%
    as.matrix(.) %>%
    as.data.frame %>%
    select(any_of(raw_param_names))
  
  # Pull out only growth mean vals
  draws.df <- draws.df %>%
    select(contains("growparam")|contains("odds")) %>%
    select(!contains("sigma")) %>%
    select(!contains("sd_odds"))
  
  centers <- draws.df %>%
    select(contains("mean")|contains("mu_growparam_CL")) %>%
    unlist()
  
  # normal values being drawn
  reg_growDraws.df <- draws.df %>%
    pivot_longer(cols = everything(),
                 names_to = "param",
                 values_to = "val") %>%
    mutate(tempID = str_extract(string = param, pattern = "\\[([^()]+)\\]"),
           tempID = str_sub(tempID, start = 2,-2)) %>%
    select(-param) %>%
    left_join(., tempID.df, by = "tempID")
  
  
  
  # creating a column that has the comparison to the mean value
  # to get the distribution to the "intercept" of the model. 
  # can be used to compare the relative growth to other species at the site
  dlt_growDraws.df <- draws.df %>%
    mutate(across(everything(), ~ .x + centers)) %>% 
    pivot_longer(cols = everything(),
                 names_to = "param",
                 values_to = "val") %>%
    mutate(tempID = str_extract(string = param, pattern = "\\[([^()]+)\\]"),
           tempID = str_sub(tempID, start = 2,-2)) %>%
    select(-param) %>%
    left_join(., tempID.df, by = "tempID")
  
  
    
  # names are backwards, don't overthink it
  growDraws.df <- reg_growDraws.df %>%
    mutate(type = "delta") %>%
    bind_rows(mutate(dlt_growDraws.df, type = "reg"))
  
  growL1data.df <- growDraws.df %>%
    mutate(site = site) %>%
    rbind(growL1data.df, .)
}


# species at all site
summarize(growL1data.df, .by = "Latin", n_d)
  
growL1data.df %>%
  filter(type == "reg",) %>%
  
  ggplot(., aes(x = val, y = Latin, fill = site)) +
  geom_density_ridges(alpha = .5) +
  ggtitle(label = "mean CL1 odds_survival values (reg)")






  ggplot(reg_growDraws.df, aes(x = val, y = param)) +
    geom_density_ridges() +
    ggtitle()
  
  
  
  
  
  # extract the species code
  
  
  
  
  draws.df <- All_params_SpeciesID.df %>%
    filter(param == param_fit) %>%
    
  
  # add in the species specific data
  
  
  
  
}

# ### Quick correlations of posteriors
# g1 <- dlt_growDraws.df
# s1 <- dlt_growDraws.df

s1.m <- s1 %>% select(c(val, Latin))
g1.m <- g1 %>% select(c(val, Latin))

#g1.df <- rstan::extract(param_fit, permuted = T)
g1.df <- g1.df$growparam_sp_CL
colnames(g1.df) <- tempID.df %>%
  arrange(desc(tempID)) %>%
  pull(Mnemonic)

#s1.df <- rstan::extract(param_fit, permuted = T) 
s1.df <- s1.df$odds_rat
colnames(s1.df) <- tempID.df %>%
  arrange(desc(tempID)) %>%
  pull(Mnemonic)

# combine the two sets of data
# We need to filter out the taxa that aren't estimated in both
s1.df <- s1.df %>% as.data.frame() %>% select(any_of(colnames(g1.df)))
g1.df <- g1.df %>% as.data.frame() %>% select(any_of(colnames(s1.df)))

# We need to make sure that column names are aligned properly
g1.m <- g1.df[,sort(colnames(g1.df))] %>% t()
s1.m <- s1.df[,sort(colnames(s1.df))] %>% t()

# verify matching of rownames
identical(rownames(g1.m), rownames(s1.m))

ary.a <- array(data = NA, dim = c(dim(g1.m),2))
ary.a[,,1] <- g1.m
ary.a[,,2] <- s1.m
ary.L <- lapply(seq(from = 1, to = dim(ary.a)[2], by = 1),
                function(x){
                  ary.a[,x,]
                })
cor.v <- lapply(ary.L, FUN = function(x){cor(x = x[,1], x[,2])}) %>% unlist
hist(cor.v, main = site)





m1 <- matrix(1:15, nrow = 3); m2 <- m1*10
a1 <- array(data = NA, dim = c(dim(m1), 2))
a1[,,1] <- m1; a1[,,2] <- m2
ary.a
a.L <- lapply(seq(from = 1, to = dim(a1)[2], by = 1),
              function(x){
                a1[,x,]
              })

lapply(a.L, FUN = function(x){cor(x = x[,1], x[,2])})


