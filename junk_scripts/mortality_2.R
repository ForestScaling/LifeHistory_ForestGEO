# Checking simple correlations
# Running data with simple correlations, varying the methods to see how they 
# impact results, correlating growth and survival rates, using simple
# bootstraps of the data to asses CI.

# ### Setup ###################################################################
# Options
FILTER_DBH.lgl <- T       # logical if to subset for min stem size
DBH_min <- 10             # cutoff in cm for minimum stem size
FILTER_MAINSTEM.lgl <- F  # Logical. if True only use main stem data.
ppn_grow.Lgl <- F         # use proportional growth, if F, use linear
grow_logTransf.Lgl <- F   # log transform growth values?
min_nobs  <- 5            # min number of obs for taxon inclusion 0 for none
surv_logit_transf <- F    # Logit transform survival
CI_a <- 0.05              # alpha for bootstrap Conf Intervals
n_bstps <- 2000           # number of bootstraps

# locations
if(!file.exists("./scripts/AA_LocationManager.R")){setwd("..")}
source("./scripts/AA_LocationManager.R")

# Data 
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")
loc_data1 <- paste0(loc_Gdr, "/data/ForestGEO/processed/canopy_assign/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/taxa/")
loc_out1 <- paste0(loc_Gdr, "/outputs/")

# Bringing in the packages:
for(x in c("magrittr", "rstudioapi", "tidyverse",
           "rlang")){
  library(x, character.only = T)
}

# sites to use
site.v <- c("HVDF", "SCBI", "SERC")
data.L <- vector("list",  length(site.v))
names(data.L) <- site.v

# Loading in the data
for(site in site.v){
  # Load in the data
  load(paste0(loc_data, site, "_STANdata.r"))
  load(paste0(loc_data1, paste0(site, "_AllStemsCL.r")))
  
  data.L[[site]] <- AllStems
}
 
# process the data into a dataframe
for(site in site.v){
  # First census
  AllStems <- data.L[[site]]
  Census.df <- AllStems[[1]] %>%
    select(c(TreeID, StemID, Mnemonic, Date, Stem, 
             DBH, CanopyLvl, RametStatus))
  
  # Second Census
  laterCensus.df <- AllStems[[2]] %>%
    select(c(TreeID, StemID, Mnemonic, Date, 
             Chng_HOM, DBH, CanopyLvl, RametStatus)) %>%
    rename_with(.fn = ~paste0(.x, "_t2"),
                .cols = -c(TreeID, StemID,
                           Mnemonic, Chng_HOM))
  
  # Group the censuses
  Census.df <- Census.df %>%
    left_join(laterCensus.df, by = c("TreeID", "StemID", "Mnemonic")) %>%
    as_tibble(.)
  
  # --- Growth Data -------------------------------------------------------------
  # Filter the data
  growth.df <- Census.df %>%
    filter(RametStatus == "A" & RametStatus_t2 == "A") %>% # must be alive t1 & t2
    filter(!FILTER_MAINSTEM.lgl | # verify the main stem
             Stem == "Main" |
             is.na(Stem)) %>%
    filter(!is.na(DBH) & !is.na(DBH_t2)) %>% # DBH can't be NA
    filter(!FILTER_DBH.lgl | DBH >= DBH_min) %>% # DBH above minimum
    filter(!as.logical(Chng_HOM)) # can't change position of DBH measurement
  
  # Calculate growth
  growth.df <- growth.df %>%
    mutate(across(contains("Date"), ~ymd(.x))) %>% # convert Date format
    mutate(diff_yr = Date_t2 - Date, .before = 1) %>% # get diff in observ
    mutate(diff_yr = diff_yr/dyears(1)) %>% # change duration into years
    mutate(linear_growth = (DBH_t2 - DBH)/diff_yr, # linear growth
           propn_growth = (DBH_t2/DBH)^(1/diff_yr)) # proportional growth
  
  
  # set the variable for growth
  grow_var <- ifelse(ppn_grow.Lgl, "propn_growth", "linear_growth")
  
  # summarize by species
  growth.df <- growth.df %>%
    rename(growth = !!grow_var) %>%
    summarise(., .by = c(Mnemonic),
              n_grow = n(),
              mean_grow = mean(growth)) %>%
    filter(n_grow >= min_nobs)
  
  if(grow_logTransf.Lgl){
    growth.df <- growth.df %>%
      mutate(mean_grow =log(mean_grow))
  }
    
    
  # --- Mortality Data ------------------------
  mortality.df <- Census.df %>%
    filter(RametStatus == "A") %>% # must be alive t1 & t2
    filter(!FILTER_MAINSTEM.lgl | # verify the main stem
             Stem == "Main" |
             is.na(Stem)) %>%
    filter(!FILTER_DBH.lgl | DBH >= DBH_min) # DBH above minimum
    
  # calculate death
  mortality.df <- mortality.df %>%
    mutate(survived = case_when(RametStatus_t2 == "D" ~ 0,
                                RametStatus_t2 == "A" ~ 1,
                                TRUE ~ NA)) %>%
    filter(!is.na(survived)) # remove NA values (missing at t2)
  
  # count number of survivors
  mortality.df <- mortality.df %>%
    summarise(.by = Mnemonic,
              n_surv = sum(survived),
              n_Total = n()) %>%
    mutate(ppn_surv = n_surv/n_Total) %>%
    filter(n_Total >= min_nobs)
  
  # Logit transforming survival
  if(surv_logit_transf){
    mortality.df <- mortality.df %>%
      mutate(ppn_surv = log(ppn_surv/(1-ppn_surv))) %>%
      filter(!(is.na(ppn_surv)|is.infinite(ppn_surv) ))
    
  }
  
  # --- Join the data -----------------------------------------------------------
  data_all.df <- mortality.df %>%
    inner_join(., growth.df, by = "Mnemonic") %>%
    filter(!is.na(ppn_surv) & !is.na(mean_grow))
  
  
  data.L[[site]] <- data_all.df %>%
    mutate(site_loc = site)

}

# turn into a single dataframe
data.df <- bind_rows(data.L)

# extract site specific coefficients for correlation lines
coef.df <- data.frame(site = site.v, b0 = NA,
                      b1 = NA, LCL = NA,
                      UCL = NA, r = NA,
                      n = NA)

# bootstraps for cor CI
bootstrap.f <- function(df){
  df <- df[sample(1:nrow(df),
                  size = nrow(df),
                  replace = T), ]
  return(df)
}

# loop through the sites
var_y <- "ppn_surv"
var_x <- "mean_grow"

for(site in site.v){
  # Calcualte coefficients for the line
  coef.df[which(coef.df$site == site), c("b0", "b1")] <-
    data.df %>%
    select(c(site_loc, !!var_x, !!var_y)) %>%
    filter(site_loc == site) %>%
    #map_df(.f = scale) %>%  # scale the data (at the level of each col)
    as.data.frame() %>% # fixes names of columns
    lm(data = ., formula =  formula(paste(var_y, "~", var_x))) %>%
    coef() # extract the coefficients
  
  # get pearsons's r
  coef.df[which(coef.df$site == site), "r"] <-
    cor(data.df[which(data.df$site_loc == site),var_y],
        data.df[which(data.df$site_loc == site),var_x])
  
  # n
  coef.df[which(coef.df$site == site), "n"] <- data.df %>%
    filter(site_loc == site) %>% nrow()
  
  # Bootstrap the CI for pearson's r
  coef.df[which(coef.df$site == site),
          c("LCL", "UCL")] <-
    map(1:n_bstps, ~bootstrap.f(df = data_all.df)) %>%
    sapply(., FUN = function(X){
      return(cor(x = X[[var_y]], y = X[[var_x]]))
      }) %>%
    quantile(., probs = c(CI_a/2, 1 - CI_a/2)) %>%
    sprintf("%.2f", .)

}



# Generating plot
# aggregating text for the CI label
CI <- paste0(coef.df$site, ", r = ",
             sprintf("%.2f",coef.df$r),
             " (",coef.df[,"LCL"],", ",
             coef.df[,"UCL"], ")",
             " n = ", coef.df[,"n"]) %>%
  paste(., collapse = "\n")
  
# Plot for corr between variables
cor.plt <-
  data.df %>% 
  ggplot(mapping = aes(x = mean_grow,
                       y = ppn_surv,
                       colour = site_loc,
                       group = as.factor(site_loc))) +
  geom_point(alpha = 1) +
  # add lines to plot to show linear correlation
  geom_abline(data = coef.df,
              aes(slope = b1,
                  intercept = b0,
                  color = site)) +
  # Add the graph's axis labels
  annotate("label",
           x = min(data.df[,var_x]),
           y = min(data.df[,var_y]),
           label = CI,
           hjust = 0,
           vjust = 0) +
  # Labelling the caption of the plot
  labs(caption = paste0("Minimum n: ",
                        min_nobs, "\n",
                        "Mainstem only: ",
                        FILTER_MAINSTEM.lgl,"\n",
                        "Growth: ",
                        ifelse(grow_logTransf.Lgl, "Log - ",""),
                        ifelse(ppn_grow.Lgl,"Proportional","Linear"), "\n",
                        "Survival-logit transf.: ",
                        surv_logit_transf,"\n",
                        "Minimum DBH:",
                        ifelse(!FILTER_DBH.lgl,
                               "No-filter",
                               paste(DBH_min,"cm")),
                        "\n",
                        "Pearson's correlation coefficient with bootstrapped ",
                        100*round(1 - CI_a, 2),
                        "% CI, (",n_bstps," bootsraps)")) +
  theme_bw()

# End of script



# 
# 
# 
# 
# # Bootstrap preserving site
# 
# # functionfor generating a bootstrap
# bootstrap.f <- function(df){
#   df <- df[sample(1:nrow(df), size = nrow(df), replace = T),]
#   return(df)}
# 
# CI_a <- 0.05 # alpha
# n_bstps <- 100
# 
# bootvals.v <- map(1:n_bstps, ~bootstrap.f(df = data_all.df)) %>%
#   sapply(., FUN = function(X,
#                            a = "ppn_surv",
#                            b = "mean_grow"){
#     return(cor(x = X[[a]], y = X[[b]]))
#   })
# 
# 
# 
# # bootstrap_cor.v <- vector(mode= "numeric", length = n_bstps)
# # for(b in 1:n_bstps){
# #   i.v <- sample(1:nrow(data_all.df),
# #               size = nrow(data_all.df),
# #               replace = T)
# #   bootstrap_cor.v[b] <- cor(data_all.df$ppn_surv[i.v],
# #                             data_all.df$mean_grow[i.v])
# # }
# 
# 
# 
# CI <- quantile(bootvals.v, probs = c(CI_a/2, 1 - CI_a/2) )
# cat("Correlation between growth and mortality = ",
#     round(cor(data_all.df$ppn_surv, data_all.df$mean_grow),2),
#     " with a bootstrapped ", 100*(1-CI_a),
#     "% CI at (", paste(round(CI, 2), collapse = ", "), ")\n",
#     sep = "")
# plot(y = data_all.df$ppn_surv,
#      x = data_all.df$mean_grow,
#      xlab = "mean growth",
#      ylab = "proportion survived")

