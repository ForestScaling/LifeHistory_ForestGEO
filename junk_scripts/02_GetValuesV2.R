# File 02_GetValuesV2.R
# Author: JMR
# Last Update: 12023-06-28
# Description: Takes in Standardized/harmonized raw-stem data and outputs the
# data needed to estimate the individual parameters.
# ### Setup ###################################################################
set.seed(123)

loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data1 <- paste0(loc_Gdr, "/data/ForestGEO/processed/canopy_assign/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/taxa/")
loc_out1 <- paste0(loc_Gdr, "/outputs/")

# Bringing in the packages:
packages <- c("MASS", "magrittr", "lubridate", "knitr",
              "fitdistrplus", "tidyverse", "kableExtra",
              "dplyr")

invisible(lapply(packages, library, character.only = T))

# Setting/finding location of the script, so that the WD is automatically
# within GITHub local copy location.
setwd("C:/Users/juan.m.rodriguez/Documents/ProgramFiles/Github/NASA_ForestScaling/LifeHistory/forestGEO")


# Bring in commonly used functions
source("scripts/00_Functions.R")

# Site vector
site <-  "HVDF"


# Loading in the data
load(paste0(loc_data1, paste0(site, "_AllStemsCL.r")))
AllStems <- get(paste0(site, "_AllStemsCL"))
# load the taxa table
load(paste0(loc_data2, site, "_TaxaTable.r"))
TaxaTable <- get(paste0(site, "_TaxaTable"))

# Table of OTUs to estimate values for:
TaxaTable_2est <- TaxaTable %>%
  filter(DropBegin == 0, DropParamtrz == 0) %>%
  dplyr::select(., "Mnemonic") %>%
  as_tibble() %>%
  mutate(Mnemonic = as.factor(Mnemonic))


verify.df <- AllStems[[1]] %>%
  dplyr::select(., c("Mnemonic", "TreeID", "StemID", "Census"))

# loop through all of the censuses and rbind them together
for(i in 2:length(AllStems)){
  verify.df <- AllStems[[i]] %>%
    dplyr::select(., c("Mnemonic", "TreeID", "StemID", "Census")) %>%
    rbind(., verify.df)}

# Get the number of censuses (BEING NA COUNTS AS A DIFFERENT CENSUS!!!)
n_census <- length(unique(verify.df$Census))
if(NA %in% unique(verify.df$Unique_Census)){
  warning("NA values are present in the census information!")
}

# Summary of the StemIDs across the censuses:
verify.df %<>% group_by(., StemID) %>%
  summarize(Unique_mnem = n_distinct(Mnemonic), # number unique mnemonics
            Unique_TreeID = n_distinct(TreeID), # number unique TreeID
            Unique_Census = n_distinct(Census), # number of unique censuses
            n_record = n()) # number unique records


# Dropping stems that are of OTUs that aren't going to be estimated for
# parameters (e.g. genus only assignments)
dropme <- TaxaTable$Mnemonic[which(as.logical(TaxaTable$DropParamtrz))]
for(i in 1:length(AllStems)){
  AllStems[[i]] <- AllStems[[i]][!(AllStems[[i]]$Mnemonic %in% dropme),]}

# ### Generate values #########################################################
# --- Canopy Level ------------------------------------------------------------
# choice for grouping by canopy level will be put here, currently for HVDF, we
# are manually setting it to drop canopy Level 4 (merge with level 3).
if(site == "HVDF"){
  n_Lvls <- 3
}else{
  warning("Number of canopy layers not picked for this site!")}


AllStems <- lapply(AllStems, function(df){
  df %<>% mutate(CanopyLvl = replace(CanopyLvl, CanopyLvl > n_Lvls, n_Lvls))
  return(df)
})


# --- Stature -----------------------------------------------------------------
stature.df <- AllStems[[1]]
# Question for considering the largest diameter trees, should we drop all
# those with their diameter taken at a different height? Or is this just
# a negligible effect?
stature.df %<>% select(., c("Mnemonic", "TreeID", "StemID",
                            "Census", "DBH","Codes","RametStatus"))

# loop through all of the censuses and rbind them together
for(i in 2:length(AllStems)){
  stature.df <- AllStems[[i]] %>%
    select(., c("Mnemonic", "TreeID", "StemID", "Census",
                "DBH","Codes","RametStatus")) %>%
    rbind(., stature.df)}

# Only keep those that have DBH measures and are alive
stature.df %<>% filter(., RametStatus == "A", !is.na(DBH)) %>%
  mutate_at(., c("TreeID", "StemID"), as.factor) %>%
  group_by(., StemID)

# Get the maximum DBH for the StemID, also checks that their is only 1 taxon
# and TreeID assigned per live-non NA DBH StemID
stature.df %<>% summarise(., maxDH = max(DBH, na.rm = T))

# Re-add the TreeID and taxon (mnemonic) information
stature.df <- AllStems[[1]] %>% select(c("StemID", "TreeID", "Mnemonic")) %>%
  mutate_all(., as.factor) %>%
  right_join(., stature.df, by = "StemID") %>%
  as_tibble(.)

# Estimate the Parameter here:
# For each unique mnemonic find the largest 6 individual live stems
# and average their DBH.
# future options could be using the X-percentile (95? or 90?)
stature.df %<>% group_by(., Mnemonic) %>%
  summarise(., statureDBH = mean(sort(maxDH, decreasing = T)[1:6],
                                 na.rm = T),  number = n())

# Add back in the taxa that couldn't be estimated.
stature.df %<>% right_join(., TaxaTable_2est, by = "Mnemonic")
# --- Condense data -----------------------------------------------------------
# Get binary variables denoting growth and survival of individual of class
# i from interval t to t + 1 for the largest stem in a group
condensed.L <- rep(list(NULL), length(AllStems) - 1)

# Loop aggregating the data into one spot
for(c in 1:(length(AllStems) - 1)){
  # Put in the census information
  condensed.L[[c]] <- AllStems[[c]] %>%
    select(., c("StemID", "TreeID",
                "Mnemonic", "Census",
                "Stem", "CanopyLvl",
                "RametStatus", "GenetStatus",
                "DBH", "Codes",
                "Date")) %>%
    as_tibble(.) %>%
    mutate(Date = ymd(Date)) %>%
    rename(RametStatus_bef = RametStatus,
           GenetStatus_bef = GenetStatus,
           DBH_bef = DBH,
           Codes_bef = Codes,
           Date_bef = Date) # before
  
  # Put the subseqent census information right next to it.
  condensed.L[[c]] <- AllStems[[c + 1]] %>%
    select(., c("StemID", "RametStatus",
                "GenetStatus", "DBH",
                "Codes", "RecruitStatus",  "Date")) %>%
    as_tibble(.) %>%
    mutate(Date = ymd(Date)) %>%
    rename(RametStatus_aft = RametStatus, # after
           GenetStatus_aft = GenetStatus,
           DBH_aft = DBH,
           Codes_aft = Codes,
           Date_aft = Date) %>%
    left_join(condensed.L[[c]], ., by = "StemID")
  
  
  # Add a column denoting the duration that elapsed between the two measurements
  condensed.L[[c]] %<>%
    mutate(Date_dur = int_length(ymd(Date_bef) %--% ymd(Date_aft))) %>%
    mutate(Date_dur = as.duration(Date_dur)/dyears(1))
}

# --- Growth ------------------------------------------------------------------
# Using growth-2
# Data pruning: if >=100 growth values are available, top 1% of values are
# dropped
# regardless, all values of growth <1 (negative growth) are removed.
# All subsequent values then are -0.999, which changes this to a percent
# increase per-year  with an additional constant to force all non-growing stems
# to be positive for the lognormal.


# Calculate growth2
growth.l <- lapply(condensed.L, function(x){
  # Filtering out all of the stems that we do NOT want to take growth
  # measures from, like those with POM changes, dead stems, missing stems,
  # etc.
  
  t_rf <- 3 # The final time across which we want to consider the growth rate
  # note that since the re-measurement intervals are ~5y this centers
  # it along that interval
  
  t_ri <- 2 # The initial time across we want the linear growth rate, note 
  # t_rf - t_ri is 1 year, so this is a linear growth rate per year
  # within the middle of the expected interval
  
  # Currently using the formulation in the BCI paper (change in DBH across 
  # the ~midpoint year)
  # Another option is to treat it as a (DBH_aft/DBH_bef)^(1/Date_dur)
  
  # Who to include
  x %<>% mutate(DBH_use = RametStatus_bef == "A" & # stem living at t_i
                  RametStatus_aft == "A" &  # stem must be living at t_f
                  !CodesLogic(Codes_aft, "A") & # stems not have POM shift
                  (!(Stem == "Secondary") | is.na(Stem)) & # only main stems
                  !is.na(DBH_bef) & # must have a DBH at both t_f and t_i
                  !is.na(DBH_aft))
  
  # Actual estimation
  x %<>% mutate(
    growth1 = case_when(DBH_use ~ DBH_bef*(
      (DBH_aft/DBH_bef)^(t_rf/Date_dur) - (DBH_aft/DBH_bef)^(t_ri/(Date_dur))),
      TRUE ~ NA_real_),
    growth3 = case_when(DBH_use ~ (DBH_aft - DBH_bef)/Date_dur,
                        TRUE ~ NA_real_),
    growth2 = case_when(DBH_use ~ (DBH_aft/DBH_bef)^(1/Date_dur),
                        TRUE ~ NA_real_)) %>%
    select(c("StemID", "TreeID", "Mnemonic", "Census", "Stem", "CanopyLvl",
             "Date_dur" ,  "growth2", "Codes_bef", "Codes_aft",
             "DBH_bef", "DBH_aft"))  
  
  return(x)
  
})

# Go by species and remove the negative growth values, and top 1% (if there
# are >=100 obs)
  n_vals2keep <- 100 # maximum number of observations used, if NA, use all
uprOutlier <-  1 - 0.01 # percentile for removing upper values. If NA
                        # don't remove any right tail values.
                        # It only will remove these values if we have more 
                        # than the minimum n_vals2keep

growthobs.L <- list(NULL)


for(l_i in 1:length(growth.l)){
  growth.df1 <- growth.l[[l_i]]
  
  # Make a table to keep track of the number of observations for each OTU
  # (number lost to same genet, number dead, etc.)
  growthobs.df <- growth.df1 %>%
    group_by(Mnemonic, CanopyLvl) %>%
    summarize(n_stem = n(),
              n_g2 = length(which(!is.na(growth2))),
              n_neggrow = length(which(growth2 < 1))) %>%
    mutate(n_used = NA_integer_,
           MnLvl = paste(Mnemonic, CanopyLvl, sep = "_"))
  growthobs.df %<>% filter(!is.na(CanopyLvl))
  
  # Add in entries for OTU-canopy level combinations not found (fill
  # values with NA
  df1 <- data.frame(Mnemonic = sort(
    rep(unique(TaxaTable_2est$Mnemonic), n_Lvls)),
    CanopyLvl = rep(1:n_Lvls, length(unique(TaxaTable_2est$Mnemonic))),
    n_stem = 0,
    n_g2 = NA,
    n_neggrow = NA,
    n_used = 0)
  df1 %<>% mutate(MnLvl = paste(Mnemonic, CanopyLvl, sep = "_"))
  
  # Add in the unrepresented combinations
  growthobs.df %<>% bind_rows(., df1[!(df1$MnLvl %in% growthobs.df$MnLvl),])
  
  for(Lvl in 1:n_Lvls){
    growth.df2 <- growth.df1 %>% filter(CanopyLvl == Lvl)
  for(OTU in unique(growth.df2$Mnemonic)){
    #if(OTU == "vaccco" & Lvl == 3){stop()}
    #print(paste(OTU, Lvl, sep = "_"))
    OTU.df <- growth.df2 %>% 
      filter(Mnemonic == OTU, # Must be in the OTU of interest
             !is.na(growth2)) # has to have a value for growth # growth can't be negative
    
    # Number of stems with growth2 values
    n_g2 <- nrow(OTU.df)
    
    # filter out so growth can't be negative
    OTU.df %<>% filter(growth2 > 1) # HERE
    
    # Verify that there are observations available
    if(!(nrow(OTU.df) >= 1)){
      growthobs.df %<>%
        mutate(., n_used = case_when(Mnemonic == OTU &
                                       CanopyLvl == Lvl ~ as.numeric(0),
                                     TRUE ~ n_used))
    }

    # Check whether to remove right tail outliers
    if((is.na(n_vals2keep)|nrow(OTU.df) >= n_vals2keep) & !is.na(uprOutlier)){
      OTU.df %<>% arrange(growth2)
      
      # If 
      nrow(OTU.df) - floor((1 - uprOutlier)*n_g2)
      
      OTU.df <- OTU.df[1:max(n_vals2keep,
                             nrow(OTU.df) - floor((1 - uprOutlier)*n_g2),
                             na.rm = T), ]
                            # discards up to 1% of rows that are on the
                            # right tail, will not discard more than to keep

    }
    
    # Check and randomly choose values
    if(nrow(OTU.df) >= n_vals2keep & !is.na(n_vals2keep)){
      
      OTU.df <- OTU.df[sample(1:nrow(OTU.df),
                              size = n_vals2keep,
                              replace = T), ]
      
    }
    
    # add the number of observations 
    growthobs.df <- growthobs.df %>%
      mutate(.,n_used = case_when(Mnemonic == OTU &
                                    CanopyLvl == Lvl ~ as.numeric(
                                      nrow(OTU.df)),
                                  TRUE ~ n_used)) 
    
    growth.df1 %<>% filter(Mnemonic != OTU | CanopyLvl != Lvl) %>%
      bind_rows(OTU.df)
    
  }}
  growth.l[[l_i]] <- growth.df1 %>%
    filter(!is.na(growth2)) %>%
    mutate(growth2 = growth2 - 1 + 0.000) # force to be %increase per year +0.001 #HERE
                                          # a constant so all values are >0
  
  
  growthobs.L[[l_i]] <- growthobs.df
}

# --- Survival ----------------------------------------------------------------
# simple survival values at both the genet and ramet level. Excludes (NA)
# for all missing stems

survive.l <- lapply(condensed.L, function(x){
  # Generate the value
  x %<>% mutate(ram_survived = case_when(RametStatus_bef == "A" &
                                       RametStatus_aft == "A" ~ as.logical(1),
                                     RametStatus_bef == "A" &
                                       RametStatus_aft == "D" ~ as.logical(0),
                                     TRUE ~ NA),
                gen_survived = case_when(GenetStatus_bef == "A" &
                                           GenetStatus_aft == "A" ~
                                           as.logical(1),
                                         GenetStatus_bef == "A" & 
                                           GenetStatus_aft == "D" ~
                                           as.logical(0),
                                         TRUE ~ NA
                                         ))
  
  # cut the unnecessary columns
  x %<>% 
    select(., c("StemID", "TreeID", "Mnemonic",
                     "Census", "Stem", "CanopyLvl", "Date_dur",
                     "ram_survived", "gen_survived"))
  
})


# --- Recruitment -------------------------------------------------------------
recruit.l <- lapply(condensed.L, function(x){
  
  # Which census is this?
  census_name <- unique(x$Census)
  
  # We need to know the "adult" size
  # cutoff of >50% of stature measure
  # other potential options
  # - 50% maximum
  # - median (or some other quantile)
  adult_cutoff <- stature.df %>%
    mutate(cutoff = statureDBH/2) %>%
    select(Mnemonic, cutoff) 
  
  
  # calculate the number of recruits
  n_rect.df <- x %>%
    filter(., RecruitStatus %in% c("GR", "RR")) %>%
    mutate(dummy_GR = as.logical(RecruitStatus == "GR")) %>% 
    group_by(., Mnemonic) %>%
    summarise(., n_GR = n_distinct(TreeID[dummy_GR]), n_recrt = n())
  
  # add in those that had 0 recruitment
  n_rect.df <- TaxaTable_2est %>%
    mutate(n_GR = 0, n_recrt = 0) %>%
    filter(!(Mnemonic %in% n_rect.df$Mnemonic)) %>% rbind(n_rect.df, .)
  
  
  
  # Going through and calculating rates accounting for abundance
  x %<>% filter(RametStatus_bef == "A") # stem was alive at t_i
  
  # calculate recuitment rates
  # recrt_am <- recruits per basal area (stems above DBH adult cutoff) at t1
  # recrt_ae <- recruits per basal area of all conspecifics alive at t1
  # recrt_sm <- recruits per stem (stems above DBH adult cutoff) at t1
  # recrt_se <- recruits per stem of all conspecifics alive at t1
  
  x_all <- x %>%
    group_by(Mnemonic) %>%
    summarise(n_stem_all = n(),
              n_gnt_all = n_distinct(TreeID),
              area_stem_all = sum( pi*(((DBH_bef/100)/2)^2),
                                   na.rm = T))
  
  x_mat <- x %>%
    left_join(., adult_cutoff, by = "Mnemonic") %>%
    filter(., DBH_bef >= cutoff) %>%
    group_by(Mnemonic) %>%
    summarise(n_stem_mat = n(),
              n_gnt_mat = n_distinct(TreeID),
              area_stem_mat = sum(pi*(((DBH_bef/100)/2)^2),
                                  na.rm = T))
  
  denom.df <- full_join(x_mat, x_all, by = "Mnemonic") %>%
    left_join(TaxaTable_2est, .)
  
  n_rect.df <- n_rect.df %>% left_join(., denom.df)
  
  # Making decisions of how recruitment is calculated
  n_rect.df %<>% mutate(recrt_stemsperMBA =
                          n_recrt/area_stem_mat, # nrec per mature BasArea
                        recrt_stemsperMInd =
                          n_recrt/n_stem_mat, # nrec per number mature stems
                        recrt_p1_stemsperMBA =
                          (1 + n_recrt)/area_stem_mat, # nrec+1 per mature BA
                        recrt_p1_stemsperMInd =
                          (1 + n_recrt)/n_stem_mat) # nrec+1 per # matr stems
  
 
  return(n_rect.df)
  
})


# ### Output the data for STAN ################################################
data.L <- list(stature = stature.df,
               recruitment = recruit.l,
               growth_data = growth.l,
               growth_obssum = growthobs.L,
               survival = survive.l,
               n_intervals = length(recruit.l),
               site = site)

assign(x = paste0(site,"_STANdata.L"), data.L)

save(list = paste0(site,"_STANdata.L"),
     file = paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/",
                   site,"_STANdata.r"))