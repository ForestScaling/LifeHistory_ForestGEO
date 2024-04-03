# Load in the packages
for(i in c("MASS", "magrittr", "lubridate", "knitr",
           "fitdistrplus", "tidyverse", "kableExtra")){
  library(i, character.only = T)
}

set.seed(123)

# --- Data import -----------
Gdr_loc <- paste0("G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats/",
                  "data/ForestGEO/processed/example/")
double.df <- read.csv(paste0(Gdr_loc, "wide_HFdata_example.csv"))
# format of the data

# --- Stature -----------------------------------------------------------------
# Calculating Stature:
# I am not sure if this is the best way to do it across the FIA data but here
# was my approach for the ForestGEO data.
stature.df <- double.df %>%
  # pivots longer the data for the before and after census
  dplyr::select(c(StemID, Mnemonic, contains("RametStatus"), contains("DBH"))) %>%
  pivot_longer(cols = -c(StemID, Mnemonic),
               names_to = c(".value","census"),
               names_sep = "_")

# calculation
stature.df <- stature.df %>%
  filter(RametStatus == "A") %>% # stem must be alive
  filter(!is.na(DBH)) %>% # must have measured DBH (no NA)
  # get the largest observation of eah live stem (no double counting
  # a stem over multiple censuses)
  summarize(.by = StemID,
            max_DBH = max(DBH),
            Mnemonic = first(Mnemonic),
            vrfy_Mnmc = n_distinct(Mnemonic)) %T>%
  # verifty that only one Mnemonic was assigned to each stem
  with(., if(any(vrfy_Mnmc > 1)){
    warning("Multiple Mnemonics per stem-ID in stature data!")
    }) %>%
  select(-c(vrfy_Mnmc))

# how wmany stems per species
n.v <- stature.df %>% summarise(.by = Mnemonic, n = n())

stature.df <- stature.df %>% # We don't need the verifcation column
  arrange(desc(max_DBH)) %>% # sort by DBH
  slice_head(., n = 6, by = Mnemonic) %>% # slice off top 6 obs
  summarise(.by = Mnemonic, stature = mean(max_DBH)) %>% # take their average
  left_join(., n.v, by = "Mnemonic")



# --- Growth ------------------------------------------------------------------
growth.df <- double.df

# Verify appropriate data var
growth.df <- growth.df %>%
  # create a column to quickly verify that the data are appropriate for growth
  mutate(DBH_use = # go through the criteria:
           RametStatus_bef == "A" & # stem living at t_i
           RametStatus_aft == "A" &  # stem must be living at t_f
           !as.logical(Codes_aft) & # stems not have POM shift
           (!(Stem == "Secondary") | is.na(Stem)) & # only main stems
           !is.na(DBH_bef) & # must have a DBH at both t_f and t_i
           !is.na(DBH_aft)) %>%
  # filter out those that don't make the cut or don't have a canopy level
  filter(DBH_use & !is.na(CanopyLvl))

# Estimation
growth.df <- growth.df %>%
  mutate(growth = case_when(
    DBH_use ~ (DBH_aft-DBH_bef)/(Date_dur) - 0, # relative is 1, abs is 0
    TRUE ~ NA_real_)) %>%
  select(c("StemID","TreeID","Mnemonic","CanopyLvl","Date_dur","growth"))  

# --- Subsample the data --- ####
p_outliers <- 0.01 # for each species:canopyLvl, drop the p highest values
n_vals2keep <- 100 # number of growth values for subsampling

# dropping values
smy.df <- growth.df %>%
  filter(growth > 0) %>% # relative is 1, abs is 0
  summarise(.by = c(Mnemonic, CanopyLvl), 
            n = n(), # number of stems
            drop = n*(1-p_outliers)) # number to drop

# dataframe of continuous growth calculations
cont_grow.df <- growth.df %>%
  # growth must be positive
  filter(growth > 0) %>% 
  arrange(desc(growth)) %>%
  group_by(Mnemonic, CanopyLvl) %>%
  mutate(id = row_number()) %>% # assign rown number so we can drop largest 
  left_join(., smy.df, by = c("Mnemonic", "CanopyLvl")) %>%
  filter(n < n_vals2keep | id <= drop) %>% # must have at least 100 vals
  # select columns of interest
  select(c(Mnemonic, CanopyLvl, growth, Date_dur, StemID))

# Randomly sample 100 obs per species:canopyLvl combination
cont_grow.df <-  cont_grow.df %>%
  group_by(Mnemonic, CanopyLvl) %>%
  slice_sample(n = n_vals2keep)

# --- Binary growth -----------------------------------------------------------
bin_grow.df <- growth.df %>%
  # binary variable for whether or not the stem grew
  mutate(bin_growth = as.integer(growth > 0)) %>% # absolute > 0, relative > 1
  # only keep the useful columns
  select(c(StemID, Mnemonic, CanopyLvl, bin_growth, Date_dur))

# --- Survival ----------------------------------------------------------------
surv.df <- double.df %>%
  # must be a main stem (or singular stem)
  filter(is.na(Stem) | Stem == "Main") %>%
  # must have a canopy Level assignement # removed those that were excluded
  filter(!is.na(CanopyLvl)) %>%
  # must have been alive at first interval
  filter(RametStatus_bef == "A") %>%
  # must be recorded as dead or alive afterwards (no missing values)
  filter(RametStatus_aft %in% c("A","D")) %>%
  # binary value for survival (1 - survived, 0 - died)
  mutate(survived = as.integer(RametStatus_aft == "A")) %>%
  # filter out other unnecessary data
  select(c(StemID, Mnemonic, CanopyLvl, Date_dur, survived))

# --- recruitment -------------------------------------------------------------
# calculated as both per basal area recruitment and per individuals

# Kambach defined "adult" as 50% maximum, here I use 50% stature height
adult_cutoff <- stature.df %>%
  mutate(cutoff = stature/2) %>%
  select(c(Mnemonic, cutoff))

# calculate mature basal area 
recruit.df <- double.df %>%
  # must be in estimated portions of the plot
  left_join(., adult_cutoff, by = "Mnemonic") %>%
  # calculate mature live Basal area
  mutate(live_MBA = pi*(DBH_bef/200)^2) %>% # basal area, area in m^2
  # dbh > cutoff & is alive
  mutate(live_MBA = live_MBA*as.integer(DBH_bef >= cutoff)) %>%
  # stem must be live in first census
  mutate(live_MBA = live_MBA*as.integer(RametStatus_bef == "A"))

# count the number of mature stems
recruit.df <- recruit.df %>%
  mutate(M_stem = as.integer(DBH_bef >= cutoff &
                               RametStatus_bef == "A"))

# sumamrize values by species
recruit.df <- recruit.df %>%
  # recruitment status of genet or ramet recruit
  mutate(recrt = as.integer(RecruitStatus %in% c("GR", "RR"))) %>%
  # calculate recruitment per species
  summarise(., .by = c(Mnemonic),
            n_recrt = sum(recrt),
            live_MBA = sum(live_MBA, na.rm = T),
            M_stem = sum(M_stem),
            n = n())
  
# calculate the rates
recruit.df <- recruit.df %>%
  mutate(recrt_p1_stemsperMInd = (n_recrt + 1)/M_stem,
         recrt_p1_stemsperMBA = (n_recrt + 1)/live_MBA) %>%
  # to remove weird ones we need to filter out Inf values
  # weird ones due no mature/live individuals previously
  # (can't base rate with an undefined denominator)
  filter(is.finite(recrt_p1_stemsperMInd) &
           is.finite(recrt_p1_stemsperMBA))



# End of calculations
# --- ouput the formatted data (for example scripts) -----------------------
save(list = c("cont_grow.df", "recruit.df", "surv.df",
              "stature.df", "bin_grow.df"),
     file = paste0(Gdr_loc, "example_formatData.R"))

# --- reformat (again) the data for the regular pipeline --------------
stature.df <- stature.df %>%
  rename(statureDBH = stature, number = n) %>%
select(Mnemonic, statureDBH, number)

cont_grow.df <- cont_grow.df %>%
  mutate(growth = exp(growth)) %>%
  rename(growth2 = growth)

bin_grow.df <- bin_grow.df %>%
  mutate(bin_growth = as.logical(bin_growth)) %>%
  rename(grow_Bin = bin_growth)

surv.df <- surv.df %>%
  rename(ram_survived = survived) %>%
  mutate(Stem = NA)

STANdata.L <- list(stature = stature.df,
                   recruitment = list(recruit.df),
                   growth_data = list(cont_grow.df),
                   growth_bin = list(bin_grow.df),
                   growth_obssum = list(NULL),
                   survival = list(surv.df),
                   n_intervals = 1,
                   site = "EXHF")

source("./scripts/AA_LocationManager.R")
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")

# output
save(list = c("STANdata.L"),
     file = paste0(loc_data, "EXHF", "_STANdata.r"))


