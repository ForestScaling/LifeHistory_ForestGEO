# Merge in Trait data

# File Locations
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"

# Location of scripts
# get directory of the script itself (WARNING MAY NOT WORK FOR BATCH JOBS)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
print(getwd())

# Load Packages
for(i in c("magrittr", "rstudioapi", "tidyverse")){
  library(i, character.only = T)
}

# read in the data
NEON.df <- read.csv(paste0(loc_Gdr, "/originalData/FunctionalTraits/", 
                          "NEON_recoveredtraits.csv"), row.names = 1)
FIA.df <- read.csv(paste0(loc_Gdr, "/originalData/FunctionalTraits/",
                         "fiaconf_intervals.csv"), row.names = 1)
# --- Processing the NEON data ------------------------------------------------
# removing the duplicate columns
NEON.df <- NEON.df %>%
  select(., contains("sp.")) %>% # duplicate columns found
  rename(Latin = sp.rowname) %>% 
  mutate(Latin = gsub(pattern = "_", # fixing latin names
                      replacement = " ",
                      Latin),
         Source = "NEON")

# outputing the NEON data
write.csv(NEON.df,
          file = paste0(loc_Gdr,
                        "/data/FunctionalTraits/",
                        "NEON_FoliarTraits.csv"),
          row.names = F)

# --- FIA trait data ----------------------------------------------------------
colnames(FIA.df) <- gsub(pattern = "_upper",
                         replacement = "-U",
                         colnames(FIA.df)) %>%
  gsub(pattern = "_lower",
       replacement = "-L",
       .)

# Retrieve the functional traits
trts.v <- FIA.df %>%
  select(-c(Scientific.Name)) %>%
  colnames(.) %>%
  gsub("-(.)", replacement = "", . ) %>%
  unique(.)

