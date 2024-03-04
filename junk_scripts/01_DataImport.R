# Script 01, Importing the data,

# Set locations so that the results are replicable
# Setting the GDrive as the location for data storage
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"

# Setting/finding location of the script, so that the WD is automatically within GITHub local copy location.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

# Bringing in the packages:
library("tidyverse")

# Get the commonly used functions
source("scripts/00_Functions.R")

# Site Name:
site <- "HVDF"

# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/",site))[str_which(list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/",site)), pattern = "_AllStem_")]


# Create a file for storing all of the Stem files
All_Stems_raw <- lapply(stem_files, FUN = reader1, loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))

# Read in the species table
Spp_table <- read.delim(paste0(loc_Gdr, "/data/ForestGEO/raw/", site,"/",site,"_TaxaGenusSp_Report.txt"), sep =  "\t", header = T)

# To bring in manually for(i in 1:length(All_Stems_raw)){x = assign(paste0(site,as.character(i)), value = All_Stems_raw[[i]])}

# Creates a vector containing all of the unique StemID values across the entirety of the plot's history.
uniq_StemID <- unique(unlist(lapply(All_Stems_raw, function(x){x$StemID})))

# Checking for consistency
# fill in the other rows in here so we can easily compare across different time intervals
All_Stems <- lapply(All_Stems_raw, function(x){return(data.frame(x,Record = 1))})

# Adding in dummy rows for Not represented 
All_Stems <- lapply(All_Stems, function(x){
  # create a dummy df with NA values for those Stems not found in that particular census
  df = data.frame(matrix(NA,
                         nrow = length(uniq_StemID[which(!(uniq_StemID %in% x$StemID))]),
                         ncol = ncol(x) - 2))
  
  # change the column names so that they match in this "NA" dataset
  colnames(df) = colnames(x)[-which(colnames(x) %in% c("StemID", "Record"))]
  
  # add the extra rows to the dataframe such that all census dataframes are the same size, 
  return(rbind(x, data.frame(StemID = uniq_StemID[which(!(uniq_StemID %in% x$StemID))],
                             Record = 0,
                             df)))})

# Reorder the rows so that they all are arranged the same
All_Stems <- lapply(All_Stems, function(x){return(x[order(x$StemID, decreasing = F),])})

# At this point all of the rows now correspond 1:1 with the other censuses. If a stem has no record for the particular census Record is set to 0,
# otherwise this value is set to 1

# Create some more standardized variables so that we can conduct analysis
All_Stems <- lapply(All_Stems, function(x){return(mutate(x, Alive = na_if(y = 9, case_when(is.na(x$Record) ~ 9,
                                                                       Status %in% c("dead","stem dead") ~ 0,
                                                                       Status %in% c("alive","broken below","alive-not measured") ~ 1,
                                                                       Status == "missing" ~ 9))))})


# check for unique codes
unique(unlist(lapply(All_Stems, function(x){str_split(x$Codes, pattern = ";")})))

# Table with distribution of species
res1 <- All_Stems[[2]] %>% group_by(Latin) %>% summarize(n=n(),percent = 100*n()/nrow(All_Stems[[2]]))
View(res1[order(res1$n,decreasing = T),])
print(nrow(All_Stems[[1]]))


# Weeding through the data 
# How many stems that changed their point of measurement
sapply(All_Stems,)












load(paste0(loc_Gdr,"/data/ForestGEO/raw/HVDF/harvardforest.stem2.rdata"))






list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/HVDF"))













x <- as_tibble(All_Stems[[1]])

x2<-mutate(x, Alive = case_when(Record == 0 ~ NA, TRUE ~ 2))
           






case_when(x$Record == 0 ~ 1, x$Record != 0 ~ 2)







for(i in 2:length(All_Stems)){
  print(which((All_Stems[[i-1]]$Record == All_Stems[[i]]$Record)&(All_Stems[[i-1]]$DBH != All_Stems[[i]]$DBH)))
  
  
  
}


i <- 2

df <- data.frame(StemID = All_Stems[[1]]$StemID,
                 hist = paste(All_Stems[[1]]$Status, All_Stems[[2]]$Status, All_Stems[[3]]$Status, sep = "_"),
                 TreeID = All_Stems[[1]]$TreeID,
                 Latin = All_Stems[[1]]$Latin)

View(df)


# What is "stem dead" vs "dead"
lapply(All_Stems_raw,function(x){str_detect(x$Status, "dead")})


# df of unique stem IDs
All_StemID <- data.frame(StemID = unique(unlist(lapply(All_Stems_raw, function(x){x$StemID}))))

for(i in 1:length(All_Stems_raw)){
  
  All_StemID <- left_join(All_StemID, All_Stems_raw[[i]], by = "StemID")
  
}





# Start with the SERC data:
HVDF_AllStem_01 <- read.table(paste0(loc_Gdr,"/data/ForestGEO/raw/HVDF/HVDF_AllStem_01.txt"), sep = "\t", head = T)
HVDF_AllStem_02 <- read.table(paste0(loc_Gdr,"/data/ForestGEO/raw/HVDF/HVDF_AllStem_02.txt"), sep = "\t", head = T)
SCBI_AllStem_02 <- read.table(paste0(loc_Gdr,"/data/ForestGEO/raw/SCBI/SCBI_AllStem_03.txt"), sep = "\t", head = T)

SERC_AllStem_01 <- read.table(paste0(loc_Gdr,"/data/ForestGEO/raw/SERC/SERC_AllStem_02.txt"), sep = "\t", head = T)

SERC_AllStem_01 <- read.table(paste0(loc_Gdr,"/data/ForestGEO/raw/SERC/SERC_AllStem_02.txt"), sep = "\t", head = T)
# List with all of the Stem files of interest
l <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/SCBI"))[str_which(list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/SCBI")), pattern = "_AllStem_")]


SCBI_AllStem_Allraw <- lapply(l, FUN = dumf1, loc = paste0(loc_Gdr, "/data/ForestGEO/raw/SCBI"))



dumf1 <- function(f, loc){
  read.table(sep = "\t", header = T, paste0(loc,"/", f))
  
}

dumf1(loc = paste0(loc_Gdr, "/data/ForestGEO/raw/SCBI"), f ="SCBI_AllStem_02.txt" )


HRVF_AllStem_01



lapply(X = list.files())

