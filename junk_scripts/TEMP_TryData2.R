loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_Cld <- paste0("C:/Users/juan.m.rodriguez/",
                  "OneDrive - University of Maine System/",
                  "GraduateProjects/FIAData/",
                  "TreeLifeHistoryStrats/FIALifeHistory")
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/processed/canopy_assign/")
loc_figOut <- "C:/Users/juan.m.rodriguez/Downloads"


# Attempt to get the current location of the script
setwd("C:/Users/juan.m.rodriguez/Documents/ProgramFiles/Github/NASA_ForestScaling/LifeHistory/forestGEO/scripts")
# setwd(getSrcDirectory(function(){}))

# Load in the packages
for(i in c("magrittr", "tidyverse")){
  library(i, character.only = T)}

# --- Compile a list of all species across all sites --------------------------
# loop through the different data files
j <- 1
Data.L.L <- list()
for(i in grep("STANdata\\.r", list.files(loc_data))){
  file.i <- list.files(loc_data)[i]
  load(paste(loc_data, file.i, sep = "/"))
  Data.L.L[[j]] <- STANdata.L
  names(Data.L.L)[j] <- STANdata.L$site
  
  j <- j+1
}

# create a dataframe to hold all the taxa from all the sites
taxaData_AllSites.df <- data.frame()

# loop throguh the different taxa files
for(i in list.files(paste0(loc_Gdr, "/data/ForestGEO/taxa"))){
  # get the site's 4 letter code
  site <- substr(i, start = 1, stop = 4)
  # Load in the data
  load(paste0(loc_Gdr, "/data/ForestGEO/taxa/",i))
  load(paste0(loc_data2, site, "_AllStemsCL.r"))
  
  stems.df <- AllStems[[1]] %>%
    filter(RametStatus == "A" & !is.na(DBH)) %>%
    mutate(BA = pi*(DBH/2)^2) %>%
    summarise(.by = c(Mnemonic), n_liveStem = n(), Total_BA = sum(BA)) %>%
    arrange(Mnemonic)
  
  # add in column to have site information.
  Spp_table <- Spp_table %>%
    mutate(Site = site, .before = 1)
  
  # add in the number of live stems used for stature
  Spp_table <- stems.df %>%
    left_join(Spp_table, ., by = "Mnemonic") %>%
    mutate(n_liveStem = replace_na(n_liveStem, 0),
           Total_BA = replace_na(Total_BA,0 ),
           p_liveStem = n_liveStem/sum(n_liveStem, na.rm = T),
           p_BA = Total_BA/sum(Total_BA),
           site = site) %>% 
    arrange(Mnemonic)
  
  # Plot rank table
  plt <- Spp_table %>%
    filter(Total_BA > 0 & !is.na(Total_BA)) %>%
    arrange(desc(Total_BA)) %>%
    mutate(rank = 1:nrow(.)) %>%
    ggplot(aes(x = rank, y = p_BA*100, label = Mnemonic)) +
    geom_point() +
      #geom_label_repel(aes(label = Mnemonic))+
    geom_text(hjust = 0, vjust = 0, angle  = 45) +
    theme_bw()+
    scale_y_log10(name= "Living basal area (%)") +
    ggtitle(label = paste(site, "rank Basal Area"))

  ggsave(filename = paste0(loc_figOut,
                           "/",
                           gsub(" ", "_",x =  plt$labels$title),
                           ".jpeg"),
         plot = plt, scale = 2)
  
  plt <- Spp_table %>%
    filter(n_liveStem > 0 & !is.na(n_liveStem)) %>%
    arrange(desc(n_liveStem)) %>%
    mutate(rank = 1:nrow(.)) %>%
    ggplot(aes(x = rank, y = p_liveStem*100, label = Mnemonic)) +
    geom_point() +
    #geom_label_repel(aes(label = Mnemonic))+
    geom_text(hjust = 0, vjust = 0, angle  = 45) +
    theme_bw() +
    scale_y_log10(name= "Living stems (%)") +
    ggtitle(label = paste(site, "rank Stem Count"))
  
  ggsave(filename = paste0(loc_figOut,
                           "/",
                           gsub(" ", "_",x =  plt$labels$title),
                           ".jpeg"),
         plot = plt, scale = 2)
  
  plt <- Spp_table %>%
    filter(Total_BA > 0 & !is.na(Total_BA)) %>%
    ggplot(aes(x = log(p_liveStem/(1-p_liveStem)),
               y = log(p_BA/(1-p_BA)),
               label = Mnemonic)) +
    geom_point() +
    geom_text(hjust = 0, vjust = 0, angle  = 0) +
    geom_abline(slope = 1, intercept = 0) +
    theme_bw() +
    scale_x_continuous(name = "logit(ppn stems)") +
    scale_y_continuous(name = "logit(ppn BA)") +
    ggtitle(label = paste(site, "Live Stem Count v. Basal Area"))
  ggsave(filename = paste0(loc_figOut,
                           "/",
                           gsub(" ", "_",x =  plt$labels$title),
                           ".jpeg"),
         plot = plt, scale = 2)
  
  
  # Bind to the running list of species (NOT FILTERED FOR DUPLICATES)
  taxaData_AllSites.df <- bind_rows(taxaData_AllSites.df, Spp_table)
}

# Discard the mnemonics that weren't included in the analysis
taxaData_AllSites.df <- taxaData_AllSites.df %>% 
  filter(!(as.logical(DropBegin) | as.logical(DropParamtrz)))

# filter so that duplicates don't show up
# table(paste(taxaData_AllSites.df$Genus,
#             taxaData_AllSites.df$Species, sep = "_"))



# taxaData_AllSites.df <- taxaData_AllSites.df %>%
#   summarise(., .by = c(Family, Genus, Species), nSites = n()) %>% 
#   arrange(Genus)


# --- Reading in the TRY data -------------------------------------------------
allTRY <- read.table(paste0(loc_Cld,"/data/procTRY/all_raw_try.txt"),
                     sep = "\t",
                     header = T,
                     quote = "",
                     fill = T,
                     skipNul = T)

# Species/OTU names in the dataset for analysis
Data_sppNames.v <- paste(taxaData_AllSites.df$Genus,
                         taxaData_AllSites.df$Species) %>%
  unique() %>%
  sort

# vector of all species in TRY (regardless of they are in ForestGEO plots)
TRY_sppNames.v <- unique(allTRY$AccSpeciesName) %>% sort

# Traits tied with TraitID
TraitID.df <- allTRY %>%
  filter(!is.na(TraitID) & AccSpeciesName %in% Data_sppNames.v) %>%
  summarise(.by = c(TraitID, TraitName), n_obs = n())

# Datasets tied with DatasetID 
Dataset.df <- allTRY %>%
  filter(!is.na(TraitID) & AccSpeciesName %in% Data_sppNames.v) %>%
  summarise(.by = c(DatasetID, Dataset), n_traits = n_distinct(TraitID))


# Try data must be of sp found in site(s), and have valid trait and dataset
TRY2 <- allTRY %>% 
  filter(!is.na(DatasetID) & AccSpeciesName %in% Data_sppNames.v) %>%
  filter(!is.na(TraitID))

# Join the try data to proportion coun
TRY2 <- TRY2 %>%
  rename("Latin" = AccSpeciesName) %>%
  left_join(., taxaData_AllSites.df,
            by = "Latin",
            relationship = "many-to-many")

n_spp.df <- taxaData_AllSites.df %>%
  summarize(.by = Site, n_Tspp = n_distinct(Mnemonic))

TRY3 <- TRY2 %>%
  summarise(.by = c(Latin, DatasetID, TraitID, Site),
            p_liveStem = unique(p_liveStem)[1],
            p_BA = unique(p_BA)[1],
            Total_BA = unique(Total_BA)[1],
            n_obs = n(),
            n_liveStem = unique(n_liveStem)[1],
            validate = n_distinct(n_liveStem)*n_distinct(p_liveStem))

# Now group by Dataset (collapse species)
TRY3.1 <- TRY3 %>%
  summarise(.by = c(TraitID, Site, DatasetID),
            p_BA = sum(p_BA),
            n_spp = n_distinct(Latin),
            p_LS = sum(p_liveStem)) %>%
  summarize(p_BA_max = max(p_BA),
            p_LS_max = max(p_LS),
            n_spp_max = max(n_spp),
            .by = c(TraitID, Site)) %>%
  left_join(., n_spp.df, by = "Site") %>%
  mutate(., p_spp_max = n_spp_max/n_Tspp) %>%
  select(-c(n_spp_max, n_Tspp))

TRY4 <- TRY3 %>%
  summarise(.by = c(TraitID, Site, Latin),
            n_liveStem = unique(n_liveStem)[1],
            p_liveStem = unique(p_liveStem)[1],
            # p_LS_max = max(p_liveStem),
            p_BA = unique(p_BA)[1],
            Total_BA = unique(Total_BA)[1],
            # p_BA_max = max(p_BA),
            n_obs = sum(n_obs)) %>%
  summarise(.by = c(Site, TraitID),
            n_LS = sum(n_liveStem),
            p_LS = sum(p_liveStem),
            p_BA = sum(p_BA),
            n_spp = n_distinct(Latin)) %>%
  left_join(., n_spp.df, by = "Site") %>%
  mutate(p_spp = n_spp/n_Tspp) %>%
  select(-c(n_spp, n_Tspp)) %>%
  arrange(desc(p_LS))



TRY5 <- left_join(TRY4, TRY3.1, by = c("TraitID", "Site")) %>%
  left_join(., y = TraitID.df, by = c("TraitID")) %>%
  select(c(TraitID, Site, p_LS, p_BA, p_LS_max, p_BA_max, TraitName, p_spp, p_spp_max))
  
write.csv(TRY5, file= paste0(loc_figOut,"/longsummry.csv"))

TRY5 %>%
  select(-c(TraitName)) %>%
  pivot_wider(names_from = Site,
              values_from = c(p_BA, p_LS, p_BA_max, p_LS_max, p_spp, p_spp_max)) %>%
  left_join(., y = TraitID.df[,-3], by = c("TraitID")) %>%
    write.csv(x = ., file = paste0(loc_figOut,"/widesummry.csv"))





# dataframe of the datasets
TRY_datasetIDs.df <- allTRY %>%
  filter(!is.na(DatasetID) & AccSpeciesName %in% Data_sppNames.v)

TRY_datasetIDs.df <- TRY_datasetIDs.df %>%
  summarize(.by = DatasetID, n = n(), n_spp = n_distinct(AccSpeciesName),
            name = unique(Dataset)[1]) %>%
  arrange(DatasetID)

# df of all the traits from the pulled TRYData
TRY_traits.df <- allTRY %>%
  filter(!is.na(TraitID) & AccSpeciesName %in% Data_sppNames.v) %>%
  arrange(TraitID) %>%
  select(c(TraitID, TraitName, AccSpeciesID)) %>%
  summarise(.by = TraitID,
            TraitName = unique(TraitName)[1],
            n_spp = n_distinct(AccSpeciesID)) %>%
  arrange(TraitID)

# Array of species in the datasets, row = OTU, col = TraitID, dim3 = Dataset,
DsCoverage.ary <- array(data = 0, dim = c(length(Data_sppNames.v),
                                          nrow(TRY_traits.df) ,
                                          nrow(TRY_datasetIDs.df)),
                        dimnames = list(Data_sppNames.v,
                                        TRY_traits.df$TraitID,
                                        TRY_datasetIDs.df$DatasetID))

# Functional trait data (ONLY, no metadata, etc.) of species at the sites
# (filtering to reduce computation time in the loop)
usedataTRY.df <- allTRY %>%
  filter(!is.na(TraitID) & AccSpeciesName %in% Data_sppNames.v)

# Go through each dataset and pick out the data
for(i in 1:nrow(TRY_datasetIDs.df)){
  DatasetID.i <- TRY_datasetIDs.df$DatasetID[i] # current dataset in the loop
  subsetTRY.df <- usedataTRY.df %>% filter(DatasetID == DatasetID.i) # filter
  
  # summarize the data so we can g
  dummy1.df <- subsetTRY.df %>%
    summarize(.by = c(AccSpeciesName, TraitID),
              n_vals = n()) %>%
    pivot_wider(names_from = TraitID,
                values_from = n_vals,
                values_fill = 0) %>%
    as.data.frame()
  
  # Match rownames
  rownames(dummy1.df) <- dummy1.df$AccSpeciesName
  
  # Remove the species name column 
  dummy1.df <- dummy1.df %>% select(-c(AccSpeciesName))
  dummy1.df <- as.array(as.matrix(dummy1.df),
                        dim = c(nrow(dummy1.df),
                                ncol(dummy1.df),
                                1))
  
  # fill in the values
  DsCoverage.ary[rownames(dummy1.df),colnames(dummy1.df),i] <- dummy1.df 
  
}









# Coverage matrix
round(100*(apply(DsCoverage.ary[,,] > 0,
                 FUN = sum,
                 MARGIN = c(2,3))/(dim(DsCoverage.ary)[1])), 2) %>%
  View(., "CoverageAry")

# Coverage Trait for each trait 
myfun <-  function(x, k) {
  u <- unique(x)
  sort(u, decreasing = TRUE)[k]
}

# temporary data 
dummy2 <- usedataTRY.df %>%
  summarize(.by = c(TraitID, DatasetID),
            n_spp_wT = n_distinct(AccSpeciesID)) %>%
  summarise(.by = c(TraitID),
            nspp_max = max(n_spp_wT),
            nspp_2nd = myfun(x= n_spp_wT, k = 2))

usedataTRY.df %>%
  summarise(.by = c(TraitName, TraitID),
            n_spp = n_distinct(AccSpeciesID)) %>%
  left_join(y = dummy2, by = "TraitID") %>%
  View(., "Traits")

usedataTRY.df


Traits_keep <- 
  
  
  # Proportion covered by each
  DsCoverage.df <- data.frame(DatasetID = unique(allTRY$DatasetID), cvg = 0)


traitIDs <- paste("TraitID", unique(allTRY$TraitID), sep = "-")
d <- matrix(0, nrow = nrow(DsCoverage.df), ncol = length(traitIDs))
d <- as.data.frame(d)
names(d) <- traitIDs

DsCoverage.df <- cbind(DsCoverage.df, d)


for(i in 1:nrow(DsCoverage.df)){
  DsID <- unique(allTRY$DatasetID)[i]
  allTRY_sbst <- allTRY %>% filter(DatasetID == DsID)
  DsCoverage.df$cvg[i] <- 100*sum(spp_names %in% allTRY_sbst$AccSpeciesName)/length(spp_names)
  
  for(j in unique(allTRY_sbst$TraitID)){
    DsCoverage.df[i, paste0("TraitID-", j)] <-
      unique(allTRY_sbst$AccSpeciesName[
        allTRY_sbst$TraitID == j &
          allTRY_sbst$AccSpeciesName %in% spp_names]) %>%
      unique %>%
      length
  }
  
}
DsCoverage.df <- D

allTRY %>% filter(DatasetID == 439) %>%
  filter( AccSpeciesName %in% spp_names) %>% View




DsCoverage.df