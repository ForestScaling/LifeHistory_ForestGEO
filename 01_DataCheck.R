# Check for consistency


# Read Packages
packages <- c("tidyverse", "magrittr", "rstudioapi")
for(p in packages){
  library(p, character.only = T)
  }

loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
site <- "SCBI"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
source("scripts/00_Functions.R")

# read in the species table.
Spp_table <- read.delim(paste0(loc_Gdr,
                               "/data/ForestGEO/raw/",
                               site,"/",
                               site,"_TaxaGenusSp_Report.txt"),
                        sep =  "\t", header = T)




# Find all the Stem files from this site
stem_files <- list.files(paste0(loc_Gdr, "/data/ForestGEO/raw/",
                                site))[
                                  str_which(list.files(paste0(loc_Gdr,"/data/ForestGEO/raw/",
                                                              site)), pattern = "_AllMeas_")]

# Create a file for storing all of the Stem files
Raw_Stems <- lapply(stem_files,
                    FUN = reader1,
                    loc = paste0(loc_Gdr, "/data/ForestGEO/raw/", site))


# get rid of extraneous columns
Raw_Stems <- lapply(Raw_Stems, FUN = function(X){
  
  # recode Status
  X <- X %>% mutate(Status = case_when(Status == "alive" ~ "A",
                                       Status == "missing" ~ "M",
                                       Status %in% c("dead",
                                                "broken below",
                                                "stem dead") ~ "D",
                                       TRUE ~ NA_character_))
  
  # demarcate those that had a change in HOM
  X <- X %>%
    mutate(Chng_HOM = as.integer(CodesLogic(Codes, c = "A", sep = ";")))
  
  X <- X %>% mutate(DBH = ifelse(DBH == 0, NA, DBH))
  
  # Filter the columns
  X <- X %>%
    as_tibble(.) %>%
    dplyr::select(., all_of(c("StemID", "TreeID", "Mnemonic", "Chng_HOM",
                       "PX", "PY", "DBH", "Status", "Census")))
  
  

})


# ### Verify Format ###########################################################
# --- Column Names ------------------------------------------------------------
# Verify column names:
cat("Checking that column names are formatted correctly...")
vfy <- lapply(Raw_Stems, function(x){
  # vector of column names
  cols <- c("StemID", "TreeID", "Mnemonic",
            "PX", "PY", "Chng_HOM",
            "DBH", "Status", "Census")
  # These should match
  out <- identical(sort(cols), sort(colnames(x)))
  return(out)
  })

# Error message check
if(!all(as.logical(vfy))){
  # fail 
  cat("\n")
  stop("Please verify all data-frames have the proper column names.")
}else{
  # pass
  cat("done.\n")
  }

# --- Status ------------------------------------------------------------------
# Verify proper values for Status
cat("Checking that 'Status' column is formatted correctly...")
vfy <- lapply(Raw_Stems, function(x){
  # All stem status should be A, M, D or P
  out <- all(x$Status %in% c("A","M","D","P"))
  return(out)
  })

# Print error message if fails
if(!all(as.logical(vfy))){
  cat("\n")
  stop(paste("Please verify that status column values are restricted to ",
             "'M', 'P', 'A', and 'D'."))
}else{
  # pass
  cat("done.\n")
  }

# --- Mnemonics ---------------------------------------------------------------
# Check that mnemonics are all present in the taxonomic table
cat(paste("Checking that 'Mnemonic' column only contains values",
          "found in the species table..."))

# unique Mnemonics 
mnemonic.v <- Spp_table$Mnemonic

vfy <- lapply(Raw_Stems, function(x){
  # All mnemonics should be in the spp table
  out <- all(x$Mnemonic %in% mnemonic.v)
  return(out)
})

# Print error message if fails
if(!all(as.logical(vfy))){
  cat("\n")
  stop(paste("Please verify that all values are present", 
             "in the species table."))
  
}else{
  # Pass
  cat("done.\n")
  }

# ### Begin adding the other vals #############################################

#  add in "records" for all stems (put P if before recruitment, put M is 
# missing and after recruitment)
Raw_Stems <- MP_adder(Raw_Stems) %>%
  lapply(., FUN = function(x){
    return(dplyr::select(as_tibble(x), -RametStatus))
  })

# Collapse the status over time into a df where each ro
chckSurv.df <- Raw_Stems %>%
  lapply(X = ., FUN = function(x){
    dplyr::select(arrange(x, StemID), Status) %>%
      pull(Status) %>%
      return(.)
    }) %>%
    do.call(paste0, .)

# criteria for correcting errors
ft1 <- "(?<=A)[M](?=M*A)" # all missing vals with flank 'A' (fill with A)
ft2 <- "[D](?=M*A)"       # 'D' vals with later 'A' with 'M' and M between +A
ft3 <- "[D](?=D*A)"       # 'D' vals with later 'A' with 'D' in between +A
ft4 <- "(?<=D)[M](?=M*D)" # 'M' with flanking 'D' values, (fill with D)
ft5 <- "(?<=P{0,1000})[M](?=A)" # M val after P val, make A 

# Pre and post checker
j <- "j"   # Check if we are still making edits (pre)
j2 <- "j2" # Check if we are still making edits (post)

# While pre =/= post
while(!identical(j, j2)){
  # new pre
  j <- chckSurv.df
  # Make edits
  chckSurv.df %<>% str_replace(., pattern = ft1,"A") %>% 
    str_replace(., pattern = ft2, "A") %>% 
    str_replace(., pattern = ft3, "A") %>% 
    str_replace(., pattern = ft4, "D") %>%
    str_replace(., pattern = ft5, "A")
  # new post
  j2 <- chckSurv.df
}

# Replace all the values, if it differs put a correction note.
for(i in 1:length(Raw_Stems)){
  # New values from corrections
  vals <- str_sub(chckSurv.df, start = i, end = i)
  # Edit vals
  Raw_Stems[[i]] <- Raw_Stems[[i]] %>%
    mutate(Edit = as.integer(Status != vals), # Add column to see if edited
           Status = vals) # update the values
}

# --- Add in whether the stem is the primary or secondary stem ----------------
Raw_Stems <- lapply(Raw_Stems, FUN = function(X){
  # Add if main or secondary stem
  X <- X[sample(1:nrow(X), # Randomize order so max DBH ties pick on at rndm
              size = nrow(X),
              replace = F),] %>%
    filter(Status == "A" & !is.na(DBH)) %>% # mains must be alive and non-NA
    group_by(TreeID) %>%
    summarize(Stem = StemID[which.max(DBH)]) %>% # get the max dbh stemID
    left_join(X, ., by = "TreeID") %>%
    # check against stemID if match then it is Main, else 2nd, ignore dead
    mutate(Stem = case_when(is.na(Stem) | Status != "A"  ~ NA_character_,
                            Stem == StemID ~ "Main",
                            Stem != StemID ~ "Secondary",
                            TRUE ~ "CHECK")) 
})

# --- Add commonality in coordinates (also check for changes) -----------------
# get the unique positions for X and Y, it should fail if there are 
# multiple values for a single stem (not including NA)
PX_uniq.v <- Raw_Stems %>%
  # Extract the PX column
  sapply(., FUN = function(X){
  return(dplyr::pull(X, PX))}) %>% 
  # take the unique value (ignore NA)
  apply(MARGIN = 1, FUN = function(Y){
    if(all(is.na(Y))){
      out <- NA
    }else{
      out <- Y[which(!is.na(Y))]
    }
    return(unique(out))
  })

# repeat for Y
PY_uniq.v <- Raw_Stems %>%
  # Extract the PY column
  sapply(., FUN = function(X){
    return(dplyr::pull(X, PY))}) %>% 
  # take the unique value (ignore NA)
  apply(MARGIN = 1, FUN = function(Y){
    if(all(is.na(Y))){
      out <- NA
    }else{
      out <- Y[which(!is.na(Y))]
    }
    return(unique(out))
  }) 

# Verify that only one PX and PY value per stem (exclude NA)
if(!(is.atomic(PY_uniq.v) && is.atomic(PX_uniq.v))){
  stop(paste("Multiple coordinate (PX and/or PY, for a single",
             "stem over the censuses (ignoring NA)! Please fix."))
} 

# Now fill in the values
Raw_Stems <- lapply(Raw_Stems, FUN = function(X){
  X %>% mutate(PX = PX_uniq.v, PY = PY_uniq.v) %>%
    return(.)
  
})



  # apply(MARGIN = 1, FUN = function(X){return(ifelse(,NA, unique(X)))})  %>% 
  # lapply(., FUN = function(X){ifelse(identical(NA, X)length(X))}) %>%
  # unlist %>%
  # {which( . == 2)}

# Ready to output
