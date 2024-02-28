# File: 00_Functions.R
# Creator: JMR
# Date last modified: 12023-01-06
# Description: File of commonly used functions for this project.
# ### Functions ###################################################################################
# --- reader1 -------------------------------------------------------------------------------------
# Reading function to read the ForestGEO data, used for the initial imports
reader1 <- function(f, loc){
  read.table(sep = "\t", header = T, paste0(loc,"/", f))
  }

# --- CodesLogic ----------------------------------------------------------------------------------
CodesLogic <- function(col, c, sep = ";"){
  parttn <- str_split(col, pattern = sep) # partition into the different codes
  return(sapply(parttn, function(x){c %in% x}))
}

# --- List widener --------------------------------------------------------------------------------
ListWiden <- function(List, by = "StemID", kcols){
  # Make a new data frame with a single row for each unique StemID
  tb <- data.frame(y = unique(unlist(lapply(List, function(x){
    return(x[, by])}))))
  
  colnames(tb) <- by
  
  # Separate out each separate time interval
  for(i in 1:length(List)){
  z <- List[[i]]
  z <- z[, which(colnames(z) %in% unique(c(by, kcols)))]
  
  # Rename the columns so that each time can be easily identified
  colnames(z)[-which(colnames(z) == by)] <-
    paste(colnames(z)[-which(colnames(z) == by)], i, sep = "_")
  
  # Join back with the main table
  tb <- left_join(tb, z, by = by)
  }
  return(tb)
}

# ---- any1b4 -------------------------------------------------------------------------------------
# Function to check if any previous df in the list have the value, best for determining if a 
# genet or a ramet is a recruit, or was previously present. Note that it does NOT like "previous"/
# "prior" values, and requires that the entries be in order, and preferably no NA for census
any1b4 <- function(data_list,
                   out_col,
                   check_col = "StemID",
                   t_name = "Census",
                   val_r = c("GR", "RR")){
  
  # Column determination
  col_c <- which(colnames(data_list[[1]]) == check_col)
  col_t <- which(colnames(data_list[[1]]) == t_name)
  col_o <- which(colnames(data_list[[1]]) == out_col)
  
  # Error checking for column name consistency
  lapply(data_list, function(a){
    if(colnames(a)[col_c] != check_col |
       colnames(a)[col_t] != t_name |
       colnames(a)[col_o] != out_col){
      stop("Column names are not in the same order across the list!")
    }
  })
  
  # get a vector of all the censuses represented
  v_censi <- unique(unlist(lapply(data_list, function(y){
    return(y$Census[which(!is.na(y$Census))])})))
  
  out <- lapply(data_list, function(x){
    # Check that their is only one 't' per dataframe
    if(length(unique(x[, col_t])) != 1){stop("More than one census per dataframe!")}
    
    # Get the current census
    curr_census <- as.integer(unique(x[, col_t])[which(!is.na(unique(x[, col_t])))])
    
    # Checks that the censuses are in order and no duplicates across multiple dfs
    if(!identical(v_censi[order(v_censi)], v_censi)){
      stop("Census(es) are not in order, please rearrange!")
    }else if(length(data_list) != length(v_censi)){
      stop("Census(es) are spread out across more than one df!")
    }
    
    # Go through the list
    if(curr_census == min(v_censi)){
      # If this is the first item in the list, skip it since it can't have been recorded previously
      return(x)
    }else{
      # Otherwise, check all those that were before in the list
      prevs <- unique(unlist(lapply(data_list[1:(which(v_censi == curr_census) - 1)], function(y){
        return(y[, col_c])})))
      
      # find all rows that have the identifier in a previous df
      x[which(!(x[, col_c] %in% prevs)), col_o] <- val_r
      
      return(x)
      }
  })
  
  # Function output
  return(out)
}

# --- stemRankAsgn --------------------------------------------------------------------------------
# Function to assign Secondary vs main stems, 
stemRankAsgn <- function(data_list, col_dbh = "DBH", col_out = "Stem",
                         consider_if = "A", col_con = "RametStatus", col_grp = "TreeID",
                         col_ind = "StemID"){
  
  # Check if previous StemAssignments exist
  if(!is.null(bind_rows(data_list)$Stem)){
  check1 <- sum(!is.na(bind_rows(data_list) %>% dplyr::select(all_of(col_out))))
    }else{check1 <- 0}
  # Denote if overriding a previous one
  if(check1 > 0){warning("Overriding previous data inside the", col_out, "column.")}
  
  output <- lapply(data_list, function(x){
    
    # Check that all of those column names are actually represented in the df
    if(length(which(!c(col_dbh, col_out, col_ind, col_con, col_grp) %in% colnames(x))) != 0){
      stop("Columns are not all present in dataframe(s)!")
    }
    
    # Only get those columns that have the value(s) in 'consider_if'
    x_red <- x %>%
      dplyr::filter(get(col_con) %in% consider_if) %>%
      group_by(get(col_grp)) %>%
      summarize(n = n(), 
                biggest = get(col_ind)[order(get(col_dbh), decreasing = T)][1]) %>%
      dplyr::filter(n > 1)
    
    
    
    x %<>% mutate(mystem = case_when(get(col_ind) %in% x_red$biggest ~ "Main",
                                           get(col_grp) %in% x_red$`get(col_grp)` &
                                             get(col_con) %in% consider_if ~ "Secondary",
                                           TRUE ~ NA_character_))
  
    x[, which(colnames(x) == col_out)] <- x[, ncol(x)]
    x <- x[, -ncol(x)]


    return(x)
    
  })
  
  

  
  return(output)
}

# --- QC_Check ------------------------------------------------------------------------------------
# checks for "Zombie ramets, those that come alive 
QC_Check <- function(data_list,
                     do.spp_chck = TRUE,
                     do.chng_HOM = TRUE,
                     do.zomb_ram = TRUE){
  
  # Reorder rows by StemID
  data_list %<>% lapply(., function(x){
    x <- x[order(x$StemID),]
    return(x)
    })
  
  # Verify that rows are in the same order
  for(i in 2:length(data_list)){
    if(!identical(data_list[[i]]$StemID, data_list[[i - 1]]$StemID)){
      stop("StemIDs are not in order!")
    }
  }
  
  
  # Check for more than one census per df & that all rows are arranged the same acro
  lapply(data_list, function(a){
    if(length(unique(a$Census)) != 1){
      stop("More than one census per df in the list, check for NA values!")
      }
    })
  
  
  
  
  # Verify that the list is ordered properly & that their are no duplicate censuses
  # Placeholder with all of the censuses
  t <- unlist(lapply(data_list, function(y){return(unique(y$Census))}))
  # Check that reordering and only calling uniques doesn't do anything.
  if(!identical(t, unique(t)[order(unique(t), decreasing = F)])){
    stop("Censuses are represented across multiple df OR df are not in order!")
    } 
  
  # Check for change in HOM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(do.chng_HOM){
    
    # Verify that rows are in the same order
    for(i in 2:length(data_list)){
      if(!identical(data_list[[i]]$StemID, data_list[[i - 1]]$StemID)){
        stop("StemIDs are not in order!")
      }
    }
    
    chng_HOM <- data.frame(StemID = NULL,
                           Census = NULL,
                           Codes = NULL)
    
    zl <- data_list
    
    for(i in 2:length(zl)){
      
      zl[[i]]$HOM[which(is.na(zl[[i]]$HOM))] <- zl[[i - 1]]$HOM[which(is.na(zl[[i]]$HOM))]
      
      zl_i <- zl[[i]]
      zl_im1 <- zl[[i - 1]]
      
      zl[[i]]$dif <- zl_i$HOM != zl_im1$HOM & !CodesLogic(data_list[[i]]$Codes, "A")
      
      chng <- which(zl[[i]]$dif)
      
      if(length(chng) > 0){
        data_list[[i]]$Codes[chng] <- paste("Edit_HOM", "A",
                                            data_list[[i]]$Codes[chng],
                                            sep = ";")
        
        chng_HOM <- rbind(chng_HOM, dplyr::select(data_list[[i]][chng,], StemID, Census, Codes) )
      }
      
    }
    
    if(nrow(chng_HOM) == 0){
      cat("Passed check for alternate HOM records!\n")
    }else{
      cat(paste(nrow(chng_HOM),
"instances of HOM changing without a recorded 'A' code. Adding codes now, along
  with an 'Edit_HOM' code  to denote that this change in codes was made."))
    }
    
  }else{
    chng_HOM <- NULL
    }
  # End of check for unmarked changes in HOM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Potentially useful parameters, may scrap later
  uniq_stems <- unique(unlist(lapply(data_list, function(z){
    return(z$StemID)})))
  l_l <- length(data_list)
  l_r <- length(uniq_stems)
  
  # Tibble containing all of the data in one spot
  df <- bind_rows(data_list) %>%
    arrange(Census) %>%
    group_by(StemID)

  # Needs to group observations by StemID so that we can look at their histories
  df_r <- df %>%
    group_by(StemID) %>%
    summarise(n_SppIDs = length(unique(ignore_NA(Mnemonic))),
              h_ram = paste(RametStatus, collapse = ""),
              A = paste(collapse = "", as.numeric(CodesLogic(Codes, "A", sep = ";"))),
              v_HOM = paste(collapse = ";", as.character(HOM)))

  # Check for multiple Mnemonic assignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(do.spp_chck){
    leng <- length(which(df_r$n_SppIDs > 1))
    if(leng > 0){
      # Write message if failed check
      cat(paste("A total of", leng,
                "StemID have more than one Mnemonic assigned to them. Please check output!\n"))
      spp_chck <- df_r %>% filter(n_SppIDs > 1) %>% dplyr::select(StemID, n_SppIDs)
    }else{
      # Write message if passed check
      cat("Passed check for multiple Mnemonic assignments to a StemID over time.\n")
      spp_chck <- "PASS"
    }
    # Give NULL if check not conducted
  }else{spp_chck <- NULL}
  # End of check for multiple Mnemonics, begin check for zombie ramets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(do.zomb_ram){
    
    
    df_r2 <- df_r %>% dplyr::select(StemID, h_ram)
    # This is a pain in the rear:
    
    
    ft1 <- "(?<=A)[M](?=M*A)" # all missing values with flanking 'A' (fill with A)
    ft2 <- "[D](?=M*A)"       # 'D' values with subsequent 'A' with 'M' in between (fill with A)
    ft3 <- "[D](?=D*A)"       # 'D' values with later 'A' with other 'D' in between (fill with A)
    ft4 <- "(?<=D)[M](?=M*D)" # 'M' with flanking 'D' values, (fill with D)
    ft5 <- "(?<=P{0,1000})[M](?=A)"
    
    j <- "j"   # Check if we are still making edits (pre)
    j2 <- "j2" # Check if we are still making edits (post)
    
    # While pre =/= post
    while(!identical(j, j2)){
      # new pre
      j <- df_r2$h_ram
      # Make edits
      df_r2$h_ram %<>% str_replace(., pattern = ft1,"A") %>% 
        str_replace(., pattern = ft2, "A") %>% 
        str_replace(., pattern = ft3, "A") %>% 
        str_replace(., pattern = ft4, "D") %>%
        str_replace(., pattern = ft5, "A")
      # new post
      j2 <- df_r2$h_ram 
    }


    
    
    
    
    
    
    
    zomb_ram <- df_r %>%
      dplyr::select(StemID, h_ram) %>%
      rename(h_ram_i = h_ram) %>%
      left_join(., rename(df_r2, h_ram_f = h_ram), by = "StemID")
    
    chng <- 0
    
    # Now we need to re-add them back into the data_list
    for(i in 1:length(data_list)){
      
      if(identical(data_list[[i]]$StemID, df_r$StemID)){
        print("PASS")
      }else{print("FAIL")
          }
      
      # Count the number of different ones
      a <- str_sub(zomb_ram$h_ram_f,i,i) != data_list[[i]]$RametStatus
      # The number of Deads changed to alive.
      chng <- chng + length(which(a & (data_list[[i]]$RametStatus == "D")))
      
      data_list[[i]] %<>% mutate(Codes = case_when((RametStatus == "D") &
          (RametStatus != str_sub(zomb_ram$h_ram_f,i,i)) ~ paste("Edit_zram", Codes, sep = ";"),
        TRUE ~ Codes),
        RametStatus = str_sub(zomb_ram$h_ram_f,i,i))
      }
    
    cat(paste("A total of", chng, "RametStatus records were changed since they represent",
              "'zombies'. Added an 'Edit_zram' status to reflect the update.\n"))
    
    cat(paste0("Filling in 'M' statuses  with 'D' or 'A' if flanked by values. E.g. if",
               " a ramet's history was 'D' -> 'M' -> 'D', the middle time is assumed to also be",
               " dead. In the case of 'A' -> 'M' -> 'D', the 'M' ramet status is maintained",
               " since when death ocurred cannot be inferred.\n"))
    
  }else{zomb_ram <- NULL}
  # end of checking for zombie ramets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Make sure the updated codes don't do anything like ';A' or 'A;B;'
  data_list %<>% lapply(., function(x){
    x$Codes <- str_replace(x$Codes, pattern = "^;", replacement = "") %>%
      str_replace(., pattern = ";$", replacement = "")
    
    return(x)
  })
  
  # Updating genet Statuses
  data_list %<>% lapply(., function(z){
    return(genetAssign(z, evid_rcrt = T))
    })
  

  
  
  # result
  result <- list(QCdata = data_list,
                 chng_HOM = chng_HOM,
                 spp_chck = spp_chck,
                 zomb_ram = zomb_ram)
  
  return(result)
    
  }

# --- MP_adder ------------------------------------------------------------------------------------

MP_adder <- function(data_list){
  
  
  
  # Check for more than one census per df
  lapply(data_list, function(x){
    if(length(unique(x$Census)) != 1){
      stop("More than one census per df in the list, check for NA values!")
      return()
    }else{
      return()}
  })
  
  # Verify that the list is ordered properly & that their are no duplicate censuses
  # Placeholder with all of the censuses
  t <- unlist(lapply(data_list, function(y){return(unique(y$Census))}))
  # Check that reordering and only calling uniques doesn't do anything.
  if(!identical(t, unique(t)[order(unique(t), decreasing = F)])){
    stop("Censuses are represented across multiple df OR df are not in order")
  } 
  
  uniq_stems <- unique(unlist(lapply(data_list, function(z){return(z$StemID)})))
  
  # Add in TreeID, check that it has one TreeID over time.
  tr <- tibble(StemID = uniq_stems)   # Table with all unique stems
  tr <- data_list[[1]] %>%# Add in the TreeID recorded from the earliest census
    dplyr::as_tibble() %>% 
    dplyr::select(c(StemID, TreeID, Mnemonic)) %>%
    {left_join(tr, ., by = "StemID")}
  
  
  # ensure that their are no StemIDs that are only missing one.
  if(nrow(tr) != length(which((is.na(tr$TreeID) & is.na(tr$Mnemonic)) |
                              (!is.na(tr$TreeID) & !is.na(tr$Mnemonic))))){
    warning("Verify that you do not have any rows entires missing either TreeID or Mnemonic")
    
  }
  
  # Fill in all known TreeID, piece meal as they show up.
  for(i in 2:(length(data_list))){
    
    # Already has a TreeID recorded
    nottoadd <- tr %>% filter(!(is.na(TreeID) & is.na(Mnemonic)) ) 
    # No TreeID recorded yet
    toadd <- tr %>% filter(is.na(TreeID) & is.na(Mnemonic)) %>% dplyr::select(., StemID)
    
    # Add the known TreeID at census and filling in missing NA, and reconnect with the knowns
    tr <- as_tibble(data_list[[i]]) %>%
      dplyr::select(c(TreeID, StemID, Mnemonic)) %>%
      left_join(toadd, ., by = "StemID") %>%
      bind_rows(nottoadd, .) %>%
      arrange(., "StemID")
  }
  
  
  # Check that TreeID assignments are consistent over time. If StemID 1 at t=1 is assigned
  # TreeID 1, make sure that (if present) its TreeID is still 1 in future censuses. 
  lapply(data_list, function(x){
    # Only look at the StemID and TreeID columns, and rearrange so they are ordered by StemID
    x2 <- x %>%
      dplyr::select(c(StemID, TreeID, Mnemonic)) %>%
      arrange(., StemID)
    
    # Subset the "entire" set of StemID by those present at current census
    x3 <- tr[which(tr$StemID %in% x2$StemID),] %>% arrange(StemID)
    
    # Logical, provides the warning
    if(sum(x3$TreeID != x2$TreeID, na.rm = T) != 0 |
       sum(x3$Mnemonic != x2$Mnemonic, na.rm = T) != 0 ){
      warning("TreeID and/or taxon are NOT conserved, please address")
      
      # Start pulling out the problematic ones
      #x4 <- x2 %>%
      #rename(TreeID_i = TreeID) %>%
      #left_join(., x3, by = "StemID") %>%
      #filter(TreeID != TreeID_i & !is.na(TreeID))
    }
  })
  
  
  data_list %<>% lapply(., function(x){
    
    # Create a vector of the missing StemID from this particular df
    msng <- uniq_stems[which(!(uniq_stems %in% x$StemID))]
    
    x %<>% arrange(StemID)
    
    # If there are none missing go to the next df
    if(length(msng) == 0){
      x %<>% mutate(Recorded = 1)
      
      # Otherwise 
    }else{
      x %<>% mutate(Recorded = 1)
      x <- tibble(StemID = msng,
                  Census = unique(x$Census),
                  RametStatus = "_",
                  Recorded = 0) %>% 
        bind_rows(x, .)
      
    }
    x %<>% arrange(StemID)
    return(x)
  })
  
  tb <- bind_rows(data_list) %>%
    group_by(StemID) %>%
    summarize(h = paste0(RametStatus, collapse = ""))
  
  tb %<>% mutate(h = str_replace_all(string = h, "P", "_")) %>%
    mutate(h = str_replace_all(string = h, "^_+", "")) %>% 
    mutate(h = str_pad(string = h, pad = "P", side = "left",
                       width = max(str_length(tb$h)))) %>%
    mutate(h = str_replace_all(h, "_", "M"))
  
  
  # And now we add them back in 
  for(i in 1:length(data_list)){
    tb_i <- tb %>% mutate(h_i = str_sub(h, i, i )) %>% dplyr::select(StemID, h_i)
    
    tb_i <- tb_i[match(tb_i$StemID, data_list[[i]]$StemID),]
    
    data_list[[i]]$RametStatus <- tb_i$h_i
  }
  
  # add the TreeID to the missing values.
  data_list %<>% lapply(., function(x){
    
    # skip the df iff no missing TreeID values
    if(length(which(is.na(x$TreeID) & is.na(x$Mnemonic))) == 0){
      return(x)}
    
    # Fill them in if there are missing ones
    x %<>% mutate(TreeID = tr$TreeID[match(.$StemID, tr$StemID)],
                  Mnemonic = tr$Mnemonic[match(.$StemID, tr$StemID)])
    
    return(x)
    
  })
  
  
  return(data_list)
}




# --- ignore_NA -----------------------------------------------------------------------------------
ignore_NA <- function(v){
  return(v[which(!is.na(v))])
}
# --- genetAssign ---------------------------------------------------------------------------------
genetAssign <- function(data, evid_rcrt = TRUE, resprt_code = "R"){
  print(paste("1:", "GenetStatus" %in% colnames(data)))
  cx <- FALSE
  
  if(!is.null(data$GenetStatus)){
    
    # Check if there already exists a column with the name 'GenetStatus'
    warning("A column of GenetStatus already exists, overriding.")
    # hold onto the names and order of the columns
    cx <- TRUE
    cn <- colnames(data)
    data %<>% dplyr::select(-c(GenetStatus))
    }

  data %<>% 
    mutate(r_resprt = CodesLogic(.$Codes,
                                 resprt_code, sep = ";")) %>%      # column to check if resprt
    group_by(TreeID) %>%
    summarise(uniqs = paste(unique(RametStatus), collapse = ""),   # Check if resprt + list r_sts
              g_resprt = TRUE %in% r_resprt) %>%
    mutate(GenetStatus = case_when(str_detect(uniqs, "A") ~ "A",   # Alive if single live stem
                                   g_resprt ~ "A",                 # Alive if incl resprt
                                   uniqs == "P" ~ "P",             # prior if only priors
                                   str_detect(uniqs, "P") &        # If at least one future recruit
                                     evid_rcrt ~ "A",              # if all 'P' then genet is 'P'
                                   str_detect(uniqs, "D") ~ "D",   # Dead, if includes dead or miss
                                   str_detect(uniqs, "M") ~ "M",   # Missing if only miss or priors
                                                
                                   TRUE ~ "FAILURE")) %>%          # Check that all were assigned
    dplyr::select(TreeID, GenetStatus) %>% left_join(data, ., by = "TreeID")  # merge with original
  print(paste("2:", "GenetStatus" %in% colnames(data)))
  # If there was already a 'GenetStatus' column before the function, put back the column names
  if(cx){
    data %<>% dplyr::select(all_of(cn))
    }
  print(paste("3:", "GenetStatus" %in% colnames(data)))
  return(data)
}
