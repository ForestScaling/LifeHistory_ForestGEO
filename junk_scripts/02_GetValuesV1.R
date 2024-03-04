# File 02_GetValues.R
# Author: JMR
# Last Update: 12023-01-19
# Description: Takes in Standardized/harmonized raw-stem data and outputs the
# data needed to estimate the individual parameters.
# ### Setup ###################################################################

loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data1 <- paste0(loc_Gdr, "/data/ForestGEO/processed/canopy_assign/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/taxa/")
loc_out1 <- paste0(loc_Gdr, "/outputs/")

#' Be sure to undue this hard-setting of working directory afterwards, I am
#' not sure where there issues compiling the report when I don't do this

# Bringing in the packages:
packages <- c("tidyverse", "magrittr", "lubridate", "knitr")
invisible(lapply(packages, library, character.only = T))

# Setting/finding location of the script, so that the WD is automatically
# within GITHub local copy location.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")


# Bring in commonly used functions
source("scripts/00_Functions.R")



# Site vector
site.v <- c("HVDF","SCBI","SERC")

# Default plotting dimension for diagnostic plots
plot_w <- 7*2.54
plot_h <- 7*2.54
plot_u <- "cm"



#+ echo = FALSE

# ### Looping through the different sites. ####################################
for(site in site.v){
  # --- Setup and verify formatting -------------------------------------------
  # load the stem data
  load(paste0(loc_data1, paste0(site, "_AllStemsCL.r")))
  AllStems <- get(paste0(site, "_AllStemsCL"))
  # load the taxa table
  load(paste0(loc_data2, site, "_TaxaTable.r"))
  TaxaTable <- get(paste0(site, "_TaxaTable"))
  
  TaxaTable_2est <- TaxaTable %>%
    filter(DropBegin == 0, DropParamtrz == 0) %>%
    select(c("Mnemonic")) %>%
    as_tibble() %>%
    mutate(Mnemonic = as.factor(Mnemonic))
    
  
  # Start off with verifying the following for each stem
  # 1. it has a record across all censuses, even if just to note that it is
  # missing in a record (status M for genet and/or ramet)
  # 2. There are no changes in identification of the nmnemonic over time (only)
  # one unique value per StemID
  # 3. Its TreeID does not change over time, a stemID can't be assigned to
  # different TreeIDs across different censuses.
  # 4. verify that all StemID are unique within a single census
  verify.df <- AllStems[[1]] %>%
    select(., c("Mnemonic", "TreeID", "StemID", "Census"))
  
  # loop through all of the censuses and rbind them together
  for(i in 2:length(AllStems)){
    verify.df <- AllStems[[i]] %>%
      select(., c("Mnemonic", "TreeID", "StemID", "Census")) %>%
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
  
  # Generate warnings
  if(!identical(unique(verify.df$Unique_mnem), as.integer(1))){
    warning(paste0("For site ", site,
                   " check that mnemonics are constant and without NAs!"))}
  if(!identical(unique(verify.df$Unique_TreeID), as.integer(1))){
    warning(paste0("For site ", site,
                   " check that all TreeID are unique, and do not use NAs!"))}
  if(!identical(n_census,
                unique(verify.df$Unique_Census),
                unique(verify.df$n_record))){
    warning(paste0("For site ", site,
                   " check that StemIDs have exactly one record per census!"))}
  
  # Will need to add in lines removing taxa not of interest to the study i.e.
  # OTUs that were used for canopy level assignments but are not going to
  # have parameters estimated for.
  
  # Mnemonics to drop
  dropme <- TaxaTable$Mnemonic[which(as.logical(TaxaTable$DropParamtrz))]
  
  # Drop them
  for(i in 1:length(AllStems)){
    AllStems[[i]] <- AllStems[[i]][!(AllStems[[i]]$Mnemonic %in% dropme),]}
  
  # --- Stature data ----------------------------------------------------------
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
  
  # --- PARAMETER FOR SIZE/STATURE --------------------------------------------
  # For each unique mnemonic find the largest 6 individual live stems
  # and average their DBH.
  # future options could be using the X-percentile (95? or 90?)
  stature.df %<>% group_by(., Mnemonic) %>%
    summarise(., statureDBH = mean(sort(maxDH, decreasing = T)[1:6],
                                   na.rm = T),  number = n())
  
  # Add back in the taxa that couldn't be estimated.
  stature.df %<>% right_join(., TaxaTable_2est, by = "Mnemonic")
  
    
  # generate diagnostic/summary plots for the data
  # Generate Log histogram of statures
  plot_out <- ggplot(data = stature.df, aes(x = statureDBH)) +
    geom_histogram(breaks = 10^seq(0,2.5, by = 0.2) ) +
    scale_x_log10() +
    scale_y_continuous() + 
    xlab("Stature (cm)") +
    ylab("Number of OTUs") +
    ggtitle(paste0(site, " log-sizes")) + 
    theme_bw()
  # Save the output plot
  print(plot_out)
  ggsave(paste0(loc_out1, "summaryplots/",
                "StatureHistLog_", site, ".png"),
         plot = plot_out,
         device = "png",
         width = plot_w,
         height = plot_h,
         units = plot_u)
  
  # Generate histogram of statures
  plot_out <- ggplot(data = stature.df, aes(x = statureDBH)) +
    geom_histogram(breaks = seq(0, 200, by = 5) ) +
    scale_x_continuous() +
    scale_y_continuous() + 
    xlab("Stature (cm)") +
    ylab("Number of OTUs") +
    ggtitle(paste0(site, " sizes")) + 
    theme_bw()
  # save the histogram
  print(plot_out)
  ggsave(filename = paste0(loc_out1, "summaryplots/",
                           "StatureHist_", site, ".png"),
         plot = plot_out,
         device = "png",
         width = plot_w,
         height = plot_h,
         units = plot_u)

  # Plot to show the number of live stems vs the stature parameter
  plot_out <- ggplot(data = stature.df,
         aes(x = statureDBH, y = number, label = Mnemonic)) +
    geom_point() +
    geom_text(hjust = -0.1, vjust = -0.1) +
    scale_y_log10() +
    ylab("Number of stems") + 
    scale_x_log10(lim = c(1, 200)) +
    xlab("Stature (cm)") +
    ggtitle(paste0(site, " stature vs count")) +
    theme_bw()
  # Save the output plot
  print(plot_out)
  ggsave(filename = paste0(loc_out1, "summaryplots/",
                           "StatureVCount_", site, ".png"),
         plot = plot_out,
         device = "png",
         width = plot_w,
         height = plot_h,
         units = plot_u)
  
  # --- data for other parameters ---------------------------------------------
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
  
  # Calculations for individual stems
  # --- getting individual data for growth ------------------------------------
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
                    (!(Stem == "Secondary") | is.na(Stem)) & # only primary stems/single stemed
                    !is.na(DBH_bef) & # must have a DBH at both t_f and t_i
                    !is.na(DBH_aft))
    
    # Actual estimation
    x %<>% mutate(growth1 = case_when(DBH_use ~ DBH_bef*(
      (DBH_aft/DBH_bef)^(t_rf/Date_dur) - (DBH_aft/DBH_bef)^(t_ri/(Date_dur))),
      TRUE ~ NA_real_),
      growth2 = case_when(DBH_use ~ (DBH_aft/DBH_bef)^(1/Date_dur),
                            TRUE ~ NA_real_)) %>%
      select(c("StemID", "TreeID", "Mnemonic", "Census", "Stem", "CanopyLvl",
               "Date_dur" , "growth1", "growth2", "Codes_bef", "Codes_aft",
               "DBH_bef", "DBH_aft"))  
    
    return(x)
    
  })
  
  # --- Making plots for growth -----------------------------------------------
  # plotting growth-1 parameters
  for(c in 1:length(growth.l)){
    use_data <- growth.l[[c]] %>%
      filter(., !is.na(growth1)) %>%
      mutate(CanopyLvl = as.factor(CanopyLvl))
      
    lwr <- -0.25
    hgr <- 1
    
    # Keeping track of which ones got dropped off 
    nremv <- length(which(use_data$growth1 < lwr | use_data$growth1 > hgr))
    census_name <- unique(use_data$Census)
    
    # create the plot for output
    plot_out <- ggplot(data = use_data,
                       aes(growth1, fill = CanopyLvl, colour = CanopyLvl)) +
      geom_histogram(binwidth = .02, alpha = 0.25, position = "identity" ) +
      scale_x_continuous(name ="Growth-I (cm)") +
      scale_fill_hue(name = "Canopy level") +
      scale_color_hue(name = "Canopy level") +
      coord_cartesian(xlim = c(lwr, hgr)) +
      theme_bw() +
      ylab("Count") +
      xlab("Growth-I (cm)") +
      theme(plot.caption = element_text(hjust=0, margin = margin(15,0,0,0)),
            legend.title = element_text()) +
      ggtitle(paste0(site,"-", census_name, " growth-I histogram"))
    
    
    # If there are values out of range, make sure to have it be recorded
    if(nremv > 0 ){
      plot_out <- plot_out + labs(caption = paste0("* ", nremv, " values are not represented."))}
    
    # Save the file as a png
    print(plot_out)
    ggsave(filename = paste0(loc_out1, "summaryplots/",
                             "growth1Hist_", site, census_name, ".png"),
           plot = plot_out,
           device = "png",
           width = plot_w,
           height = plot_h,
           units = plot_u)
    
    
    
    
    
    
    
    
    
    
    
 
}
  
  # plotting growth-2 parameters
  for(c in 1:length(growth.l)){
    use_data <- growth.l[[c]] %>%
      filter(., !is.na(growth1)) %>%
      mutate(CanopyLvl = as.factor(CanopyLvl),
             growth2 = growth2)
    
    lwr <- 0.85
    hgr <- 1.35
    
    # Keeping track of which ones got dropped off 
    nremv <- length(which(use_data$growth2 < lwr | use_data$growth2 > hgr))
    census_name <- unique(use_data$Census)
    
    # create the plot for output
    plot_out <- ggplot(data = use_data,
                       aes(growth2, fill = CanopyLvl, colour = CanopyLvl)) +
      geom_histogram(binwidth = 0.01, alpha = 0.25, position = "identity" ) +
      scale_x_continuous(name ="Growth-II", labels = scales::percent) +
      scale_fill_hue(name = "Canopy level") +
      scale_color_hue(name = "Canopy level") +
      coord_cartesian(xlim = c(lwr, hgr)) +
      theme_bw() +
      ylab("Count") +
      xlab("Growth-I (cm)") +
      theme(plot.caption = element_text(hjust = 0, margin = margin(15,0,0,0)),
            legend.title = element_text()) +
      ggtitle(paste0(site, "-", census_name, " growth-II histogram"))
    
    # If there are values out of range, make sure to have it be recorded
    if(nremv > 0 ){
      plot_out <- plot_out + labs(caption = paste0("* ", nremv, " values are not represented."))}
    
    # Save the file as a png
    print(plot_out)
    ggsave(filename = paste0(loc_out1, "summaryplots/",
                             "growth2Hist_", site, census_name, ".png"),
           plot = plot_out,
           device = "png",
           width = plot_w,
           height = plot_h,
           units = plot_u)
  }
  
  
  
  
  # --- getting individual survival data --------------------------------------
  survive.l <- lapply(condensed.L, function(x){
    
    # Criteria for consideration
    x %<>% mutate(surv_use =
                    RametStatus_bef == "A" & # stem found alive at t_i
                    RametStatus_aft %in% c("A", "D") & # stem found at t_f
                    Stem != "Secondary")               # Stem is primary
    
    # Generate the value
    x %<>% mutate(survived = case_when(RametStatus_bef == "A" &
                                 RametStatus_aft == "A" ~ as.logical(1),
                               RametStatus_bef == "A" &
                                 RametStatus_aft == "D" ~ as.logical(0),
                               TRUE ~ NA))
    
    # cut the unnecessary columns
    x %<>% select(., c("StemID", "TreeID", "Mnemonic",
                       "Census", "CanopyLvl", "Date_dur",
                      "survived"))
    
    })
  
  # --- Diagnostic plots for survival -----------------------------------------
  survive.l2 <- survive.l
  for(c in 1:length(survive.l)){
    use_data <- survive.l[[c]] 
    # Which census is here
    census_name <- unique(use_data$Census)
    
    use_data %<>%
      filter(!is.na(survived)) %>%
      mutate(mnem_CnpyL = paste0(Mnemonic, "-", CanopyLvl),
             survived = as.numeric(survived)) %>%
      group_by(., mnem_CnpyL) %>%
      summarise(., n = n(),
                n_surv = sum(survived),
                Mnemonic = unique(Mnemonic)[1],
                CanopyLvl = unique(CanopyLvl)[1])  %>% 
      mutate(p_surv = n_surv/n, CanopyLvl = as.factor(CanopyLvl))
    
    use_data2 <- use_data %>% group_by(., Mnemonic) %>%
      summarize(., n = sum(n),n_surv = sum(n_surv), p_surv = sum(n_surv)/sum(n))
    
    
    
    survive.l2[[c]] <- left_join(TaxaTable_2est, use_data2, by = "Mnemonic")
    
    
    plot_out <- ggplot(data = use_data,
           aes(x = n, y = p_surv, label = Mnemonic, colour = CanopyLvl)) + 
      scale_color_discrete(name = "Canopy level") +
      geom_point() +
      scale_x_log10() +
      ylab("Proportion survived") +
      xlab("Number of stems") +
      geom_text(hjust = - 0.1, vjust = - 0.1, show.legend = F) +
      ggtitle(paste0(site, "-", census_name, " survival vs. stem count")) +
      theme_bw()
    
    # Save the plot
    print(plot_out)
    ggsave(filename = paste0(loc_out1, "summaryplots/",
                             "PpnSurv_CnpyL_", site, census_name, ".png"),
           plot = plot_out,
           device = "png",
           width = plot_w,
           height = plot_h,
           units = plot_u)
    
    # Make the same plot without the canopy levels
    plot_out <- ggplot(data = use_data2,
                       aes(x = n, y = p_surv, label = Mnemonic)) + 
      geom_point() +
      scale_x_log10() +
      ylab("Proportion survived") +
      xlab("Number of stems") +
      geom_text(hjust = - 0.1, vjust = - 0.1, show.legend = F) +
      ggtitle(paste0(site, "-", census_name, "survival vs. stem count")) +
      theme_bw()
    
    # Save the plot
    print(plot_out)
    ggsave(filename = paste0(loc_out1, "summaryplots/",
                             "PpnSurv_all_", site,
                             census_name, ".png"),
           plot = plot_out,
           device = "png",
           width = plot_w,
           height = plot_h,
           units = plot_u)
    
    
    
  }
  
  # --- getting recruitment data ----------------------------------------------
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
    n_rect.df <- x %>% filter(., RecruitStatus %in% c("GR", "RR")) %>%
      mutate(dummy_GR = as.numeric(RecruitStatus == "GR")) %>% 
      group_by(., Mnemonic) %>% summarise(., n_GR = sum(dummy_GR), n_recrt = n()) 
    
    # add in those that had 0 recruitment
    n_rect.df <- TaxaTable_2est %>%
      mutate(n_GR = 0, n_recrt = 0) %>%
      filter(!(Mnemonic %in% n_rect.df$Mnemonic)) %>% rbind(n_rect.df, .)
    
   
    
    
    x %<>% filter(RametStatus_bef == "A") # stem was alive at t_i
    
    # calculate recuitment rates
    # recrt_am <- recruits per basal area (stems above DBH adult cutoff) at t1
    # recrt_ae <- recruits per basal area of all conspecifics alive at t1
    # recrt_sm <- recruits per stem (stems above DBH adult cutoff) at t1
    # recrt_se <- recruits per stem of all conspecifics alive at t1
    
    x_all <- x %>% group_by(Mnemonic) %>%
      summarise(n_stem_all = n(),
                n_gnt_all = n_distinct(TreeID),
                area_stem_all = sum( pi*(((DBH_bef/100)/2)^2), na.rm = T))
    
    x_mat <- x %>% left_join(., adult_cutoff, by = "Mnemonic") %>%
      filter(., DBH_bef >= cutoff) %>% group_by(Mnemonic) %>%
      summarise(n_stem_mat = n(),
                n_gnt_mat = n_distinct(TreeID),
                area_stem_mat = sum(pi*(((DBH_bef/100)/2)^2), na.rm = T))
    
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

    # Generate plots
    # plot of # recruits +1 / number of mature stems
    plot_out <- ggplot(data = n_rect.df, aes(recrt_p1_stemsperMInd)) +
       geom_histogram(bins = 20) +
       scale_x_log10(name = "(New stems + 1)/indiv (mature)") +
       ylab("Count") +
       ggtitle(paste0(site, "-", census_name, " per-adult recruitment")) +
       theme_bw()
     
    # Save the file
    print(plot_out)
     ggsave(filename = paste0(loc_out1, "summaryplots/",
                              "recrt_perMIndiv_hist_",
                              site, census_name, ".png"),
            plot = plot_out,
            device = "png",
            width = plot_w,
            height = plot_h,
            units = plot_u)
     
     # repeat but with using the basal arera as the denominator instead.
     plot_out <- ggplot(data = n_rect.df, aes(recrt_p1_stemsperMBA)) +
       geom_histogram(bins = 20) +
       scale_x_log10(name = "(New stems + 1)/m^2 (mature)") +
       ylab("Count") +
       ggtitle(paste0(site, "-", census_name," per-basal area recruitment")) +
       theme_bw()
     
     # Save the file
     print(plot_out)
     ggsave(filename = paste0(loc_out1, "summaryplots/",
                              "recrt_perMBA_hist_",
                              site, census_name, ".png"),
            plot = plot_out,
            device = "png",
            width = plot_w,
            height = plot_h,
            units = plot_u)
     
     
     return(n_rect.df)
     
     })
  
  
  # --- Output ----------------------------------------------------------------
  loc_2 <- "/data/ForestGEO/processed/useValues/"
  # Not the best solution to use assign, but I wanted to make sure I could
  # individually retrieve data for each site.
  
  Vals_out <- list(site = site,
                   n_census = n_census,
                   taxa = TaxaTable_2est,
                   stature = stature.df,
                   growth = growth.l,
                   survival = survive.l,
                   recruitment = recruit.l)
  
  assign(paste0(site,"_AllStems_Vals.L"), Vals_out)
  save(list = paste0(site,"_AllStems_Vals.L"),
       file = paste0(loc_Gdr, loc_2, site, "_AllStems_Vals.r"))
  
  
}



