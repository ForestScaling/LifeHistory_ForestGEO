# File: 01A_TALLOAllometries.R
# Author: JMR
# Last Update: 12023-01-06
# Description: Take data from the TALLO database and output the allometries
# used to estimate the crown radius for analysis.
# ### Set up and data import ##################################################
# Packages
library("tidyverse")
library("lme4")
library("nlme")
library("magrittr")

# Setting/finding location of the script, so that the WD is automatically
#within GITHub local copy location.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

# Finding the location of the Gdrive for the Tallo data
loc_Gdr <-
  "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats/data/ForestGEO/TALLO"

# Bring in the data from TALLO (Jucker et al. 12022 Global Change Biology
# doi: 10.1111/gcb.16302)
Tallo_data <- read.csv(paste0(loc_Gdr, "/Tallo.csv"),
                       header = T, na.strings = c("NA"))
Tallo_envr <- read.csv(paste0(loc_Gdr,"/Tallo_environment.csv"),
                       header = T, na.strings = c("NA"))

# Combine into a single dataset
Tallo_data <- inner_join(Tallo_data, Tallo_envr, by = "tree_id")

# ### Models ##################################################################
# We are going to roughly follow Jucker et al. 12022 Global Change Biology,
# start with removing the outliers of DBH (we don't care about height
# allometries for this analysis so we will leave those outliers in).
Tallo_data <- Tallo_data[-which(Tallo_data$crown_radius_outlier == "Y"), ]

# Remove entries with missing data for crown radius (CR) or lack taxonomic
# assignment to at least the large clade
Tallo_data <- Tallo_data[-which(is.na(Tallo_data$division) |
                                  is.na(Tallo_data$crown_radius_m)),]

# add a variable for biome-division combination
Tallo_data <- Tallo_data %>%
  mutate(biome_div = paste(biome, division, sep = "_"))

# --- Allometries Mixed effects Method 1 --------------------------------------
# Fit the mixed-effects model
fit <- lmer(log(crown_radius_m) ~ log(stem_diameter_cm) +
              (log(stem_diameter_cm) | biome_div),
                   data = Tallo_data)

# Multiplying it by the correction factor? I don't think this is technically
# correct since this is the correction factor using the SEE of the entire model
# and their are definitely more than two parameters estimated here in this
# situation but this is just temporary.
# Note that we used log_e, so we can just use the equation at the bottom of 209
# in Sprugel 11983, Ecology
k <- ncol(coef(fit)$biome_div)* nrow(coef(fit)$biome_div)
n <- length(summary(fit)$residuals) # each observation should have a residual
SEE <- (sum((summary(fit)$residuals)^2)/(n - k))^0.5
CF <- exp((SEE^2)/2)

# Adding in the estimates for crown size, we will do two groups 
Tallo_rslt_LME <- list(fit = fit, CF = CF)

# --- Method 2, Separately linear models for each biome-division --------------
forest_types <- c("Temperate broadleaf forest_Angiosperm",
                  "Temperate broadleaf forest_Gymnosperm",
                  "Temperate conifer forest_Gymnosperm",
                  "Temperate conifer forest_Angiosperm")

forest_types_short <- c("TempBL_A", "TempBL_G", "TempCon_A","TempCon_G")

TalloLM_sep_fits <- list(rep(NA, length(forest_types)))



for(i in 1:length(forest_types)){
  
  f <- forest_types[i]
  
  # Fitlering the data, and fitting the model
  dat <- Tallo_data %>% filter(biome_div == f)
  fit <- dat %>% lm(data = .,
                    formula = log(crown_radius_m) ~ log(stem_diameter_cm))
  
  plot(y = log(dat$crown_radius_m), x = log(dat$stem_diameter_cm), main = f)
  
  # Calculating the correction factor
  k <- 2
  n <- length(summary(fit)$residuals)
  SEE <- (sum((summary(fit)$residuals)^2)/(n - k))^0.5
  CF <- exp((SEE^2)/2)
  
  plt <- dat %>%
    ggplot(data = ., aes(x = log(stem_diameter_cm), y = fit$residuals)) +
    geom_point(size = 1.25, shape = 1) +
    ggtitle(f) + 
    scale_x_continuous(name = "DBH in Ln(m)") +
    scale_y_continuous(name = "Residuals in Ln(m)") + theme_bw()
      
  print(plt)
  
  # Coefficients
  b0 <- fit$coefficients[1]
  b1 <- fit$coefficients[2]
  
  TalloLM_sep_fits[[i]] <- list(type = f, fit = fit, n = n,
                                CF = CF, b0 = b0, b1 = b1)
  names(TalloLM_sep_fits)[i] <- forest_types_short[i]
}

# ### Trying to find the best model for each biome-division ###################





# ### Outputs #################################################################
# Output the list of the results from Tallo
save(Tallo_rslt_LME, file = paste0(loc_Gdr,"/Tallo_MEModelFits.r"))
save(TalloLM_sep_fits, file = paste0(loc_Gdr,"/Tallo_sepLMModelFits.r"))














#b1 <- coef(fit)$biome_div$`log(stem_diameter_cm)`
#b0 <- coef(fit)$biome_div$`(Intercept)`
#bx <- rownames(coef(fit)$biome_div)

#i <- 3199
#j <- which(bx == Tallo_data$biome_div[i])
#predict(fit)[i]
#b1[j]*log(Tallo_data$stem_diameter_cm[i]) + b0[j] 




#smy <- summary(fit)
#sum((smy$residuals)^2)



# compare to each one estimated separately
#df1 <- Tallo_data[which(Tallo_data$biome_div == "Boreal/montane forest_Gymnosperm"),]
#lm(data = df1, log(crown_radius_m ) ~ log(stem_diameter_cm ))











