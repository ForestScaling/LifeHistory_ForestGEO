# Script: 00_PackageInstall
# Author: JMR
# Last update: 12023-09-12
# Description: 
# Installing necessary packages for the project. Should check if they already exist, will not
# go through and redownload them again.

# List of packages to use:
use_packages <- c("ggplot2",
                  "rlist",
                  "tidyverse","cowplot",
                  "nlme", "ggpubr", "fitdistrplus",
                  "magrittr","rlist", "rstudioapi","colormap",
                  "stringr", "MASS", "lubridate", "fitdistrplus", "kableExtra",
                  "lme4")

cur_packages <- as.data.frame(installed.packages())$Package # vector of current packages
add_packages <- use_packages[!(use_packages %in% cur_packages)] # check which to add

install.packages(add_packages)

# End Script