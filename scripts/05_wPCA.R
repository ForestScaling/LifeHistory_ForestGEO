# 06_wPCA.R
# ### Setup ###################################################################
set.seed(123)
# read in Packages
for(x in c("tidyverse",
           "readr",
           "dplyr",
           "magrittr",
           "lmodel2",
           "smatr",
           "ggplot2",
           "RColorBrewer")){
  library(x, character.only = T)
}

use_FunTrt <- F

# set locations
loc_Gdr <- "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats"
loc_data <- paste0(loc_Gdr, "/data/ForestGEO/processed/useValues/")
loc_data2 <- paste0(loc_Gdr, "/data/ForestGEO/taxa/")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")
cat(paste0("Set working directory for scripts at:\n", getwd(),"\n"))

# load the PCA code from Kambach
source("./scripts/wPCA_code.R")

# Load in the NEON foliar trait data
NEON.df <- read.csv(paste0(loc_Gdr,
                           "/data/FunctionalTraits/",
                           "NEON_FoliarTraits.csv"))

# Logit functions
logit_inv <- function(x){
  return( 1/(exp(-1*x) + 1) )
}

logit <- function(p){
  return(log(p/(1 - p)))
}

# --- Site specific setup -----------------------------------------------------
# Load in the site's data
site <- "HVDF"
for(site in c("HVDF","SCBI","SERC")){
# Data for PCA
load(paste0(loc_Gdr, "/data/ForestGEO/processed/param_matrix/", # path
            site, "_vals_wPCA.r")) # file name

# load the taxa table
load(paste0(loc_data2, site, "_TaxaTable.r"))
TaxaTable <- Spp_table






brokenstick <- function(n, xth = 1:n){
  # Try catching problems
  if(!(n %% 1  == 0)){
    stop("n isn't an integer.")}
  if(!((is.logical(xth) & length(xth) == n) |
     all(xth %% 1 == 0) & all(xth > 0 & xth <= n))){
    stop(paste0("xth must either be a logical vector of length n,",
         " or  be integers between 0 and n inclusive."))
  }
  
  # Calculations
  pieces <- vector(mode = "numeric", length = n) 
  
  for(i in 1:n){
    pieces[i] <- sum(1/(n - 0:(n-i)))/n
    }

  return(pieces[xth])
}
  
# ### Data separation ##########################################################
# Overview of how much data is there:
smy_byParam <- data_All %>% # By Parameter, i.e. % of species with a "real" val
  group_by(param) %>%
  summarize(prop_w_estimates = sum(as.integer(param_est))/n()) 
smy_byMnemo <- data_All %>% # By species/OTU, i.e. % of "real" vals per species
  arrange(Mnemonic) %>%
  group_by(Mnemonic) %>% 
  summarize(prop_w_estimates = sum(as.integer(param_est))/n())
# OTUs to keep, must have at least k% of their parameters actuall estimated
# other parameters recieve the "avg", with a very low weight assigned.

# Cutoff for including the OTU
k <- 0.75

keep_OTU <- smy_byMnemo$Mnemonic[which(smy_byMnemo$prop_w_estimates >= k)]

# run though different parameter combinations:
# for(data_All in list(data_All,
#                      filter(data_All, !grepl("growsd", x = param)),
#                      filter(data_All, !(grepl("growBin", x = param)|
#                             grepl("growsd", x = param))))){

  




# Number of parameters (to check that each OTU has all its params listed)
n_params <- length(unique(data_All$param))


data_All <- data_All %>% # drop the other OTUs
  filter(Mnemonic %in% keep_OTU) %>%
  arrange(Mnemonic)

data_All$val[data_All$param %in%
               c("recrt_p1_stemsperMBA", "recrt_p1_stemsperMInd")] <-
  log(data_All$val[data_All$param %in%
                     c("recrt_p1_stemsperMBA", "recrt_p1_stemsperMInd")])


# verify all species have all parameters
if(nrow(data_All)/length(keep_OTU) != n_params){
  stop(paste0("data_All's, Number of rows does not correspond to the number",
              " of OTU x number of parameters!"))
}

# --- Retrieve the two matricies ----------------------------------------------
# Weights
W_min <- 10^(-6)
W.df <- data_All %>%
  select(-c(param_est, val)) %>%
  pivot_wider(., names_from = param, values_from = W)

# Vals 
Vals.df <- data_All %>%
  select(-c(param_est, W)) %>%
  pivot_wider(., names_from = param, values_from = val)



# Add in functional Trait data

FuncTrt.df <- Spp_table %>% select(c(Latin, Mnemonic)) %>%
  left_join(NEON.df, by = "Latin") %>%
  select(-c(Source, Latin)) %>%
  filter(Mnemonic %in% Vals.df$Mnemonic)

# Add in the trait data
if(use_FunTrt){
  cat("Adding functional trait data from NEON, gives missing species minimum",
      "weight and average value...")
Vals.df <- FuncTrt.df %>%
  mutate(across(-c(Mnemonic),
                ~case_when(is.na(.x) ~ mean(.x, na.rm = T),
                           TRUE ~ .x))) %>%
  right_join(Vals.df, by = "Mnemonic")

W.df <- FuncTrt.df %>%
  mutate(across(-c(Mnemonic),
                ~case_when(is.na(.x) ~ W_min,
                           TRUE ~ 1))) %>%
  right_join(W.df, by = "Mnemonic")

cat("done.\n")
}

W.df <- W.df %>% arrange(., Mnemonic)
Vals.df <- Vals.df %>% arrange(., Mnemonic)
keep_OTU <- sort(keep_OTU)


#if(site != "HVDF"){stop()}
# check that the names are lined up:
if(!identical(W.df$Mnemonic, keep_OTU) | !identical(Vals.df$Mnemonic, keep_OTU)){
  stop("Mnemonic names not all aligned!")
}else if(!identical(colnames(W.df), colnames(Vals.df))){
  stop("Parameter columns names not all aligned!")
}else{
  # Generate the matricies as such, keep the rows and columns in order!
  OTU_names <- Vals.df$Mnemonic
  Vals.df %<>% select(-c(Mnemonic))
  W.df %<>% select(-c(Mnemonic))
  Param_names <- colnames(Vals.df)
}

n_params <- ncol(W.df)

# Run the wPCA
wPCA_res <- PCA.weighted(dat.for.PCA = Vals.df, weights.for.PCA = W.df)

#

# diagnostics
# histogram of Euclidean distances
dist_diag.plt <-
  wPCA_res$pca_coordinates %>%
  dist(x = ., method = "euclidean") %>%
  data.frame(EucDist = .) %>%
  mutate(EucDist = as.numeric(EucDist)) %>%
  ggplot(., mapping = aes(x = EucDist)) +
  geom_histogram(mapping = aes(y = after_stat(density)), bins = 25) +
  geom_density(linewidth = 1) +
  ggtitle(paste0("wPCA_EucDist_", site)) +
  theme_bw() +
  ylab("Density") + 
  xlab("Euclidean distance") +
  scale_y_continuous(expand = expansion(mult = c(0,.1)))
print(dist_diag.plt)

brkn_stick <- brokenstick(n = n_params)
# screeplot 
scree_diag.plt <- data.frame(PC = 1:n_params,
           ExplnVar = wPCA_res$expl_var,
           BrknStick = brkn_stick) %>%
  pivot_longer(cols = -PC,
               names_to = "Variance",
               values_to = "Ppn_var") %>%
  ggplot(aes(x = PC, y = Ppn_var, group = Variance, color = Variance)) +
  geom_line() +
  geom_point() +
  ylab("Proportion of variance") +
  scale_x_continuous(n.breaks = n_params) +
  scale_color_manual(labels = c("Broken-stick", "Explained"),
                     values = c("navy","firebrick4")) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
  ggtitle(paste0("Scree_", site)) +
  theme_bw()
print(scree_diag.plt)

# Index of loadings 
IndexofLoadings.f <- function(loadings, eigenvals){
  return(t(t(as.matrix(loadings)^2)*eigenvals^2))
}

IL <- IndexofLoadings.f(as.matrix(wPCA_res$factor_loadings[,-1]),
                                  wPCA_res$expl_var)

# Function for running the permutation test
wPCA_perm.f <- function(data, weights){
  # check that data dim = weight dims
  if(any(dim(data) != dim(weights))){
    stop("Dimensions of data and weights don't match.")
    }
  
  perm_data <- data # output for the data
  perm_w <- weights # output for the weights
  
  for(c in 1:ncol(data)){
    perm_pos <- sample(1:nrow(data),
                       size = nrow(data),
                       replace = F)
    
    perm_data[,c] <- data[perm_pos, c]
    perm_w[,c] <- weights[perm_pos, c]
    
  }
  
  # Run the weighted PCA
  wPCARes <- PCA.weighted(dat.for.PCA = perm_data,
                          weights.for.PCA = perm_w)
  
  return(wPCARes)
}

# Run the permutations
n_perm <- 999
perm_res.mat <- matrix(NA, nrow = n_perm, ncol = ncol(Vals.df))
names(perm_res.mat) <- names(wPCA_perm.f(Vals.df, W.df)$expl_var)

# Permutations for index of loadings
perm_IL.a <- array(NA, dim = c(nrow(IL), ncol(IL), n_perm + 1))



# permutations
for(perm in 1:n_perm){
  # all of the permutation results
  perm_res <- wPCA_perm.f(Vals.df, W.df) 
  # permutation for the explained variance
  perm_res.mat[perm, ] <- perm_res$expl_var 
  # permutation ofr index of the loadings
  perm_IL.a[,,perm+1] <- IndexofLoadings.f(
    as.matrix(perm_res$factor_loadings[,-1]),
    perm_res$expl_var)

}



# Add in the real data
perm_res.mat <- rbind(perm_res.mat, wPCA_res$expl_var) # for explained var
perm_IL.a[,,1] <-IL

# one tailed p-value for IL
out.a <-array(NA, dim(perm_IL.a))
for(j in 1:dim(perm_IL.a)[3]){
  out.a[,,j] <- perm_IL.a[,,j] >= IL
}

# Permutations for the permutation of 
out.a <- apply(out.a,c(1,2), FUN = sum)/dim(out.a)[3] # get p-value

# format results
out.df <- out.a %>% as.data.frame()
# format names 
rownames(out.df)<- unlist(wPCA_res$factor_loadings[,1])
colnames(out.df) <- colnames(wPCA_res$factor_loadings[,-1])

cat("One-tailed permutation p-value for Index of loadings, n_perm:", n_perm,"\n")
print(out.df)


# Get the one tailed p-value for explained variance
for(i in 1:ncol(Vals.df)){
  perm_res.mat[,i] <- as.numeric(perm_res.mat[,i] >= wPCA_res$expl_var[i])
}
cat("One-tailed permutation p-value for explained variance:\n")
out.v <- perm_res.mat %>%
  apply(., MARGIN = 2, FUN = function(X)return(sum(X)/(n_perm + 1)))
print(out.v)



  
# transform loadings and permutation p_value to a single df for plotting
# pivot the df
out.df <- out.df %>%
  mutate(., var = rownames(.)) %>%
  pivot_longer(cols = -c(var),
               names_to = "PC",
               values_to = "perm_p_value")
# join with the loadings
out.df <- wPCA_res$factor_loadings %>%
  rename(var = demographic_rate) %>%
  pivot_longer(cols = -c(var),
               names_to = "PC",
               values_to = "loading") %>%
  left_join(out.df, ., by = c("var","PC")) 

loading.plt <- out.df %>%
  mutate(var = as.factor(var)) %>%
mutate(PC = fct_relevel(PC, unique(out.df$PC)),
       var = fct_relevel(var, levels(var)[nlevels(var):1]))  %>%
  filter(PC %in% paste0("PCA", 1:5)) %>%
ggplot(.,aes(x = PC, y = var, fill = log10(perm_p_value))) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "white", mid = "white", midpoint = -1) +
  geom_text(aes(label = round(loading,2))) +
  ggtitle(paste(site, "Permutation of loadings")) +
  labs(caption = paste0("Loadings of variables to the different PCs, ",
       "values denote loading, color denotes p-value for \none tailed ",
       "permutation of Index of loading (n_permutations = ",n_perm,"). ",
       "Using permutation test of \nexplained variance per PC, only the ",
       "following PC were signficant (a = 0.05):\n",
       paste(names(out.v[which(out.v <= 0.05)]), collapse = ", ")))


# perm_dat <-Vals.df
# perm
# PCA.weighted(dat.for.PCA = Vals.df, weights.for.PCA = W.df)
# 
# 
#   
# Scale for plotting
scaling <- "1"
if(scaling == "2"){
  G <- wPCA_res$pca_coordinates %>%
    as.matrix # matrix of PC-scores
  E <- wPCA_res$factor_loadings[,-1] %>%
    as.matrix # matrix of trait score
  G %*% t(E) # returns the CENTERED data matrix
  
  # Norms
  normsG <- sqrt(diag(crossprod(G)))
  normsE <- sqrt(diag(crossprod(E)))
  
  G2 <- sweep(G, 2, normsG, "/") %>%
    as_tibble()
  
  E2 <- sweep(E, 2, normsG, "*") %>%
    as_tibble() %>%
    mutate(demographic_rate = wPCA_res$factor_loadings$demographic_rate, .before = 1)
  
  wPCA_res$pca_coordinates <- G2
  wPCA_res$factor_loadings <- E2
}


scl <- max(abs(wPCA_res$pca_coordinates$PCA1),
           abs(wPCA_res$pca_coordinates$PCA2))
scl <- ceiling(scl*2)/2


param_plotstuff <- data.frame(param = Param_names) %>% 
                              mutate(set = case_when(
                                grepl(pattern = "growparam_sp_",
                                      x = Param_names) ~ "growth_avg",
                                grepl(pattern = "godds_rat",
                                      x = Param_names) ~ "growth_p",
                                grepl(pattern = "sodds_rat",
                                      x = Param_names) ~ "surv_p",
                                grepl(pattern = "recrt",
                                      x = Param_names) ~ "recrtment",
                                grepl(pattern = "stature",
                                      x = Param_names) ~ "stature",
                                grepl(pattern = "growsd_",
                                      x = Param_names) ~ "growth_sd",
                                TRUE ~ "Other")) %>%
  mutate(set = as.factor(set)) %>%
  arrange(param)

param_plotstuff$color <- NA

# colors for the different groups
param_plotstuff$color[param_plotstuff$param == "log_stature"] <-
  "black"

param_plotstuff$color[param_plotstuff$set == "growth_sd"] <-
  c("#FFC15E","#F7B05B","#F7934C", "#CC5803")[1:3]

param_plotstuff$color[param_plotstuff$set == "surv_p"] <-
  c("#FF2C55","#E2294F","#7D1128", "#C41E3D","#3C0919")[1:3]

param_plotstuff$color[param_plotstuff$set == "_"] <-
  c("#858AE3","#613DC1","#4E148C","#2C0735")[1:3]

param_plotstuff$color[param_plotstuff$set == "recrtment"] <-
  c("#06BEE1","#1768AC","#2541B2", "#03256C")[2:3]

param_plotstuff$color[param_plotstuff$set == "growth_avg"] <-
  c("#90A955","#4F772D","#31572C", "#132A13")[1:3]

param_plotstuff$color[param_plotstuff$set == "growth_p"] <-
  c("#DBFEB8","#C5EDAC","#99C2A2")[1:3]

param_plotstuff$color[param_plotstuff$set == "Other"] <-
  c("#343d46","#4f5b66","#65737e","#a7adba","#c0c5ce")
# pickcentral <- function(n, v, upper = T){
#     l <- length(v)
#     d <- l - n
#     
#     v[-1*c(1:(d/2), ceiling(1+l - (d/2)):l)]
#     
#   }
  
  
wPCA_res
  
lambda_k <- apply(wPCA_res$factor_loadings[,-1], 2, FUN = function(X){sum(X^2)})
 # for scaling loading  

wPCA_plot <- ggplot(data = wPCA_res$pca_coordinates, aes(x = PCA1, y = PCA2)) + 
  annotate(geom = "text",
           x = wPCA_res$pca_coordinate$PCA1,
           y = wPCA_res$pca_coordinate$PCA2,
           label = OTU_names) +
  theme_bw() +
  scale_x_continuous(limits = c(-1*scl, scl),
                     name = ifelse(
                       wPCA_res$expl_var[1] > brokenstick(n_params, 1),
                                  paste0("wPC1* (",
                                         100*round(wPCA_res$expl_var[1],3),
                                         "% Variation)"),
                                  paste0("wPC1 (",
                                         100*round(wPCA_res$expl_var[1],3),
                                         "% Variation)"))) +
  scale_y_continuous(limits = c(-1*scl, scl),
                     name = ifelse(
                       wPCA_res$expl_var[1:2] > brokenstick(n_params, 1:2),
                                   paste0("wPC2* (",
                                          100*round(wPCA_res$expl_var[2],3),
                                          "% Variation)"),
                                   paste0("wPC2 (",
                                          100*round(wPCA_res$expl_var[2],3),
                                          "% Variation)"))) +
  geom_segment(data = wPCA_res$factor_loadings,
               aes(x = 0, y = 0,
                   xend = PCA1*scl, yend = PCA2*scl,
                   color = demographic_rate), linewidth = 1.5, arrow = arrow()) +
  scale_color_manual(values = param_plotstuff$color) +
  # annotate("path",
  #          x=xc+r*cos(seq(0,2*pi,length.out=100)),
  #          y=yc+r*sin(seq(0,2*pi,length.out=100)))
  ggtitle(paste0(site,"-wPCA"))

#print(wPCA_plot)

library(gridExtra)
grid.arrange(wPCA_plot, loading.plt, nrow = 1)

}




# data_plot <- data_All %>%
#   filter(param_est) %>% select(-c(W)) %>%
#   pivot_wider(., names_from = param, values_from = val)
# 
# 
# 
# sma_res <- lmodel2(data = data_plot,
#                    log_stature ~ recrt_p1_stemsperMBA
#                    )$regression.results
# sma_res %<>% filter(Method == "SMA")
# ggplot(data = data_plot,
#        aes(x = recrt_p1_stemsperMBA, y = log_stature))+
# annotate(geom = "text", x = data_plot$recrt_p1_stemsperMBA,
#          y = data_plot$log_stature, label = data_plot$Mnemonic) + theme_bw()+
#   geom_abline(data = sma_res, aes(intercept = Intercept, slope = Slope))
# 
# data_plot <- data_All %>%
#   filter(param_est) %>% select(-c(W)) %>%
#   pivot_wider(., names_from = param, values_from = val)
# 
# 
# sma_res <- lmodel2(data = data_plot,
#                    growparam_sp_CL1 ~ recrt_p1_stemsperMBA)$regression.results
# 
# sma_res %<>% filter(Method == "SMA")
# ggplot(data = data_plot,
#        aes(x = recrt_p1_stemsperMBA, y = growparam_sp_CL1))+
#   annotate(geom = "text", x = data_plot$recrt_p1_stemsperMBA,
#            y = data_plot$growparam_sp_CL1, label = data_plot$Mnemonic) + theme_bw()+
#   geom_abline(data = sma_res, aes(intercept = Intercept, slope = Slope))
# 
# data_plot <- data_All %>%
#   filter(param_est) %>% select(-c(W)) %>%
#   pivot_wider(., names_from = param, values_from = val)
# 
# # Growth vs survival
# sma_res <- lmodel2(data = data_plot,
#                    growparam_sp_CL1 ~ surv_pCL1
# )$regression.results
# sma_res %<>% filter(Method == "SMA")
# ggplot(data = data_plot,
#        aes(x = surv_pCL1, y = growparam_sp_CL1))+
#   annotate(geom = "text", x = data_plot$surv_pCL1,
#            y = data_plot$growparam_sp_CL1, label = data_plot$Mnemonic) + theme_bw()+
#   geom_abline(data = sma_res, aes(intercept = Intercept, slope = Slope))
# 
# sma_res <- lmodel2(data = data_plot,
#                    growparam_sp_CL1 ~ logit(surv_pCL1)
# )$regression.results
# sma_res %<>% filter(Method == "SMA")
# ggplot(data = data_plot,
#        aes(x = logit(surv_pCL1), y = growparam_sp_CL1))+
#   annotate(geom = "text", x = logit(data_plot$surv_pCL1),
#            y = data_plot$growparam_sp_CL1, label = data_plot$Mnemonic) + theme_bw()+
#   geom_abline(data = sma_res, aes(intercept = Intercept, slope = Slope))