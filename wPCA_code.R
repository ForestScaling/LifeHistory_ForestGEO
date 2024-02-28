# wPCA_code.R
# Author: B. Rosenbaum, N. Ruger,S. Kambach
# Purpose: Includes basic functions for generating and plotting a wPCA

# 2. conduct weighted PCA -------------------------------------------------

PCA.weighted = function(dat.for.PCA,
                        weights.for.PCA){
  
  results = list()
  
  
  #####################################
  # Analyse demographic rates of x layers with a weighted PCA
  # Initial script by Benjamin Rosenbaum, Dr. Nadja Rueger (nadja.rueger@idiv.de)
  # restructured by Stephan Kambach (stephan.kambach@idiv.de)
  # see L. Delchambre. Weighted principal component analysis: a weighted covariance eigendecomposition approach. Mon. Not. R. Astron. Soc. 446(2), 3545-3555, 2014.
  #   data.estimates: table with rows=observations, cols = variables
  #   data.weights: table with the corresponding weights
  
  # scaled by plot-level sds
  X = dat.for.PCA
  W = weights.for.PCA
  
  X.save = X
  W.save = W
  
  X <- t(as.matrix(X))
  W <- t(as.matrix(W))
  
  # replace zero weights
  n <- as.numeric(ncol(X))
  d <- as.numeric(nrow(X))
  
  # substract weighted mean, save original data
  Xorig <- X
  
  # save sds and center for for back-transformation later
  centers = rowMeans(W*X)/rowMeans(W)
  sds.within.plot = centers
  
  X <- X - rowMeans(W*X)/rowMeans(W) # X = centered and scaled within plots
  
  # scale by weighted sd in each dimension
  # X sds per plot
  for (j in 1:d){
    sds.within.plot[j] <- sqrt(sum(W[j, ]^2*X[j, ]^2)/sum(W[j, ]^2))
    X[j, ] <- X[j, ] / sds.within.plot[j]}
  
  # weighted covariance matrix S
  S <- (W*X)%*%t(W*X) / (W%*%t(W))
  
  #  eigenvalues (in decreasing order) and corresponding eigenvectors (in columns)
  EV <- eigen(S)
  
  #----------------------------------------------------------------
  # step 2: loop through ranks to compute true explained variance
  #----------------------------------------------------------------
  Chi2 <- vector(length=d) # explained variance
  
  for(rank in 1:nrow(X)){
    
    # matrix of principal components (new basis)
    P <- EV$vectors[, 1:rank]
    
    # matrix of coefficients (coordinates in new basis)
    C <- matrix(data=NA, nrow=rank, ncol=n)
    
    for (i in 1:n){
      w <- diag(W[, i]^2)
      
      C[, i] <- solve(t(P)%*%w%*%P, t(P)%*%w%*%X[, i])}
    
    # matrix of projections of old data (PC approx X)
    PC <- P%*%C
    
    # explained variance
    Chi2[rank] <- sum(sum((W*PC)^2)) / sum(sum((W*X)^2))}
  
  # subtract expl variances to get the individual CHi2 of every principal component
  for(rank in 2:nrow(X)){
    Chi2[rank] = Chi2[rank] - sum(Chi2[1:(rank-1)])}
  
  #----------------------------------------------------------------
  # step 3: decide for a rank based on the explained variance
  #         and compute data for output (same steps as above)
  #----------------------------------------------------------------
  rank <- nrow(X)
  P <- EV$vectors[, 1:rank]
  C <- matrix(data=NA, nrow=rank, ncol=n)
  
  for (i in 1:n){
    w <- diag(W[,i]^2)
    C[, i] <- solve(t(P)%*%w%*%P, t(P)%*%w%*%X[, i])}
  
  PC <- P%*%C
  
  # re-order the principal components from largest to smallest explained variance
  P = P[,order(Chi2, decreasing = T)]
  C = C[order(Chi2, decreasing = T),]
  Chi2 = Chi2[order(Chi2, decreasing = T)]
  
  # Output is a list with 6 entries:
  #   chi2: explained variance
  #   principal_components:    matrix of principal components (new basis)
  #   pca_coordinates:    matrix of coefficients (coordinates in new basis)
  #   projections_old_data:   matrix of projections of old data (PC approx X)
  #   centers
  #   sds
  
  # re-arrange data for easier later analysis
  names(Chi2) = paste("PCA", 1:length(Chi2), sep = "")
  
  P = data.frame(names(dat.for.PCA), P)
  names(P) = c("demographic_rate", paste("PCA", 1:length(Chi2), sep = ""))
  P = tibble(P)
  
  C = data.frame(t(C))
  names(C) = c(paste("PCA", 1:length(Chi2), sep = ""))
  C = tibble(C)
  
  results =
    list("expl_var" = Chi2,
         "factor_loadings" = P,
         "pca_coordinates" = C,
         "centers" = centers,
         "sds" = sds.within.plot,
         "raw_data" = X.save,
         "raw_weights" = W.save)
  
  return(results)
}

# 3. plot weighted PCA ----------------------------------------------------
plot.weighted.PCA = function(dat.temp){
  plot.temp = ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(data = dat.temp$pca_coordinates, aes(x = PCA1, y = PCA2), alpha = 0.3) +
    geom_segment(data = dat.temp$factor_loadings, aes(x = 0, y = 0, xend = PCA1 *5, yend = PCA2*5),
                 arrow = arrow(length = unit(0.5,"cm"))) +
    geom_label(data = dat.temp$factor_loadings, aes(x = PCA1 * 3, y = PCA2 * 3, label = demographic_rate)) +
    xlab(paste0("PCA1 (", round(dat.temp$expl_var[1] * 100, 2), "%)")) +
    ylab(paste0("PCA2 (", round(dat.temp$expl_var[2] * 100, 2), "%)")) +
    theme_bw() + theme(panel.grid = element_blank())
  
  return(plot.temp)
}
