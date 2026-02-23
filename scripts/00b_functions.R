
#         Title: 00b_functions
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão  (adapted from Wiedenroth, 2023)
#         Last update: 31.03.2025
#         Usage: Includes all generic and specific functions for each script in 
#         the SDM workflow.
#         Assumptions: The main (analysis) project folder is the working 
#         directory; this script is located in the scripts folder
#         Run inside scripts with: source("scripts/00b_functions.R")


#==============================================================================#
####                      00 Package management functions                   ####
#==============================================================================#

### install and/or load packages
# Set up project environment and load packages
## install.load.package
install_load_package <- function(x) {
  
  if (!require(x, character.only = TRUE))
    install.packages(x, repos = 'http://cran.us.r-project.org')
  require(x, character.only = TRUE)
} # add names of the packages required as character objects

# name the required packages
packages <- c(
  "dplyr", "data.table", "terra", "broom", "tidyverse", "sf", "purrr", "furrr", 
  "sfheaders", "stringr", "stringi", "tidyr", "readxl", "R.filesets", "here", 
  "dismo", "doParallel", "foreach", "blockCV", "automap", "corrplot", "ecospat",
  "randomForest", "gbm", "mgcv", "maxnet", "viridis", "ggplot2", "gridExtra",
  "RColorBrewer", "cowplot", "tidyterra", "ggpubr")

# run the function to load the packages
sapply(packages, install_load_package)


#==============================================================================#
####                    01 Occurrence preparation functions                 ####
#==============================================================================#

### Spatial Thinning ###
# Load function for "fast" spatial thinning (e.g. 400m)
# runs and cores reduced to 1 as not running on HPC - adapt if necessary
## thin ##
thin <- function(sf, thin_dist = 400, runs = 1, ncores = 1){

  require(sf, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(furrr, quietly = TRUE)
  
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  sf_buffer <- st_buffer(sf, thin_dist)
  buff_int <- st_intersects(sf, sf_buffer) 
  buff_int <- setNames(buff_int, 1:length(buff_int))
  
  n_int <- map_dbl(buff_int, length)
  
  plan(multisession, workers = ncores)
  
  seeds <- sample.int(n = runs)
  results_runs <- future_map(seeds, function(i){
    
    set.seed(i)
    while (max(n_int) > 1) {
      
      max_neighbors <- names(which(n_int == max(n_int)))
      
      # remove point with max neighbors
      sampled_id <- sample.vec(max_neighbors, 1)
      
      pluck(buff_int, sampled_id) <- NULL
      buff_int <- map(buff_int, function(x) setdiff(x, as.numeric(sampled_id)))
      n_int <- map_dbl(buff_int, length)
    }
    
    unlist(buff_int) %>% unique()
    
  }, .options = furrr_options(seed = TRUE))
  
  lengths <- map_dbl(results_runs, length)
  
  selected_run <- results_runs[[sample.vec(which(lengths == max(lengths)), 1)]]
  
  out <- sf[selected_run,]
  
  out <- sf_to_df(out)[,3:4] %>%
    rename("lon" = "x", "lat" = "y")
  
  out
}


### Environmental data processing ###
# Sets the specified environmental predictor (n) as categorical
## envdata_cat ##
envdata_cat <- function(envdata,n){
  
  # Define n as categorical, while keeping NAs
  biotop_cat_raster <- envdata[[names(envdata)[n]]]
  biotop_cat_raster[is.nan(biotop_cat_raster)] <- NA
  values <- unique(values(biotop_cat_raster, na.rm = TRUE))
  levels(biotop_cat_raster) <- data.frame(id=values,
                                          category=as.character(values))
  envdata[[names(envdata)[n]]] <- biotop_cat_raster
  envdata[[names(envdata)[n]]][is.na(envdata[[names(envdata)[n]]])] <- NaN
  
  return(envdata)
}


#==============================================================================#
####                         02 Collinearity functions                      ####
#==============================================================================#

### Select non-correlated variables ###
# this is a block cross-validated test. It ends up with the predictors that are 
# not correlated
## select07_blockcv ##
select07_blockcv <- function(X, 
                             y, 
                             family = "binomial",
                             univar = "glm2", 
                             threshold = 0.7, 
                             method = "spearman",
                             sequence = NULL, weights = NULL, 
                             sp_block) {
  
  # selects variables based on removing correlations > 0.7, retaining those
  # variables more important with respect to y
  # Order of importance can be provided by the character vector 'sequence'
  
  # 1. step: cor-matrix
  # 2. step: importance vector
  # 3. step: identify correlated pairs
  # 4. step: in order of importance: remove less important of the correlated 
  #           variables, recalculate correlation matrix a.s.f.
  
  folds <- sp_block$folds_list
  
  weights.full <- weights
  
  ks <- sp_block$folds_ids
  
  imp <- apply(X, 2, 
               compute_univar_cv, 
               response = y, 
               family = family, 
               univar = univar, 
               ks = ks,
               weights = weights)
  
  cm <- cor(X, method = method)
  
  if (is.null(sequence)) {
    
    sort.imp <- colnames(X)[order(imp, decreasing = T)]
  } else {
    
    sort.imp <- sequence
  }
  
  pairs <- which(abs(cm) >= threshold, arr.ind = T) # identifies correlated variable pairs
  index <- which(pairs[,1] == pairs[,2])            # removes entry on diagonal
  pairs <- pairs[-index,]                           # -"-
  
  exclude <- NULL
  for (i in 1:length(sort.imp)){
    
    if ((sort.imp[i] %in% row.names(pairs)) &
        ((sort.imp[i] %in% exclude) == F)) {
      
      cv <- cm[setdiff(row.names(cm), exclude), sort.imp[i]]
      cv <- cv[setdiff(names(cv), sort.imp[1:i])]
      exclude <- c(exclude, names(which((abs(cv) >= threshold)))) 
    }
    
  }
  
  pred_sel <- sort.imp[!(sort.imp %in% unique(exclude)), drop = F]
  return( list(D2 = sort(imp, decreasing = T), 
               cor_mat = cm, 
               pred_sel = pred_sel)
  )
}


### Compute univariate models ###
# is used in select07_blockcv
## compute_univar_cv ##
compute_univar_cv <- function(variable, response, family, univar, ks, weights){
  preds <- numeric(length(response))
  univar_list <- list()  
  
  ## univar_select ----
  # is used in compute_univar.v
  univar_select <- function(m) {
    
    df <- data.frame(occ = response, env = variable)
    train_df <- df[!ks == m,]
    test_df <- df[ks == m,]
    if(length(unique(train_df[,2])) > 4) {
      
      univar <- "gam"
    }
    
    if(length(unique(train_df[,2])) < 5 & length(unique(train_df[,2])) > 2) {
      
      univar <- "glm2"
    }
    
    if(length(unique(train_df[,2])) < 3) {
      
      univar <- "glm1"
    }
    
    return(univar)
  }
  
  for (j in sort(unique(ks))) {
    
    univar_list[[j]] <- univar_select(m = j)
  }
  
  # the least complex model algorithm across all folds should always be used
  # for all folds. So if for one fold glm1 was needed we used glm1 for all 
  # folds.
  if(any(univar_list == "glm1")) {
    
    univar <- "glm1"
  } else {
    
    if(any(univar_list == "glm2")) {
      
      univar <- "glm2"
    } else {
      
      univar <- "gam" 
    }
    
  }
  
  for (n in sort(unique(ks))) {
    
    df <- data.frame(occ = response, env = variable)
    train_df <- df[!ks == n,]
    test_df <-  df[ks == n, ]
    
    m1 <- switch(univar,
                 glm1 = glm(occ ~ env, data = train_df, family = family, 
                            weights = weights[!ks == n]),
                 glm2 = glm(occ ~ poly(env, 2), data = train_df, 
                            family = family, weights = weights[!ks == n]),
                 gam = mgcv::gam(occ ~ s(env, k = 4), data = train_df, 
                                 family = family, weights = weights[!ks == n]))
    
    preds[ks == n] <- predict(m1, newdata = test_df, type = 'response')
    
    # print(summary(m1))
    
  }
  
  d2 <- expl_deviance(response, preds, weights = weights)
  ifelse(d2 < 0, 0, d2)
}



## Explained deviance ##
# Calculates the explained deviance based on the dismo package
## expl_deviance ##

expl_deviance <- function(obs, pred, family = 'binomial', 
                          weights = rep(1, length(obs))){
  
  if (family == 'binomial') {pred <- ifelse(pred < .00001, .00001, 
                                            ifelse(pred > .9999, .9999, pred))}
  
  null_pred <- rep(mean(obs), length(obs))
  
  1 - (dismo::calc.deviance(obs, pred, family = family, weights = weights) / 
         dismo::calc.deviance(obs, null_pred, family = family, 
                              weights = weights))
}


#==============================================================================#
####                       03 model validation functions                    ####
#==============================================================================#


### Evaluate SDM ###
# Evaluates SDM models by calculating different performance metrics
## evalSDM ##
evalSDM <- function(observation, predictions, thresh=NULL, 
                    thresh.method='MaxSens+Spec', req.sens=0.85, 
                    req.spec = 0.85, FPC=1, FNC=1,
                    weights=rep(1, length(observation))){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), obs = observation,
                           pred = predictions)
  
  if (is.null(thresh)) {
    
    thresh.mat <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, 
                                                      req.sens=req.sens, 
                                                      req.spec = req.spec, 
                                                      FPC=FPC, FNC=FNC)
    
    thresh <- thresh.mat[thresh.mat$Method==thresh.method,2]
  }
  
  cmx.opt <- PresenceAbsence::cmx(DATA= thresh.dat, threshold=thresh)
  
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
             TSS = TSS(cmx.opt), 
             Kappa = PresenceAbsence::Kappa(cmx.opt, st.dev=F),
             Sens = PresenceAbsence::sensitivity(cmx.opt, st.dev=F),
             Spec = PresenceAbsence::specificity(cmx.opt, st.dev=F),
             PCC = PresenceAbsence::pcc(cmx.opt, st.dev=F),
             D2 = expl_deviance(observation, predictions, weights=weights),
             thresh = thresh)
}


### True Skills Statistic ###
# Calculates the TSS for a specific model
## TSS ##
TSS = function(cmx){
  
  PresenceAbsence::sensitivity(cmx, st.dev=F) + 
    PresenceAbsence::specificity(cmx, st.dev=F) - 1
}


#==============================================================================#
####                            End of functions                            ####
#==============================================================================#
