
#         Title: 03b SDM model validation
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão (adapted Wiedenroth, 2023)
#         Last update: 23.07.2025
#         Usage: This script uses combined occurrence/enviro data 
#         "final_subset" (02b), spatial blocks "scv" (02a), and the trained 
#         models (03a) to validate a set of SDMs for each species. Calculates
#         two prediction thresholds: (1) maximum TSS; (2) optimizes threshold
#         by incorporating (averaging) both max TSS and pseudo-absences weights.
#         Algorithms: MaxEnt.
#         Output: performance metrics for each algorithm and every species. The 
#         threshold-dependent statistics (TSS, Kappa, Sens, Spec & PCC) 
#         delivered, relate to the maxTSS threshold.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# load packages and functions
source("scripts/00b_functions.R")

# results folder
resultsfolder <- file.path("results")

# list occ-env-df
pa_env_files <- list.files(file.path(resultsfolder,
                                     "02b_collinearity/env-occ-df"),
                           full.names = TRUE, pattern = ".rds")

# files spatial blocks
spatial_blocks_files <- list.files(file.path(resultsfolder,
                                             "02a_spatial-blocks"),
                                   full.names = TRUE, pattern = ".rds")

# full models
full_model_files <- list.files(file.path(resultsfolder, "03a_model-training"),
                          full.names = TRUE, pattern = "full")

# final models
simp_model_files <- list.files(file.path(resultsfolder, "03a_model-training"),
                              full.names = TRUE, pattern = "simp")

# set output_dir
output_dir <- file.path(resultsfolder, "03b_model-validation")

# loop over all the species and validate species distribution models
for(i in pa_env_files) { 
  
#==============================================================================#
####                               Load data                                ####
#==============================================================================#
  
  # load env-occ-df
  occ <- loadRDS(i)
  
  # species name
  spp <- occ$species[1]
  
  cat("\n")
  cat(paste("Start validating SDMs for selected species: ",spp,sep=""))
  cat("\n")

  # load spatial blocks (scv)
  scv <- loadRDS(spatial_blocks_files[grep(spp, spatial_blocks_files)])
  
  # load models
  full_model <- load(full_model_files[grep(spp, full_model_files)])
  
  # load full predictor variables (pred_sel)
  pred_sel_full <- names(maxent_m$samplemeans)
  
  # set condition for simplified model
  if(length(grep(spp, simp_model_files))==0){
    
    cat("    No simplified model for this species")
    cat("\n")
  } else{
    
    load(simp_model_files[grep(spp, simp_model_files)])
    
      # get simplified model predictor names
      pred_sel_simp <- names(maxent_sm$samplemeans)
  }
  
  # set weights for the presences as 1. The weights of the pseudo absences are
  # calculated so that their sum equals the sum of the weights of the presences
  wgt <- ifelse(occ$occurrence == 1, 1, 
                sum(occ$occurrence == 1) / sum(occ$occurrence == 0 ))
  
  # fix variable names
  new_names <- sub("_0-5cm_mean","", names(occ)[3:ncol(occ)])
  new_names <- sub("_0-30cm_mean","", new_names)
  occ <- occ %>% rename_with(~ new_names, .cols = 3:ncol(occ))
  
  
#==============================================================================#
####                           Model validation                             ####
#==============================================================================#
  
 # Maxent validation -----------------------------------------------------------
  
  cat("    MaxEnt model validation\n")
  cat("     > full model validation\n")

  # calculate block cross-validated predictions and (full) model performance
  maxent_cv_pred_full = numeric(length = nrow(occ))
  
  for(i in seq_len(length(scv$folds_list))){
    
    cv_train <- occ[unlist(scv$folds_list[[i]][1]),]
    cv_test <- occ[unlist(scv$folds_list[[i]][2]),]
    
    # update the model for the new training data
    modtmp_full <- maxnet(p = cv_train$occurrence, 
                          data = cv_train[,pred_sel_full])
    
    # make predictions for k-fold:
    maxent_cv_pred_full[unlist(scv$folds_list[[i]][2])] <- 
      predict(modtmp_full, cv_test, type = "logistic")
    }
  
  # calculate model performance (maximizing TSS)
  maxent_perf_max <- evalSDM(occ$occurrence, maxent_cv_pred_full,
                              weights = wgt)
  
  # calculate boyce index
  maxent_boyce_index_full <- 
    ecospat.boyce(fit = maxent_cv_pred_full,
                  obs = maxent_cv_pred_full[which(occ$occurrence == 1)],
                  nclass = 0, window.w = "default", res = 100,
                  PEplot = FALSE, rm.duplicate = TRUE, method = "kendall")
  maxent_boyce_index_full <- maxent_boyce_index_full$cor
  
  # make predictions of all models for occurrence locations
  pred_occ_full <- data.frame(maxent = predict(maxent_m, occ[, pred_sel_full],
                                               type = "logistic"))
  
  # calculate model explained deviance
  expl_dev_full <- data.frame(maxent = expl_deviance(occ$occurrence,
                                                     pred_occ_full$maxent,
                                                     weights = wgt))
  
  # calculate optimized threshold, obtained by averaging the threshold that 
  # maximizes TSS and the one that weights FNC to pseudo-absences weights
  opt_thresh <- mean(c((evalSDM(occ$occurrence, maxent_cv_pred_full,
                                thresh.method='Cost',
                                FPC=1, FNC=unique(wgt)[2])$thresh),
                       (evalSDM(occ$occurrence, maxent_cv_pred_full,
                                thresh.method='MaxSens+Spec')$thresh)))
  
  # calculate model performance (optimized threshold)
  maxent_perf_opt <- evalSDM(occ$occurrence, maxent_cv_pred_full,
                             thresh=opt_thresh,weights = wgt)
  
  # summarize performances
  maxent_perf_full <- as.data.frame(matrix(0,1,15))
  names(maxent_perf_full) <- c("AUC","Boyce","D2","Thresh_max","TSS_max",
                               "Kappa_max","Sens_max","Spec_max","PCC_max",
                               "Thresh_opt","TSS_opt","Kappa_opt","Sens_opt",
                               "Spec_opt","PCC_opt")
  maxent_perf_full[,c(1,4:9)] <- maxent_perf_max[,c(1,8,2:6)]
  maxent_perf_full[,2] <- maxent_boyce_index_full
  maxent_perf_full[,3] <- expl_dev_full
  maxent_perf_full[,c(10:15)] <- maxent_perf_opt[,c(8,2:6)]
  
  # save model performances and cross-validated predictions
  save(maxent_perf_full, maxent_cv_pred_full,
       file = paste0(output_dir, "/03b_", spp,"_MaxEntfull_perf-cv-pred.RData"))
  
  cat("      > full model performances saved\n")
  
  # set condition for no simplified model
  if(length(grep(spp, simp_model_files))!=0){
    
    cat("     > simplified model validation\n")

    # calculate block cross-validated predictions and (simplified) model 
    # performance
    maxent_cv_pred_simp = numeric(length = nrow(occ))
  
    for(i in seq_len(length(scv$folds_list))){
    
      cv_train <- occ[unlist(scv$folds_list[[i]][1]),]
      cv_test <- occ[unlist(scv$folds_list[[i]][2]),]
    
      # update the model for the new training data
      modtmp_simp <- 
        maxnet(p = cv_train$occurrence, data = cv_train[,pred_sel_simp])
    
      # make predictions for k-fold:
      maxent_cv_pred_simp[unlist(scv$folds_list[[i]][2])] <- 
        predict(modtmp_simp, cv_test, type = "logistic")
      }

    # calculate model performance (maximizing TSS)
    maxent_perf_maxsimp <- evalSDM(occ$occurrence, maxent_cv_pred_simp,
                                weights = wgt)
  
    # calculate boyce index
    maxent_boyce_index_simp <- 
      ecospat.boyce(fit = maxent_cv_pred_simp, 
                    obs = maxent_cv_pred_simp[which(occ$occurrence == 1)], 
                    nclass = 0, window.w = "default", res = 100, 
                    PEplot = FALSE, rm.duplicate = TRUE, method = "kendall")
    maxent_boyce_index_simp <- maxent_boyce_index_simp$cor
    
    # make predictions of all models for occurrence locations
    pred_occ_simp <- 
      data.frame(maxent = predict(maxent_sm, occ[, pred_sel_simp],
                                  type = "logistic"))
    
    # calculate model explained deviance
    expl_dev_simp <- data.frame(maxent = expl_deviance(occ$occurrence,
                                                       pred_occ_simp$maxent,
                                                       weights = wgt))
    
    # calculate optimized threshold, obtained by averaging the threshold that 
    # maximizes TSS and the one that weights FNC to pseudo-absences weights
    opt_thresh <- mean(c((evalSDM(occ$occurrence, maxent_cv_pred_simp,
                                  thresh.method='Cost',
                                  FPC=1, FNC=unique(wgt)[2])$thresh),
                         (evalSDM(occ$occurrence, maxent_cv_pred_simp,
                                  thresh.method='MaxSens+Spec')$thresh)))
    
    # calculate model performance (optimized threshold)
    maxent_perf_optsimp <- evalSDM(occ$occurrence, maxent_cv_pred_simp,
                               thresh=opt_thresh,weights = wgt)
    
    # summarize performances
    maxent_perf_simp <- as.data.frame(matrix(0,1,15))
    names(maxent_perf_simp) <- c("AUC","Boyce","D2","Thresh_max","TSS_max",
                                 "Kappa_max","Sens_max","Spec_max","PCC_max",
                                 "Thresh_opt","TSS_opt","Kappa_opt","Sens_opt",
                                 "Spec_opt","PCC_opt")
    maxent_perf_simp[,c(1,4:9)] <- maxent_perf_maxsimp[,c(1,8,2:6)]
    maxent_perf_simp[,2] <- maxent_boyce_index_simp
    maxent_perf_simp[,3] <- expl_dev_simp
    maxent_perf_simp[,c(10:15)] <- maxent_perf_optsimp[,c(8,2:6)]
    
    # save model performances and cross-validated predictions
    save(maxent_perf_simp, maxent_cv_pred_simp,
         file = paste0(output_dir, "/03b_", spp,
                       "_MaxEntsimp_perf-cv-pred.RData"))
    
    cat("      > simplified model performances saved\n")
  }
}  # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#