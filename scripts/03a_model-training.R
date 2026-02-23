
#         Title: 03a SDM model training workflow
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão (adapted Wiedenroth, 2023)
#         Last update: 22.07.2025
#         Usage: This script uses combined occurrence/enviro data 
#         "final_subset" (02b), spatial blocks "scv" (02a), and "pred_sel" (02b)
#         to train a set of SDMs for each species. Calculates variable 
#         importance for this model, by permutation 100 times each variables and
#         assessing its impact on the deviance explained (all importances are
#         scaled to 100). Afterwards the model is simplified to include only
#         variables with an importance greater than random (1/# of variables).
#         Algorithms: MaxEnt. 
#         Output: initial and simplified model objects for each species,
#         and variable importance values for initial (full) model.


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

# final predictor variables
pred_sel_files <- list.files(file.path(resultsfolder, "02b_collinearity"),
                             pattern = ".rds", full.names = TRUE)

# set output_dir to 03a
output_dir <- file.path(resultsfolder, "03a_model-training")

# loop over all the species and train species distribution models
for(i in pa_env_files) { 
  
  
#==============================================================================#
####                               Load data                                ####
#==============================================================================#
  
  # load env-occ-df
  occ <- loadRDS(i)
  
  # species name
  spp <- occ$species[1]
  
  cat("\n")
  cat(paste("Building SDMs for selected species: ",spp,sep=""))
  cat("\n")
  
  # load final predictor variables (pred_sel)
  pred_sel <- loadRDS(pred_sel_files[grep(spp, pred_sel_files)])
  
  # set weights for the presences as 1. The weights of the pseudo absences are
  # calculated so that their sum equals the sum of the weights of the presences
  wgt <- ifelse(occ$occurrence == 1, 1, 
                sum(occ$occurrence == 1) / sum(occ$occurrence == 0 ))
  
  # fix variable names
  new_names <- sub("_0-5cm_mean","", pred_sel)
  new_names <- sub("_0-30cm_mean","", new_names)
  pred_sel <- new_names
  occ <- occ %>% rename_with(~ new_names, .cols = 3:ncol(occ))
  
  
#==============================================================================#
####                    Species Distribution Models                         ####
#==============================================================================#
  
## Maxent ------------------------------------------------------------------
  cat("    Algorithm selected: MaxEnt")
  cat("\n")
  
  maxent_m <- maxnet(p = occ$occurrence, data = occ[,pred_sel])
  
  cat("    Initial model built")
  cat("\n")
  
  # save models
  save(maxent_m,file = paste0(output_dir,"/03a_",spp, "_MaxEnt_SDMfull.RData"))
  
  cat("    Initial model saved")
  cat("\n")  
  cat("    Calculating deviance explained")
  cat("\n")  
  
  # make predictions of all models for occurrence locations
  pred_occ <- data.frame(maxent = predict(maxent_m, occ[, pred_sel],
                                          type = "logistic"))
  
  # calculate deviance explained
  dev_expl <- data.frame(maxent = expl_deviance(occ$occurrence, pred_occ$maxent,
                                                weights = wgt))
  
  
#==============================================================================#
####                   Variable importance calculation                      ####
#==============================================================================#  
  
  cat("    Calculating variable importances")
  cat("\n")  
  
  # calculate variable importance of the Maxent models through permutation of
  # each variable 100 times, and assessing its impact on model performance 
  # (deviance explained)
  num_permutations <- 100
  
  var_importance <- sapply(names(occ[,3:ncol(occ)]), function(var) {
    
    cat(paste0("     > Permuting variable ", var, "\n"))
    importance_values <- numeric(num_permutations)
    for (perm in 1:num_permutations) {
      
      cat(paste0("      > Permutation ", perm, " of ", num_permutations, "\n"))
      permuted_env <- occ[,3:ncol(occ)]
      permuted_env[[var]] <- sample(permuted_env[[var]])
      permuted_pred <- predict(maxent_m, permuted_env,
                               type = "logistic")
      permuted_expl_deviance <- expl_deviance(occ[,2], permuted_pred,
                                              weights = wgt)
      importance_values[perm] <- dev_expl - permuted_expl_deviance
    }
    mean_importance <- mean(unlist(importance_values))
    return(mean_importance)
  })
  
  var_importance <- var_importance[order(var_importance, decreasing = TRUE)]  
  
  # scale variable importance values to 100
  perc_varImp_full <-
    (var_importance/(unlist(dev_expl)))/sum(var_importance/(unlist(dev_expl)))
  
  # save variable importance values for the full model
  write.csv(perc_varImp_full,
            file = paste0(resultsfolder, "/04a_variable-importance/04a_", spp,
                          "_MaxEntfull_varImp.csv"))
  
  cat("    Variable importance values saved")
  cat("\n")  
  
  
#==============================================================================#
####                           Model simplification                         ####
#==============================================================================#  
  
  cat("    Simplifying model")
  cat("\n")  
  
  # selecting variables
  numvars <- length(var_importance)
  criticalvalue <- 1/numvars
  selvars <- perc_varImp_full >= criticalvalue
  varnames <- names(perc_varImp_full)[selvars]
  
  if (length(varnames)<= 1){
    
    cat("     > Two few variables in initial modlel. No simplified model created")
    cat("\n")  
  } else{
    
    maxent_sm <- maxnet(p = occ$occurrence, data = occ[,varnames])
    
    cat("    Fitted a model with",length(varnames),"variables")
    cat("\n")  
    
    # save models
    save(maxent_sm,file = paste0(output_dir,"/03a_",spp, "_MaxEnt_SDMsimp.RData"))
    
    cat("    Simplified model saved")
    cat("\n")  
  }
} # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
