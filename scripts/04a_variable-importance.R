
#         Title: 04a SDM model variable importance
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão
#         Last update: 23.07.2025
#         Usage: This script uses combined occurrence/enviro data (02b),
#         and the trained and validated models (03a & 03b), to calculate the 
#         model variable importance for each species' model.
#         Algorithms: MaxEnt.
#         Output: variable importance for each model in a ".rds" and a 
#         ".csv" file.


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

# final models
model_files <- list.files(file.path(resultsfolder, "03a_model-training"),
                          full.names = TRUE, pattern = "simp")

# model performances
model_perf_files <- list.files(file.path(resultsfolder, "03b_model-validation"),
                               full.names = TRUE, pattern = "simp")

# set output_dir to 04a
output_dir <- file.path(resultsfolder, "04a_variable-importance")

# loop over all species and predict habitat suitability
# due to the high number of permutations (100) per variable, it might be 
# advisable not to run in a loop but rather one species at a time
for(i in pa_env_files) {        
  
  
#==============================================================================#
####                               Load data                                ####
#==============================================================================#
  
  # load env-occ-df
  occ <- loadRDS(i)
  
  # species name
  spp <- occ$species[1]
  
  # set condition for no simplified model
  if(length(grep(spp, model_files))!=0){
  
    cat("\n")
    cat(paste0("Calculating model variable importance for selected species: ",
               spp))
    cat("\n")
  
    # load models
    load(model_files[grep(spp, model_files)])
  
    # load model performances
    load(model_perf_files[grep(spp, model_perf_files)])
    
    # set weights for the presences as 1. The weights of the pseudo absences are
    # calculated so that their sum equals the sum of weights of all presences
    wgt <- ifelse(occ$occurrence == 1, 1, 
                  sum(occ$occurrence == 1) / sum(occ$occurrence == 0 ))
    
    # fix variable names
    new_names <- sub("_0-5cm_mean","", names(occ))
    new_names <- sub("_0-30cm_mean","", new_names)
    pred_sel <- new_names
    occ <- occ %>% rename_with(~ new_names)
    
  
#==============================================================================#
####                   Variable importance calculation                      ####
#==============================================================================#  
  
    # calculate variable importance of the Maxent models through permutation of
    # each variable 100 times, and assessing its impact on model performance 
    # (deviance explained)
    num_permutations <- 100
  
    var_importance <- sapply(names(maxent_sm$samplemeans), function(var) {
      cat(paste0("    Permuting variable ", var, "\n"))
      importance_values <- numeric(num_permutations)
      for (perm in 1:num_permutations) {
        cat(paste0("     > Permutation ", perm, " of ", num_permutations, "\n"))
        permuted_env <- occ[,3:ncol(occ)]
        permuted_env[[var]] <- sample(permuted_env[[var]])
        permuted_pred <- predict(maxent_sm, permuted_env,
                                 type = "logistic")
        permuted_expl_deviance <- expl_deviance(occ[,2], permuted_pred, 
                                                weights = wgt)
        importance_values[perm] <- 
          maxent_perf_simp$D2 - permuted_expl_deviance
      }
      mean_importance <- mean(importance_values)
      return(mean_importance)
    })
  
    var_importance <- var_importance[order(var_importance, decreasing = TRUE)]  
  
    # scale variable importance values to 100
    perc_varImp <- 
      (var_importance/maxent_perf_simp$D2)/
      sum(var_importance/maxent_perf_simp$D2)
  
    # save model performance and cross validation predictions -----------------
    perc_varImp <- as.data.frame(perc_varImp)
  
    write.csv(perc_varImp,file = paste0(output_dir,"/04a_", spp,
                                        "_MaxEntsimp_varImp.csv"))
  }
}  # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
