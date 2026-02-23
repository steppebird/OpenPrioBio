
#         Title: 04b SDM partial dependency plots (response curves)
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão
#         Last update: 23.07.2025
#         Usage: This script uses combined occurrence/enviro data (02b), the
#         trained models (03a), and respective variable importances (04a) to 
#         create partial dependency plots (response curves).
#         Algorithms: MaxEnt.
#         Output: response curves for all variables of each model, ordered by
#         respective importances.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# Load packages and functions
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

# variable importances
var_imp_files <- list.files(file.path(resultsfolder, "04a_variable-importance"),
                            full.names = TRUE, pattern = "simp")

# set output_dir to 04b
output_dir <- file.path(resultsfolder, "04b_response-curves")

# loop over all species and predict habitat suitability
for(i in pa_env_files) { 
  
#==============================================================================#
####                               Load data                                ####
#==============================================================================#
  
  
  # load env-occ-df
  occ <- loadRDS(i)
  
  # species name
  spp <- occ$species[1]
  
  if (length(model_files[grep(spp, model_files)]) != 0){
    
    cat("\n")
    cat(paste0("Plotting response curves for selected species: ",spp))
    cat("\n")
    
    # load models
    load(model_files[grep(spp, model_files)])
    
    # load variable importances
    var_importance <- read.csv(var_imp_files[grep(spp, var_imp_files)])
    
    
#==============================================================================#
####                        Partial dependency plots                        ####
#==============================================================================#  
    
    varnames <- var_importance$X
    
    png(paste0(output_dir,"/04b_",spp,"_resp-curves.png"),
        width = 1200,height = 900)    # Increase width and height if necessary
    plot(maxent_sm, vars = varnames, common.scale = TRUE, type = "logistic")
    dev.off() 
  }
}  # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
