
#         Title: 04c SDM extract response values for categorical variable
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão
#         Last update: 02.06.2025
#         Usage: This script uses combined occurrence/enviro data (02b), the
#         trained models (03a) to extract response values for a specific 
#         categorical variable.
#         Algorithms: MaxEnt.
#         Output: response value for variable of choice.


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

# set output_dir to 04c
output_dir <- file.path(resultsfolder, "04c_response-values")

# define variable to extract response value from
sel_pred <- "tree_spec"

# loop over all species and predict habitat suitability
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
    cat(paste0("Extracting response values for ",spp," to ",sel_pred))
    cat("\n")
  
    # load models
    load(model_files[grep(spp, model_files)])

    # set condition for lulc being a selected variable
    if (sel_pred %in% names(maxent_sm$samplemeans)){
  
#==============================================================================#
####                         Extract response values                        ####
#==============================================================================#
  
      # Generate response curves
      response_values <- response.plot(maxent_sm, sel_pred, type="logistic",
                                     plot = F)
      write.csv(response_values, file = paste0(output_dir,"/04c_",spp,"_",
                                               sel_pred,"_response-values.csv"))
    } else {
      
      cat(" > Variable is not included in simplified model")
      cat("\n")
    }
  }
}  # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
