
#         Title: 03c SDM model prediction
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão (adapted Wiedenroth, 2023)
#         Last update: 23.07.2025
#         Usage: This script uses combined occurrence/enviro data
#         "final_subset" (02b), and the trained models (03a), to predict a set
#         of SDMs for each species, for every year.
#         Algorithms: MaxEnt.
#         Assumptions: The species grid sizes are expressed in a ".csv" file in 
#         the main Analysis folder.
#         Output: Habitat suitability predictions for each species.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# Load packages and functions
source("scripts/00b_functions.R")

# data folder
envdatafolder <- 
  file.path("/new_nfs_shared/data/04_PROCESSED_DATA/coregistration/NRW")

# results folder
resultsfolder <- file.path("results")

# set output_dir to 03c
output_dir <- file.path(resultsfolder, "03c_model-prediction")

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

# species grid sizes
gridsizes <- read.csv(file = "sppgridsizes.csv",header = F)

# loop over all species and predict habitat suitability
for (h in 1: length(pa_env_files)) { 
  
#==============================================================================#
####                               Load data                                ####
#==============================================================================#
  
  i <- pa_env_files[h]
  
  # load env-occ-df
  occ <- loadRDS(i)
  
  # species name
  spp <- occ$species[1]
  
  if (length(model_files[grep(spp, model_files)]) != 0){
    
    cat("\n")
    cat(paste("Start predicting habitat suitability for selected species: ",spp,sep=""))
    cat("\n")
  
    if (spp == sort(gridsizes,"V1")[h,1]){
    
      grid_size <- sort(gridsizes,"V1")[h,2]
    }
  
   
    # load in the environmental data (both static and dynamic variables)  # adapt as necessary
    envdata_static_files <- 
      list.files(file.path(envdatafolder,
                           paste0("coregistration_2016_2024_",grid_size,
                                  "m/static")),full.names = TRUE, pattern = ".tif")
    envdata_dynamic_files <- 
      list.files(file.path(envdatafolder,
                           paste0("coregistration_2016_2024_",grid_size,
                                  "m/dynamic")),full.names = TRUE, pattern = ".tif")
    envdata_static <- terra::rast(envdata_static_files)
    envdata_dynamic <- terra::rast(envdata_dynamic_files)
  
    # rename (static environmental) layers
    current_names <- names(envdata_static)
    new_names <- sub("_NRW","", current_names)
    new_names[23] <- "tree_div"
    names(envdata_static) <- new_names
    rm(current_names, new_names)
  
    # set NA values for "tree_spec", "tree_div" & "canopyheight" to 0
    envdata_static$tree_spec[is.na(envdata_static$tree_spec)] <- 0
    envdata_static$tree_div[is.na(envdata_static$tree_div)] <- 0
    envdata_static$canopyheight[is.na(envdata_static$canopyheight)] <- 0
    
    # fix categorical variable
    envdata_static <- envdata_cat(envdata_static, 9)
    envdata_static <- envdata_cat(envdata_static, 11)
    envdata_static <- envdata_cat(envdata_static, 21)
    envdata_static <- envdata_cat(envdata_static, 24)
  
    # select year for prediction
    unique_years <- seq(2016,2024,1)
    
    for (selected_year in unique_years){
      
      cat(paste0("    Year of prediction: ",selected_year))
      cat("\n")
  
      envdyn_index <- grep(selected_year, names(envdata_dynamic))
      envdyn_select <- subset(envdata_dynamic, envdyn_index)
  
      # rename variables
      current_names <- names(envdyn_select)
      new_names <- sub(paste0("_",selected_year),"", current_names)
      names(envdyn_select) <- new_names
      names(envdyn_select)[16] <- "soil_moist"
      rm(current_names, new_names)
  
      envdata <- c(envdata_static,envdyn_select)
  
      # set reference system
      projcrs = crs(envdata)
      
      # load models
      load(model_files[grep(spp, model_files)])
  
      # get final predictor variables (pred_sel)
      pred_sel <- names(maxent_sm$samplemeans)
  
      # load model performances
      load(model_perf_files[grep(spp, model_perf_files)])
  
      # environmental dataframe to predict for
      preddata <- data.frame(crds(envdata[[pred_sel]]), 
                             as.points(envdata[[pred_sel]]))
  
      # free up some memory
      gc()
  
  
      # set weights for the presences as 1. The weights of the pseudo absences are
      # calculated so that their sum equals the sum of the weights of the presences
      wgt <- ifelse(occ$occurrence == 1, 1, 
                  sum(occ$occurrence == 1) / sum(occ$occurrence == 0 ))
  
  
#==============================================================================#
####                           Model prediction                             ####
#==============================================================================#
  
      cat("    Predictions for study area")
      cat("\n")
  
      # make predictions of all models
      pred <- data.frame(preddata[, c("x", "y")],
                         maxent = predict(maxent_sm, preddata,
                                          type = "logistic"))
  
      gc()
  
      # save predictions
      save(pred, file = paste0(output_dir, "/03c_", spp,"_MaxEnt_",
                               selected_year,"_pred.RData"))
  
      cat("     > RData file saved")
      cat("\n")
  
      # make a SpatRaster from predictions
      r_pred <- rast(pred, crs = projcrs)
  
      # save SpatRaster
      writeRaster(r_pred, filename = paste0(output_dir,"/03c_",spp,
                                          "_MaxEnt_",selected_year,
                                          "_pred.tif"),overwrite = TRUE)
  
      cat("     > GeoTIFF file saved")
      cat("\n")
      cat("    Predictions for occurrence locations")
      cat("\n")
  
      # fix variable names
      new_names <- sub("_0-5cm_mean","", names(occ)[3:ncol(occ)])
      new_names <- sub("_0-30cm_mean","", new_names)
      occ <- occ %>% rename_with(~ new_names, .cols = 3:ncol(occ))
      
      # make predictions of all models for occurrence locations
      pred_occ <- data.frame(maxent = predict(maxent_sm, occ[, pred_sel],
                                            type = "logistic"))
  
      gc()
  
      # save occurrence predictions
      save(pred_occ, file = paste0(output_dir,"/03c_",spp,
                                 "_MaxEnt_",selected_year,"_pred_occ.Rdata"))
  
      cat("     > RData file saved")
      cat("\n")
  
      # convert predictions to binary output using optimized threshold
      cat("    Computing binary predictions (optimized threshold)")
      cat("\n")
      
      pred_bin_ot <- data.frame(preddata[, c("x", "y")],
                             sapply(names(pred[-c(1:2)]), FUN=function(alg){
                               
                               ifelse(pred[,alg] >= 
                                        maxent_perf_simp$'Thresh_opt', 1, 0)}))
      
      gc()
      
      # make spat raster from binary predictions
      r_pred_bin_ot <- rast(pred_bin_ot, crs = projcrs)
      
      # save SpatRaster
      writeRaster(r_pred_bin_ot, filename = paste0(output_dir,"/03c_", spp,
                                                "_MaxEnt_",selected_year,
                                                "_pred_bin_ot.tif"),
                  overwrite = TRUE)
      
      cat("     > GeoTIFF file saved")
      cat("\n")
    }
  }
}  # end of loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
