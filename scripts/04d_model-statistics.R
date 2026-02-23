
#         Title: 04d SDM model statistics compilation
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão
#         Last update: 23.07.2025
#         Usage: This script uses combined occurrence/enviro data 
#         "final_subset" (02b), the trained and validated models (03a & 03b),
#         and the respective variable importance values (04a) to compile the
#         model descriptive statistic for each species.
#         Algorithms: MaxEnt. 
#         Output: Model descriptive statistics for all species.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# results folder
resultsfolder <- file.path("results")

# data folder
envdatafolder <- 
  file.path("/new_nfs_shared/data/04_PROCESSED_DATA/coregistration/NRW")

# list environmental data
envdata_static_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_100m/static"), 
             full.names = TRUE, pattern = ".tif")

envdata_dynamic_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_100m/dynamic"), 
             full.names = TRUE, pattern = ".tif")

# list occ-env-df
pa_env_files <- list.files(file.path(resultsfolder,
                                     "02b_collinearity/env-occ-df"),
                           full.names = TRUE, pattern = ".rds")

# full models
full_model_files <- list.files(file.path(resultsfolder, "03a_model-training"),
                               full.names = TRUE, pattern = "full")

# final models
simp_model_files <- list.files(file.path(resultsfolder, "03a_model-training"),
                               full.names = TRUE, pattern = "simp")

# full model performances
full_model_perf_files <- 
  list.files(file.path(resultsfolder, "03b_model-validation"),
             full.names = TRUE, pattern = "full")

# final model performances
simp_model_perf_files <- 
  list.files(file.path(resultsfolder, "03b_model-validation"),
             full.names = TRUE, pattern = "simp")

# full model variable importances
full_model_varImp_files <- 
  list.files(file.path(resultsfolder,"04a_variable-importance"),
             full.names = TRUE, pattern = "full")

# final model variable importances
simp_model_varImp_files <- 
  list.files(file.path(resultsfolder,"04a_variable-importance"),
             full.names = TRUE, pattern = "simp")

# set output_dir to 04b
output_dir <- file.path(resultsfolder, "04d_model-statistics")


#==============================================================================#
####                                Load data                               ####
#==============================================================================#

## get variable names
# load environmental data
envdata_static <- terra::rast(envdata_static_files)
envdata_dynamic <- terra::rast(envdata_dynamic_files)

# rename (static environmental) layers
current_names <- names(envdata_static)
new_names <- sub("_NRW","", current_names)
new_names[23] <- "tree_div"
names(envdata_static) <- new_names
rm(current_names, new_names)

# select latest dynamic environmental layers (2024)
envdyn_index <- grep("2024", names(envdata_dynamic))
envdyn_select <- subset(envdata_dynamic, envdyn_index)

# rename variables
current_names <- names(envdyn_select)
new_names <- sub("_2024","", current_names)
names(envdyn_select) <- new_names
names(envdyn_select)[16] <- "soil_moist"
rm(current_names, new_names)

varnames <- names(c(envdata_static,envdyn_select))
  
# create compilation of model performances accross species  
performances <- as.data.frame(matrix(0,2,length(pa_env_files)))
rownames(performances) <- c("AUCfull","AUCsimp")

# create compilation of full model variable contribution values accross species
modelconts_full <- 
  as.data.frame(matrix(0,length(varnames),length(pa_env_files)))
rownames(modelconts_full) <- varnames

# create compilation of final model variable contribution values accross species
modelconts_simp <- 
  as.data.frame(matrix(0,length(varnames),length(pa_env_files)))
rownames(modelconts_simp) <- varnames

for (h in 1: length(pa_env_files)){
  
  i <- pa_env_files[h]
  
  occ <- readRDS(i)
  
  cat("\n")
  cat("Retrieving files...")
  cat("\n")
  
  spp <- occ$species[1]
  
  cat(paste("Compiling model statistics for selected species: ",spp,sep=""))
  cat("\n")

  cat(" > Model performances")
  cat("\n")
  
  colnames(performances)[h] <- spp

  # extract full model performances
  load(full_model_perf_files[grep(spp, full_model_perf_files)])
  performances[1,h] <- maxent_perf_full$AUC
  
  # extract final model performances
  
  if(length(grep(spp, simp_model_perf_files)) == 0) {
    performances[2,h] <- 0
  } else {
    load(simp_model_perf_files[grep(spp, simp_model_perf_files)])
    performances[2,h] <- maxent_perf_simp$AUC
  }
    
  cat(" > Variable importance values")
  cat("\n")
  
  # assign species names
  colnames(modelconts_full)[h] <- spp
  
  # extract full model variable importance values
  full_model_varImp <- 
    read.csv(full_model_varImp_files[grep(spp, full_model_varImp_files)])
  
  for (j in 1:length(full_model_varImp$X)){
    
    modelconts_full[full_model_varImp$X[j],h] <- full_model_varImp$x[j]
  }
  
  # assign species names
  colnames(modelconts_simp)[h] <- spp
  
  # extract final model variable importance values
  if (length(grep(spp, simp_model_varImp_files)) == 0) {
    
    modelconts_simp[simp_model_varImp$X[j],spp] <- 0
  } else {
    
    simp_model_varImp <- 
      read.csv(simp_model_varImp_files[grep(spp, simp_model_varImp_files)])
  
    for (j in 1:length(simp_model_varImp$X)){
      
      modelconts_simp[simp_model_varImp$X[j],spp] <- 
      simp_model_varImp$perc_varImp[j]
      }
  }
} # end loop around each species

write.csv(performances, file = file.path(output_dir,
                                         "04d_Model-performances.csv"))
write.csv(modelconts_full, file = file.path(output_dir,
                                            "04d_Model-variables-full.csv"))
write.csv(modelconts_simp, file = file.path(output_dir,
                                            "04d_Model-variables-simp.csv"))

cat("\n")
cat("Model performances and variable importance values saved")
cat("\n")


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
