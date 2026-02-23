
#         Title: 00c SDM workflow
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão
#         Last update: 20.05.2025
#         Usage: Runs all scripts for the full SDM workflow
#         Assumptions: Modules are independent, meaning it can be run partially,
#         as long as all necessary inputs are available


#==============================================================================#
####                          00 Workspace set up                           ####
#==============================================================================#

## set workspace----------------------------------------------------------------

# set working directory
setwd("~/Analysis")

# create folder structure
# sets up the necessary folder structure for running the analyses
# only to be ran once, at the beginning of the project
source("scripts/00a_folder-structure.R")


#==============================================================================#
####                     01a Occurrence data preparation                    ####
#==============================================================================#

### prepare occurrence data ###
# description: cleans unnecessary fields, filters data by species and year (if
# necessary), filters by spatial accuracy, removes duplicates, sorts 
# pseudo-absences (if necessary), and does spatial thinning
# input: obsp (*.shp); envdata (*.tif)
# output: occ (01a_spp.rds)
# notes: split into 3 different grains of analysis (100m, 200m, & 1000m)

# occurrence data at 100m
source("scripts/01a_occurrence-data_prep_100m.R")

# occurrence data at 200m
source("scripts/01a_occurrence-data_prep_200m.R")

# occurrence data at 1000m
source("scripts/01a_occurrence-data_prep_1000m.R")


#==============================================================================#
####                    01b Environmental data extraction                   ####
#==============================================================================#

### extract environmental data ###
# description: extracts environmental data layers for each occurrence location
# input: envdata (*.tif); occ (01a_spp.rds); envpreds (spp_vars.csv)
# output: occ_env (01b_spp_enviro.rds)
# notes: split into 3 different grains of analysis (100m, 200m & 1000m)

# extraction at 100m
source("scripts/01b_enviro-data_extraction_100m.R")

# extraction at 200m
source("scripts/01b_enviro-data_extraction_200m.R")

# extraction at 1000m
source("scripts/01b_enviro-data_extraction_1000m.R")


#==============================================================================#
####                           02a Spatial blocks                           ####
#==============================================================================#

### define spatial blocks ###
# description: generates spatial blocks blocks for creating environmentally
# separated folds for k-fold cross-validation purposes
# input: envdata (*.tif); occ_env (01b_spp_enviro.rds)
# output: spatial-blocks (02a_spp_spatial-blocks.rds)
# plots: block size plot (spp_blocksize.png); spatial-blocks
# (spp_spatialblocks.png)
# notes: split into 3 different grains of analysis (100m, 200m & 1000m)

# Spatial blocks at 100m
source("scripts/02a_spatial-blocks_100m.R")

# Spatial blocks at 200m
source("scripts/02a_spatial-blocks_200m.R")

# Spatial blocks at 1000m
source("scripts/02a_spatial-blocks_1000m.R")


#==============================================================================#
####                            02b Collinearity                            ####
#==============================================================================#

### run collinearity analysis ###
# description: selects non-correlated variables (Spearman rho < 0.7), based on
# univariate explained deviance, and adapts predictor variable selection to 
# the occurrence data (1 predictor variable per 10 single-events)
# input: envdata (*.tif); occ_env (01b_spp_enviro.rds); spatial-blocks 
# (02a_spp_spatial-blocks.rds)
# output: pred_sel (spp_predictor-vars.rds); pa_env (02b_spp.rds)
# plots: spp_corrplot.png
# notes: split into 3 different grains of analysis (100m, 200m & 1000m)

# Collinearity at 100m
source("scripts/02b_collinearity_100m.R")

# Collinearity at 200m
source("scripts/02b_collinearity_200m.R")

# Collinearity at 1000m
source("scripts/02b_collinearity_1000m.R")


#==============================================================================#
####                            03a Model training                          ####
#==============================================================================#

### train models ###
# description: trains the models, calculates variable importance and simplifies
# the these by excluding variables with importance lower than chance
# input: pa_env (02b_spp.rds); pred_sel (spp_predictor-vars.rds)
# output: full-model (03a_spp_SDMfull.RData); simplified-model 
# (03a_spp_SDMsimp.RData); full_var_imp (04a_spp_full_varImp.csv)

source("scripts/03a_model-training.R")


#==============================================================================#
####                           03b Model validation                         ####
#==============================================================================#

### validate models ###
# description: validates the models, including the calculation of the best 
# threshold for binary output simplification
# input: pa_env (02b_spp.rds); spatial-blocks (02a_spp_spatial-blocks.rds);
# pred_sel (spp_predictor-vars.rds); full-model (03a_spp_SDMfull.RData);
# simplified-model (03a_spp_SDMsimp.RData)
# output: full-model-perf (03b_spp_full_perf-cv-pred.Rdata); 
# simplified-model-perf (03b_spp_simp_perf-cv-pred.Rdata)

source("scripts/03b_model-validation.R")


#==============================================================================#
####                           03c Model prediction                         ####
#==============================================================================#

### predict models ###
# description: predicts the models into probabilistic and binary outputs
# input: envdata (*.tif); pa_env (02b_spp.rds); pred_sel 
# (spp_predictor-vars.rds); simplified-model (03a_spp_SDMsimp.RData); 
# simplified-model-perf (03b_spp_simp_perf-cv-pred.Rdata); gridsizes 
# (sppgridsizes.csv)
# output: model-pred (03c_spp_pred.RData); model-pred-occ 
# (03c_spp_pred_occ.Rdata)
# plots: model-pred-prob (03c_spp_pred.tif); model-pred-bin 
# (03c_spp_pred_bin.tif)

# make predictions for year 2024
source("scripts/03c_model-prediction.R")

# make predictions for all years (2016-2024)
source("scripts/03c_model-prediction_allyears.R")


#==============================================================================#
####                   04a Variable importance calculation                  ####
#==============================================================================#

### calculate variable importance for the final models ###
# description: calculates variable importance for the simplified model
# input: pa_env (02b_spp.rds); simplified-model (03a_spp_SDMsimp.RData); 
# simplified-model-perf (03b_spp_simp_perf-cv-pred.Rdata)
# output: simp_var_imp (04a_spp_simp_varImp.csv)
# notes: this procedure is also included in the model-train code, where it is 
# used to calculate variable importance for the full model

source("scripts/04a_variable-importance.R")


#==============================================================================#
####                       04b Partial dependency plots                     ####
#==============================================================================#

### plot response curves for the final models ###
# description: creates the partial dependency plots for the simplified model 
# variables
# input: pa_env (02b_spp.rds); simplified-model (03a_spp_SDMsimp.RData); 
# simp_var_imp (04a_spp_simp_varImp.csv)
# plots: resp-curves (04b_spp_resp-curves.png)
# notes: this code doesn't work properly for the categorical variables, and 
# could be improved, by e.g. incorporating the next script

source("scripts/04b_response-curves.R")


#==============================================================================#
####                      04c Response values extraction                    ####
#==============================================================================#

### extract response values for categorical variable "BiotopCat" ###
# description: extracts response values for teh categorical variables
# input: pa_env (02b_spp.rds); simplified-model (03a_spp_SDMsimp.RData)
# output: resp-values (04c_spp_response-values.csv)
# notes: could be incorporated in the previous script to optimise the plotting
# of the categorical variables

source("scripts/04c_response-values.R")


#==============================================================================#
####                      04d Model statistics compilation                  ####
#==============================================================================#

### compile model statistics ###
# compiles model performances and variable importance for both the full and
# simplified models
# input: envdata (*.tif); pa_env (02b_spp.rds); full-model 
# (03a_spp_SDMfull.RData); simplified-model (03a_spp_SDMsimp.RData); 
# full-model-perf (03b_spp_full_perf-cv-pred.Rdata); 
# simplified-model-perf (03b_spp_simp_perf-cv-pred.Rdata); full_var_imp 
# (04a_spp_full_varImp.csv); simp_var_imp (04a_spp_simp_varImp.csv)
# output: model-perf (04d_Model-performances.csv); full-model-vars 
# (04d_Model-variables-full.csv); simp-model-vars (04d_Model-variables-simp.csv)

source("scripts/04d_model-statistics.R")


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
