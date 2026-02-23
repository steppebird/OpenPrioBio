
#         Title: 00a_folder-structure
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão (adapted from Wiedenroth, 2023)
#         Last update: 26.05.2025
#         Usage: Sets up the necessary folder structure for running the SDM
#         models within a specific project
#         Assumptions: The main (analysis) project folder is the working 
#         directory; this script is located in the scripts folder
#         Notes: In case the analysis needs to be split into several blocks/
#         guilds, this structure might need to be adapted accordingly


#==============================================================================#
####                             Folder structure                           ####
#==============================================================================#

folders <- c(
  # data folders
  "data/raw",
  "data/processed/envdata",
  "data/processed/biodata",
  
  # results folder
  "results/01a_occurrence-prep/100m",
  "results/01a_occurrence-prep/200m",
  "results/01a_occurrence-prep/1000m",
  "results/01b_enviro-extraction/100m",
  "results/01b_enviro-extraction/200m",
  "results/01b_enviro-extraction/1000m",
  "results/02a_spatial-blocks/plots",
  "results/02b_collinearity/corrplots",
  "results/02b_collinearity/env-occ-df",
  "results/03a_model-training",
  "results/03b_model-validation",
  "results/03c_model-prediction",
  "results/04a_variable-importance",
  "results/04b_response-curves",
  "results/04c_response-values",
  "results/04d_model-statistics",
  
  # script folders
  "scripts/"
)

# Create each folder if it does not exist yet
for (folder in folders) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
}


#==============================================================================#
####                            End of work flow                            ####
#==============================================================================#

