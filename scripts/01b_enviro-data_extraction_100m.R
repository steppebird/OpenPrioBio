
#         Title: 01b environmental data extraction - 100m
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão  (adapted Wiedenroth, 2023)
#         Last update: 22.07.2025
#         Usage: Takes processed (01a) spp data for a specific region, 
#         extracts environmental data layers for each location, and writes 
#         outputs for each spp as .RDS is (maintaining naming)
#         Assumptions: The environmental data is not currently attached to the 
#         occurrences. The species-specific variables should be available as a 
#         ".csv" in the "envdata" folder. The "final_occurrences_combined.rds"
#         file should not be in the 01a_occurrence-prep folder.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# load packages and functions
source("scripts/00b_functions.R")

# data folder
envdatafolder <- 
  file.path("/new_nfs_shared/data/04_PROCESSED_DATA/coregistration/NRW")
biodatafolder <- file.path("data/processed/biodata")

# results folder
resultsfolder <- file.path("results/")

# set output_dir to 01b
output_dir <- file.path(resultsfolder, "01b_enviro-extraction/100m")

# list the occurrence dataframes
occurrence_list <- list.files(file.path(resultsfolder,
                                        "01a_occurrence-prep/100m"), 
                              full.names = TRUE, pattern = ".rds") 

# load in the environmental data (both static and dynamic variables)
cat("\n")
cat("Retrieving environmental data...")
cat("\n")

envdata_static_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_100m/static"), 
                            full.names = TRUE, pattern = ".tif")

envdata_dynamic_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_100m/dynamic"), 
                                   full.names = TRUE, pattern = ".tif")

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


#==============================================================================#
####                    Environmental data extraction                       ####
#==============================================================================#

# create a list to store extracted species/environment occurrence records
enviro_occurrences <- list()

for (i in occurrence_list) {

  cat("\n")
  cat("Retrieving file...")
  cat("\n")
  
  # read in the .rds from 01a (presences and absences of the respective species)
  occ <- readRDS(i)
  
  # get species name
  Species <- unique(occ$species)
  
  # remove environmental columns
  occ <- occ[,-c(5:7)]

  # select (year specific) dynamic environmental variables
  selected_year <- unique(occ$year)
  
  cat(paste("Selected species: ",Species,sep=""))
  cat("\n")
  cat(paste("Selected year: ",selected_year,sep=""))
  cat("\n")
  
  envdyn_index <- grep(selected_year, names(envdata_dynamic))
  envdyn_select <- subset(envdata_dynamic, envdyn_index)
  
  # rename variables
  current_names <- names(envdyn_select)
  new_names <- sub(paste0("_",selected_year),"", current_names)
  names(envdyn_select) <- new_names
  names(envdyn_select)[16] <- "soil_moist"
  rm(current_names, new_names)
  
  cat("    Extracting environmental data for each occurrence")
  cat("\n")
  
  # extract predictor data for each occurrence cell and add it to data frame
  occ_env_year <- cbind(occ, terra::extract(x = envdata_static,
                                            y = occ[,c('x','y')], 
                                       cells=F)[,-1],
                   terra::extract(x = envdyn_select, y = occ[,c('x','y')], 
                                  cells=T)[,-1])
  
  # only retain rows that do not contain NAs in columns with env. data
  occ_noNAs_year <- subset(occ_env_year, 
                      rowSums(is.na(occ_env_year)) == 0)
  
  cat("    Compiling data")
  cat("\n")
  
  # store processed enviro-occs in "enviro_occurrences" list
  enviro_occurrences[[paste0("01b_", Species, selected_year,"_enviro")]] <- 
    occ_noNAs_year
  
} # end of loop over species in occurrence_list

# compile single species occurrences list (all years)

enviro_occurrences_spp <- list()

spp_allocc <- unique(sapply(enviro_occurrences, '[[', "species"))
spp_allocc <- unlist(spp_allocc)
unique_spp <- unique(spp_allocc)
  
for (spp in unique_spp) {
  
  cat("\n")
  cat(paste("Compiling multi-year occurrences for ",spp,sep=""))
  cat("\n")
  
  selspp_allocc <- enviro_occurrences[grep(spp, names(enviro_occurrences))]
  allocc <- as.data.frame(selspp_allocc[1][[1]])
  
  if(length(selspp_allocc)>1){
    for (o in 2:length(selspp_allocc)){
      allocc <- allocc %>% add_row(as.data.frame(selspp_allocc[o][[1]]))
    }
  }
  
  enviro_occurrences_spp[[paste0("01b_",spp,"_enviro")]] <- allocc
}


#==============================================================================#
####          Write a combined RDS file by species and year 2022            ####
#==============================================================================#

# write prepared species df to file
lapply(names(enviro_occurrences_spp), function(name) {
  saveRDS(enviro_occurrences_spp[[name]], file = file.path(output_dir, 
                                                       paste0(name,".rds")))
})


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
