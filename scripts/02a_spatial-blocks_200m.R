
#         Title: 02a spatial blocks workflow - 200m
#         Project: OpenPrioBio
#         Author:Pedro J. Leitão (adapted Wiedenroth, 2023)
#         Last update: 22.07.2025
#         Usage: Uses combined species and environmental data (outputs from 02a)
#         and produces spatial blocks for cross-validation purposes
#         Assumptions: 5 folds, at least 15 blocks per fold


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

# set output_dir to 02a
output_dir <- file.path(resultsfolder, "02a_spatial-blocks/")

# list environmental data
envdata_static_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_200m/static"), 
             full.names = TRUE, pattern = ".tif")

envdata_dynamic_files <- 
  list.files(file.path(envdatafolder, "coregistration_2016_2024_200m/dynamic"), 
             full.names = TRUE, pattern = ".tif")

# list occ-env-df
pa_env_files <- 
  list.files(file.path(resultsfolder, "01b_enviro-extraction/200m"),
                           full.names = TRUE, pattern = ".rds")


#==============================================================================#
####                                Load data                               ####
#==============================================================================#

# loop over all species
for (i in pa_env_files) {
  
  # load occ-env-df
  occ <- readRDS(i)
  
  cat("\n")
  cat("Retrieving file...")
  cat("\n")
  
  # Extract species names
  spp <- occ$species[1] # get individual file name for spp 
  
  cat(paste0("Selected species: ",spp))
  cat("\n")
  cat("Selected year for spatial blocks: 2024")
  cat("\n")
  
  # load environmental data
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
  
  # select latest dynamic environmental layers (2024)
  envdyn_index <- grep("2024", names(envdata_dynamic))
  envdyn_select <- subset(envdata_dynamic, envdyn_index)
  
  # rename variables
  current_names <- names(envdyn_select)
  new_names <- sub("_2024","", current_names)
  names(envdyn_select) <- new_names
  names(envdyn_select)[16] <- "soil_moist"
  rm(current_names, new_names)
  
  envdata <- c(envdata_static,envdyn_select)
  
  # fix categorical variable
  envdata_static <- envdata_cat(envdata_static, 9)
  envdata_static <- envdata_cat(envdata_static, 11)
  envdata_static <- envdata_cat(envdata_static, 21)
  envdata_static <- envdata_cat(envdata_static, 24)
  

#==============================================================================#
####                             Spatial blocks                             ####
#==============================================================================#
  
  # make a sf object out of the dataframe in order to determine the block_size
  # via the cv_spatial_autocor function
  sf_spec_occ <- sf::st_as_sf(occ, coords = c("x", "y"), crs = crs(envdata))
  
  cat("    Calculating spatial block size")
  cat("\n")
  
  # determine block size 
  # determine the block size based on spatial autocorrelation (minimum 15)
  cv_blocksize <- cv_spatial_autocor(x = sf_spec_occ, r = envdata,
                                     column = "occurrence", plot = FALSE)
  
  ### plot with species in title
  # open png device for writing to file
  png(filename = paste0(output_dir, "plots/" , spp, "_blocksize.png"),
      width = 2000, height = 2000, res = 240)
  
  # plot block size map
  print(cv_blocksize$plots + ggplot2::labs(title = occ$spp[1]))
  
  # close device for writing to file
  dev.off()
  
  cat("     > Block size map saved in plots folder")
  cat("\n")
  
  if(cv_blocksize$range >= 70000){                          # adapt as necessary
    cv_bl_size <- 5000
    
    cat("     > Block size too large: setting block size to 5000m")
    cat("\n")
    
  } else{
    cv_bl_size <- cv_blocksize$range
    
    cat(paste0("     > Block size: ", round(cv_bl_size, digits = 2), "m"))
    cat("\n")
  }
  
  cat("    Creating spatial blocks")
  cat("\n")
  
  # create spatial blocks
  scv <- cv_spatial(
    x = sf_spec_occ,
    column = "occurrence",
    r = envdata,
    k = 5, # number of folds
    size = cv_bl_size, # size of blocks in meters
    selection = "random", # random blocks-to-fold
    iteration = 50, # find evenly dispersed folds
    progress = FALSE, # turn off progress bar
    biomod2 = FALSE, # do not create folds for biomod2
    raster_colors = terrain.colors(10, rev = TRUE),
    plot = F, # do not plot by default
    report = F # do not report the sizes of the subsets 
  )
  
  # number of block < 15
  # automate the block size selection to have at least 15 blocks. If at first
  # we had less than 15 blocks then we should end up with exactly 15 blocks.
  blocknr <- length(scv$blocks$block_id) < 15
  
  if(blocknr){
    
    cat("    > Too few blocks created. Adapting block size to get 15 blocks")
    cat("\n")
  } else{
    
    ### plot with spp in title
    # open png device for writing to file
    png(filename = paste0(output_dir, "/plots/" , spp, "_spatialblocks.png"),
        width = 3000, height = 3000, res = 240)
    
    print(cv_plot(scv, r = envdata,
                  raster_colors = terrain.colors(10,rev = TRUE),
                  label_size = 3)+ ggplot2::labs(title = occ$spp[1]))
    
    # close device for writing to file
    dev.off()
  }
  
  cat("     > Spatial blocks map saved in plots folder")
  cat("\n")
  
  # create exactly 15 blocks
  while(blocknr) {
    
    if(length(scv$blocks$block_id) < 15) {
      
      cv_bl_size <- cv_bl_size/1.1
      
      scv <- cv_spatial(
        x = sf_spec_occ,
        column = "occurrence",
        r = envdata,
        k = 5, # number of folds
        size = cv_bl_size, # size of blocks in meters
        selection = "random", # random blocks-to-fold
        iteration = 50, # find evenly dispersed folds
        progress = FALSE, # turn off progress bar
        biomod2 = FALSE, # do not create folds for biomod2
        raster_colors = terrain.colors(10, rev = TRUE), # options from cv_plot for a better color contrast
        report = FALSE, # do not report the sizes of the subsets
        plot = FALSE) # do not plot the results
    } else
      if(length(scv$blocks$block_id) == 15) {
        
        # open png device for writing to file
        png(filename = paste0(output_dir, "/plots/" ,
                              spp, "_spatialblocks.png"),
            width = 3000, height = 3000, res = 240)
        
        print(cv_plot(scv, r = envdata,
                      raster_colors = terrain.colors(10,rev = TRUE),
                      label_size = 3)+ ggplot2::labs(title = occ$spp[1]))
        
        # close device for writing to file
        dev.off()
        
        cat("     > Spatial blocks map saved in plots folder")
        cat("\n")
        
        blocknr <- FALSE
      } else
        if(length(scv$blocks$block_id) > 15) {
          
          cv_bl_size <- cv_bl_size * 1.01
          
          scv <- cv_spatial(
            x = sf_spec_occ,
            column = "occurrence",
            r = envdata,
            k = 5, # number of folds
            size = cv_bl_size, # size of blocks in meters
            selection = "random", # random blocks-to-fold
            iteration = 50, # find evenly dispersed folds
            progress = FALSE, # turn off progress bar
            biomod2 = FALSE, # do not create folds for biomod2
            raster_colors = terrain.colors(10, rev = TRUE), # options from cv_plot for a better color contrast
            report = FALSE, # do not report the sizes of the subsets
            plot = FALSE) # do not plot the results
        } else
          if(length(scv$blocks$block_id) == 15) {
            
            # open png device for writing to file
            png(filename = paste0(output_dir, "/plots/" ,
                                  spp, "_spatialblocks.png"),
                width = 3000, height = 3000, res = 240)
            
            print(cv_plot(scv, r = envdata,
                          raster_colors = terrain.colors(10,rev = TRUE),
                          label_size = 3)+ ggplot2::labs(title = occ$spp[1]))
            
            # close device for writing to file
            dev.off()
            
            cat("     > Spatial blocks map saved in plots folder")
            cat("\n")
            
            blocknr <- FALSE
          }
  }
  
  cat(paste0("    Created ", length(scv$blocks$block_id), " blocks, with block size: ",
             round(cv_bl_size, digits = 2),"m"))
  cat("\n")
  
  # save spatial blocks to file
  saveRDS(scv, file = file.path(output_dir, paste0("02a_", spp,"_spatial-blocks.rds")))
  
} # end of spatial blocks per spp


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
