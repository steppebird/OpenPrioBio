
#         Title: 01a occurrence data preparation - 1000m
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão  (adapted Wiedenroth, 2023)
#         Last update: 22.07.2025
#         Usage: Takes bird species data for specific region, for each year, and 
#         individual species, removes duplicates, (accuracy, breeding season,
#         and spatial duplicates inside same 1000m grid cell), sorting of 
#         pseudo-absences, thinning 2000m (twice the grid cell size), then writes
#         outputs for each species and each year as .rds files
#         Assumptions: Some lines in this code may need changing for use with
#         the specific datasets


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

# output directory
output_dir <- file.path("results/01a_occurrence-prep/1000m")


#==============================================================================#
####                         Load data and prepare                          ####
#==============================================================================#

# Read in the data (in shapefile format)
obsp <- read.csv2(file.path(biodatafolder,
                            "GBIF_NRW_filtered_group3_3035.csv"),sep=",")

# Remove NULL values
obs <- obsp[complete.cases(obsp$species, obsp$month, 
                           obsp$coordinateUncertaintyInMeters), ]

# Convert vector dataset to df to start data extraction
obs <- terra::as.data.frame(obs, row.names=NULL, optional=FALSE, geom=NULL,)

obs <- select(obs, c(-speciesKey, -countryCode, -individualCount, -gbifID, 
                     -family, -taxonRank, -basisOfRecord, -institutionCode, 
                     -stateProvince, -eventDate, -day, -NameDeutsch, -NameWiss,
                     -group))

# Add environmental data to (later) remove occurrences that fall within the same
# spatial cell (duplicates) - uses one predictor variable as reference
nvdata <- 
  terra::rast(file.path(envdatafolder,"coregistration_2016_2024_1000m/dgm_NRW.tif"))

# target coordinate reference system (same as envdata)
projcrs <- crs(envdata)

# tidy up files
rm(obsp) #tidy up

# Separate out the individual species and years
unique_yrs <- unique(obs$year)
unique_yrs <- sort(unique_yrs)
unique_spp <- unique(obs$species)
unique_spp <- sort(unique_spp)


#==============================================================================#
####                           Occurrence cleaning                          ####
#==============================================================================#

cat("\n")
cat("Occurrence data preparation")

# Create a list to store final occurrence records
final_occurrences <- list()

# create a data frame to store species totals for internal analysis (e.g., 
# how many occurrences are removed during each filtering step)
species_totals_table <- data.frame()

# Loop round occurrence records for each species, within each year
for (j in unique_spp) {
  
  cat("\n")
  cat(paste("#### Selected species: ",j," ####",sep=""))
  cat("\n")
  
  for (i in unique_yrs) {

    if (nrow(obs[obs$year == i & obs$species == j,]) != 0){
      
      cat("\n")
      cat(paste("## Selected year: ",i," ##",sep=""))
      cat("\n")
      cat("    Retrieving presences")
      cat("\n")

# filter by year and species ---------------------------------------------------

      cat("    Filtering by year and species")
      cat("\n")
    
    # Create the first subset of original which separates the species and year/s
      subset_yr_spp <- obs[obs$year == i & obs$species == j,]
    

# filter by accuracy -----------------------------------------------------------
    
      cat("    Filtering by spatial accuracy")
      cat("\n")
    
      sub1_accuracy <- 
        subset_yr_spp[subset_yr_spp$coordinateUncertaintyInMeters <= 1000, ]
    

# filter by breeding season ----------------------------------------------------

      cat("    Filtering by breeding season")
      cat("\n")
    
      sub3_breed <- 
        sub1_accuracy[sub1_accuracy$month >= 
                        sub1_accuracy$Brut_beginn | sub1_accuracy$month <=
                        sub1_accuracy$Brut_ende,]

    
# project to working coordinate system -----------------------------------------

      cat("    Reprojecting to working coordinate system")
      cat("\n")
    
      sub3_breed[,7] <- as.numeric(sub3_breed[,7])
      sub3_breed[,8] <- as.numeric(sub3_breed[,8])
    
      # Vectorise the last derived subset using coordinates
      sub5_vect <- terra::vect(sub3_breed, geom=c("x", "y"), crs="EPSG:3035")
    
      # Convert the occurrence to working coordinate system
      occ_epsg <- terra::project(sub5_vect, projcrs)
    
    
# remove duplicates ------------------------------------------------------------
    
      cat("    Removing duplicates")
      cat("\n")
    
      # remove duplicate occurrence records to a single occurrence per grid cell
    
      # find duplicate occurrences within each grid cell
      # extract environmental data, coordinates, and the information
      overlay <- cbind(values(occ_epsg), crds(occ_epsg),
                       terra::extract(envdata, occ_epsg, cells = TRUE)) 
    
      # keep only the first occurrence within each cell of the grid
      Nodup <- overlay[!duplicated(overlay$cell),] 
    
      # add an "occurrence" column before combining with absences
      Nodup$occurrence <- 1
    

# create and add pseudo-absences------------------------------------------------
# absences are created so that they don't fall within cells with presences
    
      cat("    Randomly selecting background data / pseudo-absences")
      cat("\n")
    
      # make a regional mask that contains NAs in presence locations
      species_cell <- 
        terra::extract(envdata, Nodup[, c("x", "y")], cells = TRUE)$cell
      presence_mask <- envdata
      values(presence_mask)[species_cell] <- NA
    
      # Randomly select background data but excluding presence locations
      pseudo_abs <- terra::spatSample(presence_mask, size = nrow(Nodup)*10, 
                                      method = "random", na.rm=TRUE, 
                                      as.points=TRUE, cells = TRUE)
    
      # add coordinates and make it a data frame
      pseudo_abs_df <- cbind(values(pseudo_abs), crds(pseudo_abs))
    
      # prepare pseudo-absences to join with presences
      pseudo_abs_df$species <- Nodup$species
      pseudo_abs_df$year <- Nodup$year
      pseudo_abs_df$ID <- NA
      pseudo_abs_df$occurrence <- 0
    
      # remove unnecessary fields from Nodup
      Nodup <- Nodup[,-c(2,4:6)]
    
      # combine presences and absences
      pa <- rbind(Nodup, pseudo_abs_df)
    
      # remove duplicated absences
      pa_nodups <- pa[!duplicated(pa$cell),]    
    
        
# spatial thinning -------------------------------------------------------------

      cat("    Spatial thinning combined presences and pseudo-absences")
      cat("\n")
    
      # transform to sf object to apply thin function
      Nodup_sf <- st_as_sf(pa_nodups, coords = c("x", "y"), crs=projcrs)
    
      # apply thinning function (adapt distance to 2x grid size)
      thinned <- thin(Nodup_sf, thin_dist = 2000, runs = 5, ncores = 5)
    
      # update column name in thinned file to match x/y of output
      colnames(thinned) <- c("x", "y")
    
      # join the original columns back to the thinned rows
      final_subset <- inner_join(pa_nodups, thinned, by = c("x", "y"))
    

# saving the outputs -----------------------------------------------------------

      cat("    Compiling and saving the outputs")
      cat("\n")
    
      # add the final subset to store in the final outputs list: one output per 
      # species and year
      final_occurrences[[paste0("01a_",j, i)]] <- final_subset
    
      # Populate results of species totals table (species by species)
      species_totals_table <- rbind(species_totals_table, data.frame(
        Species = j,
        Year = i,
        Initial_presence = nrow(subset_yr_spp),
        accuracy = nrow(sub1_accuracy),
        Breeding_presence = nrow(sub3_breed),
        NoDups = nrow(Nodup),
        Absence = nrow(pseudo_abs_df),
        Combined_PA = nrow(pa_nodups),
        Thinned_occurrences = nrow(thinned),
        Final_pres = sum(final_subset$occurrence == 1),
        Final_abs = sum(final_subset$occurrence == 0)))
    }
  }
}


#==============================================================================#
####        Write combined RDS with all species for specific year           ####
#==============================================================================#

# Write prepared species df with geometry included
lapply(names(final_occurrences), function(name) {
  saveRDS(final_occurrences[[name]], 
          file = file.path(output_dir, paste0(name, ".rds")))
})

# save combined list (all species together)
# save a .rds file with extension .rd, to not conflict with script 01b
saveRDS(final_occurrences, 
        file = file.path(output_dir, "final_occurrences_combined.rd"))

# write species totals table to file
write.csv(species_totals_table, 
          file = file.path(output_dir, "species_results_table.csv"), 
          row.names = FALSE)


#==============================================================================#
####                           End of work flow                             ####
#==============================================================================#
