
#         Title: 02b collinearity workflow - 1000m
#         Project: OpenPrioBio
#         Author: Pedro J. Leitão (adapted Wiedenroth, 2023)
#         Last update: 22.07.2025
#         Usage: Selects non correlated variables for species distribution 
#         models. Not correlated variables have a (Spearman) correlation 
#         coefficient of < 0.7. If two variables are correlated, the most
#         important variable is kept. Variable importance is determined via the 
#         explained deviance of univariate models.


#==============================================================================#
####                            Workspace set up                            ####
#==============================================================================#

# set working directory
setwd("/home/eouser/Projects/OpenPrioBio/Analysis_NRW2/")

# load packages and functions
source("scripts/00b_functions.R")

# data folder
biodatafolder <- file.path("data/processed/biodata")

# results folder
resultsfolder <- file.path("results")

# list occ-env-df
occ_env_files <- list.files(file.path(resultsfolder, "01b_enviro-extraction/1000m"),
                           full.names = TRUE, pattern = ".rds")

# files spatial blocks
spatial_blocks_files <- list.files(file.path(resultsfolder, "02a_spatial-blocks"), 
                                   full.names = TRUE, pattern = ".rds")

# set output_dir to 02b
output_dir <- file.path(resultsfolder, "02b_collinearity")

# loop over all species and test for collinearity between the predictor variables
for (i in occ_env_files) {
  
#==============================================================================#
####                                Load data                               ####
#==============================================================================#
  
  # load load occ-env-df
  occ <- readRDS(i)
  
  cat("\n")
  cat("Retrieving files...")
  cat("\n")
  
  # Extract species names
  spp <- occ$species[1]
  
  cat(paste("Selected species: ",spp,sep=""))
  cat("\n")
  
  # load spatial blocks
  scv <- readRDS(spatial_blocks_files[grep(spp, spatial_blocks_files)])
  
  
#==============================================================================#
####                              Collinearity                              ####
#==============================================================================#
  
  # set weights for the presences as 1. The weights of the pseudo absences are
  # calculated so that their sum equals the sum of the weights of the presences
  wgt <- ifelse(occ$occurrence == 1, 1, 
                sum(occ$occurrence == 1) / sum(occ$occurrence == 0 ))
  
  # remove variables with a stdev of zero (excl. categorical variable)
  notcat <- !(envpreds %in% c("forest_type","lulc","stocked_forest","tree_spec"))
  stdev_values <- sapply(occ[, envpreds[notcat]], function(x) {
    sd(x, na.rm = TRUE)
  })
  zero_stdev_vars <- names(stdev_values)[stdev_values == 0]
  
  envpreds <- envpreds[!envpreds %in% zero_stdev_vars]
  
  if(length(zero_stdev_vars)!=0){
    
    cat("    The following variable had a stdev of zero and were removed:")
    cat("\n")
    cat(paste0("     > ",zero_stdev_vars,"\n"))
    cat("\n")
  }
  
  cat("    Calculating pairwise Spearman correlations")
  cat("\n")
  
  # check for collinearity between predictors (excl. categorical variable)
  notcat <- !(envpreds %in% c("forest_type","lulc","stocked_forest","tree_spec"))
  cor_mat <- cor(occ[, envpreds[notcat]], method = "spearman")
  
  # save collinearity plot per species
  # open png device for writing to file
  png(filename = paste0(output_dir, "/corrplots/", spp, "_corrplot.png"),
      width = 1600, height = 1600, res = 240)
  
  # plot the correlation matrix
  corrplot::corrplot.mixed(cor_mat, tl.pos = "lt", tl.cex = 0.6,
                           number.cex = 0.5, addCoefasPercent = TRUE)
  
  # close device for writing to file
  dev.off()
  
  cat("     > Correlation plot has been saved in corrplots folder")
  cat("\n")
  cat("    Removing predictors with too few non-zero values:")
  cat("\n")
  
  ## Identify predictors with too few non-zero values ----
  ## and remove from initial predictor list
  pred_meta <- occ[,envpreds[notcat]]  %>%
    tidyr::pivot_longer(everything(), names_to = "predictor", 
                        values_to = "value")%>%
    group_by(predictor) %>%
    reframe(perc.zeros = length(which(value == 0)) / n(),
            n.not.zero = length(which(value != 0))) %>%
    arrange(desc(perc.zeros))
  
  ### keep only predictors with more than 1% non-zero values   
  keep.pred.thresh <- .99
  pred.to.remove <- 
    pred_meta$predictor[pred_meta$perc.zeros >= keep.pred.thresh]   
  
  if(length(pred.to.remove)>0){
    
    cat(paste0("     > Predictor variables removed for ", spp, ":"))
    cat("\n")
    cat(paste0("     > ",pred.to.remove,"\n"))
    cat("\n")
  }   else {
    
    cat("     > No predictor variables removed")
    cat("\n")
  }
  
  pred.list = envpreds
  pred.list = pred.list[!pred.list %in% pred.to.remove]
  
  cat("    Selecting non-correlated variables:")
  cat("\n")
  
  # running function to determine non-correlated variables
  notcat <- !(pred.list %in% c("forest_type","lulc","stocked_forest","tree_spec"))
  var.sel <- select07_blockcv(X = occ[, pred.list[notcat]],
                              y = occ$occurrence,
                              threshold = 0.7,
                              univar = NULL,
                              sp_block = scv,
                              weights = wgt)
  
  num.vars <- length(var.sel$pred_sel)+4
  
  cat(paste0("     > ",num.vars, " variables selected"))
  cat("\n")
  cat("    Adapting variable selection to occurence data:")
  cat("\n")
  
  # calculate number of variables to be selected
  # 1 variable for each 10 presences or absences
  occ_num1 <- sum(occ$occurrence)/10
  occ_num2 <- sum(occ$occurrence==0)/10
  
  # select lowest number, between presences and absences 
  occ_num <- if(occ_num1 < occ_num2) {
    
    occ_num1
  } else {
    
    occ_num2
  }
  
  # identify the variables that are not correlated. 
  pred.sel <- var.sel$pred_sel
  
  # select the number of predictor variables based on the number of occurrences
  pred.sel <- pred.sel[1:(occ_num)]
  pred.sel <- na.omit(pred.sel)
  
  # include categorical variables
  pred.sel <- c(pred.sel,"forest_type","lulc","stocked_forest","tree_spec")
  
  num.selvars <- length(pred.sel)
  
  cat(paste0("     > ",num.selvars, " variables remaining"))
  cat("\n")
  cat("    List of selected variables:")
  cat("\n")
  cat(paste0("     > ",pred.sel,"\n"))
  
  # save selected variable names into a single file
  saveRDS(pred.sel, file = file.path(output_dir,
                                     paste0("02b_",spp,"_pred-vars.rds")))
  
  # save env-occ-df but remove the occurrence locations due to the sensitivity
  # of this information
  occ <- dplyr::select(occ, species, occurrence, all_of(pred.sel))
  saveRDS(occ, file = file.path(output_dir,
                                paste0("env-occ-df/02b_",spp,".rds")))
  
} # end collinearity loop around each species


#==============================================================================#
####                             End of workflow                            ####
#==============================================================================#
