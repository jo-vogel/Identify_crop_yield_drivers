# Fitting of the Lasso model
# Authors: Pauline Rivoire, Johannes Vogel, Cristina Deidda

# This code was run under version 3.6 of R

# Structure of the code
# a) Load the standardized data (meteorological variables and extremal indicators)
# b) Create training and testing dataset
# c) Run the cross validation to find lambda.min and lambda.1se
# d) Run the Lasso regression 


# Required libraries ####
library(glmnet);library(InformationValue);library(abind);library(stringr);library(ggplot2)

# Load standardized Data #####
path_data <- message("insert data directory here") # Insert path of the input data here
path_model <- message("insert model directory here") # where to store the output model
path_code <- message("insert code directory here") # Insert path of the code here

load(paste0(path_data,"extremeindices_and_monthlymeteovar_rescaled_995pix.RData"))

source(paste0(path_code,"/Data_processing.R"))

var_num <- apply(non_na_col, 1, sum)
numLevels_list <- sapply(1:pix_num, function(x){ rep(1, times = var_num[x])})
for (i in 1:pix_num){
  names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
}



# Run the CrossValidation #####
# With the cross validation, we obtain lambda.min and lambda.1se
model_cv_fitting <- list()
nbyears_final_training_data <- numeric()

for (pixel in 1:pix_num) {
  
  var_pix <- as.matrix(x1_train_list[[pixel]])
  yield_pix <- as.matrix(y1_train_list[[pixel]])
  which_na_xtrm <- which(is.na(var_pix[,1]))
  nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
  
  if (sum(yield_pix[-which_na_xtrm,])<=nbyears_final_training_data[pixel] &
      sum(yield_pix[-which_na_xtrm,])>=nbyears_final_training_data[pixel]-8){
    #For a correct fitting, we need at least 8 occurences of good years and bad years
    #If it is not the case, the model is not run
    model_cv_fitting[[pixel]] <- "Training years w/o na in extremes have less than 8 bad years"
  } else {
    if(length(which_na_xtrm)>0){ #remove potential missing data in extreme indices
      model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix[-which_na_xtrm,],
                                             y = yield_pix[-which_na_xtrm,],
                                             family = "binomial",
                                             alpha = 1, #alpha = 1 corresponds to the Lasso regression
                                             nfolds = 10) #number of folds for the cross validation
    } else {
      if(sum(yield_pix)>=nbyears_final_training_data[pixel]-8){
        #For a correct fitting, we need at least 8 occurences of good years and bad years
        #If it is not the case, the model is not run
        model_cv_fitting[[pixel]] <- "Training years have less than 8 bad years"
      } else {
        model_cv_fitting[[pixel]] <- cv.glmnet(x = var_pix,
                                               y = yield_pix,
                                               family = "binomial",
                                               alpha = 1,
                                               nfolds = 10)
      }#end if else not enough bad occurences
      
    }#end if exists na else
    
  }#end if all years kept are good else
  
  print(paste(pixel, "out of", pix_num))
}#end for pixel


save(model_cv_fitting, file = paste0(path_model, "/cv_month_xtrm_LASSO_threshbadyield005_seed",
                                     seed, "_train", train_size,"_995pixels.RData"))


# Run the model with lambda1se and lambda min just obtained #####
load(file = paste0(path_model, "cv_month_xtrm_LASSO_threshbadyield005_seed",
                   seed, "_train", train_size, "_995pixels.RData"))


lasso_model_lambdamin <- list()
lasso_model_lambda1se <- list()

for (pixel in 1:pix_num) {
  if(is.character(model_cv_fitting[[pixel]])){ #If the cross validation could not be run
    lasso_model_lambdamin[[pixel]] <- model_cv_fitting[[pixel]]
    lasso_model_lambda1se[[pixel]] <- model_cv_fitting[[pixel]]
  } else {
    var_pix <- as.matrix(x1_train_list[[pixel]])
    yield_pix <- as.matrix(y1_train_list[[pixel]])
    which_na_xtrm <- which(is.na(var_pix[,1]))
    nbyears_final_training_data[pixel] <- (dim(var_pix)[1]-length(which_na_xtrm))
    
    training_years_wo_na <- which(!is.na(x1_train_list[[pixel]]$dtr))
    
    if(length(which_na_xtrm)>0){
      
      lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
                                               y = yield_pix[-which_na_xtrm,],
                                               family = "binomial",
                                               alpha = 1,
                                               lambda = model_cv_fitting[[pixel]]$lambda.min)
      
      lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix[-which_na_xtrm,],
                                               y = yield_pix[-which_na_xtrm,],
                                               family = "binomial",
                                               alpha = 1,
                                               lambda = model_cv_fitting[[pixel]]$lambda.1se)
    } else {
      
      lasso_model_lambdamin[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
                                               family = "binomial",
                                               alpha = 1,
                                               lambda = model_cv_fitting[[pixel]]$lambda.min)
      
      lasso_model_lambda1se[[pixel]] <- glmnet(x = var_pix, y = yield_pix,
                                               family = "binomial",
                                               alpha = 1,
                                               lambda = model_cv_fitting[[pixel]]$lambda.1se)
    }#end if else there exists na in extremes
    
  } #end ifelse the cross validation could not be run
  print(paste(pixel, "out of", pix_num))
  
}#end for pixel

save(lasso_model_lambdamin, file = paste0(path_model, "/Lasso_lambdamin_month_xtrm_LASSO_threshbadyield005_seed",
                                          seed, "_train", train_size,"_995pix.RData"))
save(lasso_model_lambda1se, file = paste0(path_model, "/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",
                                          seed, "_train", train_size,"_995pix.RData"))