# This file contains functions required to run "Figures.R"
# Authors: Pauline Rivoire, Johannes Vogel

# This code was run under version 3.6 of R

# Structure of the code
# a) Functions for the cut off level
# b) Function to get the name of coefficients sorted by decresing absolute value
# c) Additional functions, useful to plot Fig 7.




# Functions for the cut off level ####

library(ROCR);library(grid);library(caret);library(dplyr);library(scales)
library(ggplot2);library(gridExtra);library(data.table)

# ------------------------------------------------------------------------------------------
# [AccuracyCutoffInfo] : 
# Obtain the accuracy on the trainining and testing dataset.
# for cutoff value ranging from .4 to .8 ( with a .05 increase )
# @train   : your data.table or data.frame type training data ( assumes you have the predicted score in it ).
# @test    : your data.table or data.frame type testing data
# @predict : prediction's column name (assumes the same for training and testing set)
# @actual  : actual results' column name
# returns  : 1. data : a data.table with three columns.
#            		   each row indicates the cutoff value and the accuracy for the 
#            		   train and test set respectively.
# 			 2. plot : plot that visualizes the data.table

AccuracyCutoffInfo <- function( train, test, predict, actual ){
  # change the cutoff value's range as you please 
  cutoff <- seq( .4, .8, by = .05 )
  
  accuracy <- lapply( cutoff, function(c)
  {
    # use the confusionMatrix from the caret package
    cm_train <- confusionMatrix( as.factor(as.numeric( train[[predict]] > c )), as.factor(train[[actual]]) )
    cm_test  <- confusionMatrix( as.factor(as.numeric( test[[predict]]  > c )), as.factor(test[[actual]])  )
    
    dt <- data.table( cutoff = c,
                      train  = cm_train$overall[["Accuracy"]],
                      test   = cm_test$overall[["Accuracy"]] )
    return(dt)
  }) %>% rbindlist()
  
  # visualize the accuracy of the train and test set for different cutoff value 
  # accuracy in percentage.
  accuracy_long <- gather( accuracy, "data", "accuracy", -1 )
  
  plot <- ggplot( accuracy_long, aes( cutoff, accuracy, group = data, color = data ) ) + 
    geom_line( size = 1 ) + geom_point( size = 3 ) +
    scale_y_continuous( label = percent ) +
    ggtitle( "Train/Test Accuracy for Different Cutoff" )
  
  return( list( data = accuracy, plot = plot ) )
} #end function AccuracyCutoffInfo


# ------------------------------------------------------------------------------------------
# [ConfusionMatrixInfo] : 
# Obtain the confusion matrix plot and data.table for a given
# dataset that already consists the predicted score and actual outcome.
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome 
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cutoff  : cutoff value for the prediction score 
# return   : 1. data : a data.table consisting of three columns
#            		   the first two stores the original value of the prediction and actual outcome from
#			 		   the passed in data frame, the third indicates the type, which is after choosing the 
#			 		   cutoff value, will this row be a true/false positive/ negative 
#            2. plot : plot that visualizes the data.table 

ConfusionMatrixInfo <- function( data, predict, actual, cutoff )
{	
  # extract the column ;
  # relevel making 1 appears on the more commonly seen position in 
  # a two by two confusion matrix	
  predict <- data[[predict]]
  actual  <- relevel( as.factor( data[[actual]] ), "1" )
  
  result <- data.table( actual = actual, predict = predict )
  
  # caculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 1, "TP",
                            ifelse( predict >= cutoff & actual == 0, "FP", 
                                    ifelse( predict <  cutoff & actual == 1, "FN", "TN" ) ) ) %>% as.factor() ]
  
  # jittering : can spread the points along the x axis 
  plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
    geom_violin( fill = "white", color = NA ) +
    geom_jitter( shape = 1 ) + 
    geom_hline( yintercept = cutoff, color = "blue", alpha = 0.6 ) + 
    scale_y_continuous( limits = c( 0, 1 ) ) + 
    scale_color_discrete( breaks = c( "TP", "FN", "FP", "TN" ) ) + # ordering of the legend 
    guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
    ggtitle( sprintf( "Confusion Matrix with Cutoff at %.2f", cutoff ) )
  
  return( list( data = result, plot = plot ) )
}


# ------------------------------------------------------------------------------------------
# [ROCInfo] : 
# Pass in the data that already consists the predicted score and actual outcome.
# to obtain the ROC curve 
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cost.fp : associated cost for a false positive 
# @cost.fn : associated cost for a false negative 
# return   : a list containing  
#			 1. plot        : a side by side roc and cost plot, title showing optimal cutoff value
# 				 	   		  title showing optimal cutoff, total cost, and area under the curve (auc)
# 		     2. cutoff      : optimal cutoff value according to the specified fp/fn cost 
#		     3. totalcost   : total cost according to the specified fp/fn cost
#			 4. auc 		: area under the curve
#		     5. sensitivity : TP / (TP + FN)
#		     6. specificity : TN / (FP + TN)

ROCInfo <- function( data, predict, actual, cost.fp, cost.fn ){
  # calculate the values using the ROCR library
  # true positive, false postive 
  # pred <- prediction( data[[predict]], data[[actual]] )
  pred <- prediction( data[["predict"]], data[["actual"]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost 
  # false postive rate * number of negative instances * false positive cost + 
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 0 ) + 
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 1 )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)   
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) + 
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )				
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    # scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )	
  
  # the main title for the two arranged plot
  sub_title <- sprintf( "Cutoff at %.2f - Total Cost = %.f, AUC = %.3f", best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2,
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ), vjust=1, hjust=-0.1 ) )
  
  return( list( plot 		  = plot, 
                cutoff 	  = best_cutoff, 
                totalcost   = best_cost, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr ) )
} #end function ROCInfo





adjust_cutoff <- function(model_vector, x1_train_list, y1_train_list, work_pix, cost_fp = 100, cost_fn = 100){
  
  # Finding appropriate cutoff level ####
  #######################################
  
  
  y1_train_list_red <- as.list(lapply(work_pix, function(work_pix){y1_train_list[[work_pix]]}))
  mypred_train <- lapply(work_pix, function(x){predict(model_vector[[x]],as.matrix(x1_train_list[[x]]), type = "response")})
  
  # Data set with actuals and predictions
  data_train_all <- pblapply(1:length(work_pix),
                             function(x){ data.frame("Actuals" = as.numeric(y1_train_list_red[[x]]),
                                                     "Predictions" = as.numeric(mypred_train[[x]]))}) # train data
  # Calculate confusion matrix
  cm_info_all <-  pblapply(1:length(work_pix),
                           function(x){ConfusionMatrixInfo( data = data_train_all[[x]], predict = "Predictions", 
                                                            actual = "Actuals", cutoff = .5 )})
  # Calculate ROC curve and cost function
  roc_info_all <- pblapply(1:length(work_pix),
                           function(x){ROCInfo( data = cm_info_all[[x]]$data, predict = "predict",
                                                actual = "actual", cost.fp = cost_fp, cost.fn = cost_fn )})
  # note: the cutoff from cm_info_all has no role here
  cutoff_avg <- pbsapply(1:length(work_pix), function(x){roc_info_all[[x]]$cutoff}) # find the cutoff value
  return(mean(cutoff_avg)) # calculate the average cutoff value
  
} # end  functionadjust_cutoff




# Function to get the name of coefficients sorted by decresing absolute value ####
# Inputs :
#           -coeff:         coefficients from Lasso regression
#           -nb_of_coeff:   number of coefficients wanted (sorted by order of importance)
#
# Outputs :
#           -coeff_name_date: a matrix of dimensions nb_of_coeff*2
#                             column 1 = name of the variable
#                             column 2 = month of the variable (NA if variable is extremal indicator)
# Line 1 corresponds to the variable with the highest absolute value of the fitted coefficient
# Line 2 corresponds to the variable with the second highest absolute value, etc.

get_firstcoeffs <- function(coeff, nb_of_coeff = 1){
  N <- nb_of_coeff
  coeff_name_date <- matrix(data = NA, nrow = N, ncol = 2)
  row.names(coeff_name_date) <- as.character(1:N)
  colnames(coeff_name_date) <- c("var", "month")
  
  for (ind in 1:N) {
    coeff_names <- rownames(coeff)[sort(abs(as.numeric(coeff)), decreasing = T,
                                        index.return=T)$ix[-which(sort(abs(as.numeric(coeff)), 
                                                                       decreasing = T,
                                                                       index.return=T)$ix==1)]][ind]
    if(substr(coeff_names, start = 1, stop = 3)=="dtr"){
      coeff_name_date[ind,1] <- "dtr"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="frs"){
      coeff_name_date[ind,1] <- "frs"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="txx"){
      coeff_name_date[ind,1] <- "txx"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="tnn"){
      coeff_name_date[ind,1] <- "tnn"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 3)=="rx5"){
      coeff_name_date[ind,1] <- "rx5"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 5)=="tx90p"){
      coeff_name_date[ind,1] <- "tx90p"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 5)=="tn10p"){
      coeff_name_date[ind,1] <- "tn10p"
      coeff_name_date[ind,2] <- NA
    }
    
    if(substr(coeff_names, start = 1, stop = 4)=="tmax"){
      coeff_name_date[ind,1] <- "tmax"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 6, stop = 11)
    }
    
    if(substr(coeff_names, start = 1, stop = 1)=="v"){
      coeff_name_date[ind,1] <- "vpd"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 5, stop = 10)
    }
    
    if(substr(coeff_names, start = 1, stop = 1)=="p"){
      coeff_name_date[ind,1] <- "pr"
      coeff_name_date[ind,2] <- substr(coeff_names, start = 4, stop = 9)
    }
    
  }#end for ind
  
  return(coeff_name_date)
  
}#end function get_firstcoeffs




# Additional functions, useful to plot Fig 7. ####

# count number of variable kept (Fig. 7a.)
number_coeff_kept <- function(coeff){ # give number of coeff !=0
  return(length(which(abs(coeff[-1,])>0)))
}#end function number_coeff_kept

# count number of selected extreme indices (Fig. 7b)
extreme_in_coeff <- function(coeff){ #function to check how many extreme indeices are selected as predictors
  extreme_indices <- c("dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p")
  if(max(abs(coeff[extreme_indices,]))==0){
    return(0)
  } else {
    return(length(which(abs(coeff[extreme_indices,])>0)))
  }
}#end function extreme_in_coeff

# Count the number of selected seasons for each grid point (Fig. 7d)
count_seasons <- function(coeff){
  if (length(which((coeff)!=0))-1>0) {
    coeff_kept <- get_firstcoeffs(coeff, nb_of_coeff = length(which((coeff)!=0))-1)
    nb_of_seas <- 0
    # Winter
    if(sum(!is.na(coeff_kept[,2])) & (sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Feb", "Dec", "Jan")))){
      nb_of_seas <- nb_of_seas + 1
    }
    # Spring
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("May", "Mar", "Apr"))){
      nb_of_seas <- nb_of_seas + 1
    }
    # Summer
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Jun", "Jul", "Aug"))){
      nb_of_seas <- nb_of_seas + 1
    }
    # Autumn
    if(sum(!is.na(coeff_kept[,2])) & sum(substr(coeff_kept[,2], start = 1, stop = 3) %in% c("Sep", "Nov", "Oct"))){
      nb_of_seas <- nb_of_seas + 1
    }
    if (nb_of_seas>0){
      return(nb_of_seas)
    } else {
      return("No met. var")
    }
  } else {
    return("No met. var")
  }#end if else
} #end function count_seasons