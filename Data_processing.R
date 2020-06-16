# This file contains the data processing, used in the Lasso regression and when plotting the figures 
# Authors: Johannes Vogel, Pauline Rivoire

# This code was run under version 3.6 of R

pix_num <- length(Data_xtrm_standardized$longitudes) # Number of all wheat growing pixels in the northern hemisphere
nb_years <- dim(Data_xtrm_standardized$yield)[2]
# An additional dimension is added to be able to combine the 2-dimensional variables (grid points and years) with the 3-dimensional variables (grid points, months and years)
yield_3dim <- array(Data_xtrm_standardized$yield,dim=c(pix_num,1,nb_years))
dtr_3dim <- array(Data_xtrm_standardized$dtr,dim=c(pix_num,1,nb_years))
frs_3dim <- array(Data_xtrm_standardized$frs,dim=c(pix_num,1,nb_years))
txx_3dim <- array(Data_xtrm_standardized$txx,dim=c(pix_num,1,nb_years))
tnn_3dim <- array(Data_xtrm_standardized$tnn,dim=c(pix_num,1,nb_years))
rx5_3dim <- array(Data_xtrm_standardized$rx5,dim=c(pix_num,1,nb_years))
tx90p_3dim <- array(Data_xtrm_standardized$tx90p,dim=c(pix_num,1,nb_years))
tn10p_3dim <- array(Data_xtrm_standardized$tn10p,dim=c(pix_num,1,nb_years))

# Combine all variables in one dataset: crop yield, climate extreme indicators, monthly mean meteorological variables
Model_data <- abind(yield_3dim,dtr_3dim,frs_3dim,txx_3dim,tnn_3dim,rx5_3dim,tx90p_3dim,tn10p_3dim
                    ,Data_xtrm_standardized$tasmax,Data_xtrm_standardized$vpd,Data_xtrm_standardized$pr,along=2)
colnames(Model_data) <- c("Yield", "dtr", "frs", "txx", "tnn", "rx5", "tx90p", "tn10p",
                          "tmax_Aug_Y1", "tmax_Sep_Y1", "tmax_Oct_Y1", "tmax_Nov_Y1", "tmax_Dec_Y1", "tmax_Jan_Y2",
                          "tmax_Feb_Y2", "tmax_Mar_Y2", "tmax_Apr_Y2", "tmax_May_Y2", "tmax_Jun_Y2", "tmax_Jul_Y2",
                          "tmax_Aug_Y2", "tmax_Sep_Y2", "tmax_Oct_Y2", "tmax_Nov_Y2", "tmax_Dec_Y2",
                          "vpd_Aug_Y1", "vpd_Sep_Y1", "vpd_Oct_Y1", "vpd_Nov_Y1", "vpd_Dec_Y1", "vpd_Jan_Y2",
                          "vpd_Feb_Y2", "vpd_Mar_Y2", "vpd_Apr_Y2", "vpd_May_Y2", "vpd_Jun_Y2", "vpd_Jul_Y2",
                          "vpd_Aug_Y2", "vpd_Sep_Y2", "vpd_Oct_Y2", "vpd_Nov_Y2", "vpd_Dec_Y2",
                          "pr_Aug_Y1", "pr_Sep_Y1", "pr_Oct_Y1", "pr_Nov_Y1", "pr_Dec_Y1", "pr_Jan_Y2",
                          "pr_Feb_Y2", "pr_Mar_Y2", "pr_Apr_Y2", "pr_May_Y2", "pr_Jun_Y2", "pr_Jul_Y2",
                          "pr_Aug_Y2", "pr_Sep_Y2", "pr_Oct_Y2", "pr_Nov_Y2", "pr_Dec_Y2")


Yield <- Data_xtrm_standardized$yield
threshold <- 0.05 # threshold for bad yields
low_yield <- apply(Yield, MARGIN = 1, FUN=quantile, probs=threshold, na.rm=T)
cy <- t(sapply(1:pix_num,function(x){ifelse(Yield[x,]<low_yield[x],0,1)})) # assign crop yield to the two class bad and normal yield years
cy_reshaped <- array(data=cy,dim=c(dim(cy)[1],1,nb_years))
Model_data[,1,] <- cy_reshaped # Replace continuous crop yield with binary categories (bad and normal yield)


# Exclude NA variable columns (columns are NA if the respective months are outside of the maximum growing season of the corresponding grid point)
na_col <- matrix(data=NA,nrow=pix_num,ncol=dim(Model_data)[2])
for (j in 1:pix_num){
  for (i in 1:dim(Model_data)[2]){
    na_col[j,i] <- all(is.na(Model_data[j,i,])) # TRUE if entire column is NA
  }
}
non_na_col <- !na_col # columns without NAs
non_na_col[,1] <- FALSE # exclude yield (it is no predictor and should therefore be ignored)

# Exclude years with NAs (years are NA if the growing season exceeds 365 days)
na_time <- vector("list",length=pix_num) # for each grid point, the positions of NAs over time
for (j in 1:pix_num){
  na_time[[j]] <- which(is.na(Model_data[j,1,])) # locations of years with NA values
}

years_with_na <- vector("logical",length=pix_num)
for (i in 1:pix_num){
  years_with_na[i] <- ifelse(length(na_time[[i]]) == 0,F,T)
}

##### Split data into training and testing data set #####
training_indices <- vector("list",length=pix_num)
testing_indices <- vector("list",length=pix_num)
seed <- 1994
train_size <- 70
set.seed(seed)
for (x in 1:pix_num) {
  if (years_with_na[x]) {
    training_indices[[x]] <- sort(sample(x=(1:nb_years)[-na_time[[x]]], size = floor((nb_years-length(na_time[[x]]))*(train_size/100)))) # Assign percentage of years according to train_size after neglecting years with NAs
    testing_indices[[x]] <- (1:nb_years)[-c(na_time[[x]], training_indices[[x]])] # All other years without NAs are assigned to the testing indices
  } else { # Simplified assignment in case there is no need to account for years with NAs
    training_indices[[x]] <- sort(sample(1:nb_years, size = floor(nb_years*(train_size/100))))
    testing_indices[[x]] <- (1:nb_years)[-training_indices[[x]]]    
  }
}

Training_Data <- lapply(1:pix_num,function(x){Model_data[x,,training_indices[[x]]]})
Testing_Data <- lapply(1:pix_num,function(x){Model_data[x,,testing_indices[[x]]]})

# Split in training and testing predictors and predictands
pix_in <- 1:pix_num
x1_train_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Training_Data[[x]][non_na_col[x,],]))}) # predictors
y1_train_list <- lapply(seq_along(pix_in), function(x){ Training_Data[[x]][1,]}) # predictand
x1_test_list <- lapply(seq_along(pix_in), function(x){ as.data.frame(t(Testing_Data[[x]][non_na_col[x,],]))}) # predictors
y1_test_list <- lapply(seq_along(pix_in), function(x){Testing_Data[[x]][1,]}) # predictand
