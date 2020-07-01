# Identify_crop_yield_drivers
This repository provides the code of the article "Identifying meteorological drivers of extreme impacts: an application to simulated crop yields"

This readme gives an overview of the structure of the code related to the submitted manuscript of “Identifying meteorological drivers of extreme impacts: an application to simulated crop yields”.
The codes were developed using R version 3.6 and Python version 3.7.

To run the code you need to assign the path where the code and the input data are stored on your computer to the variables “path_code” ,"path_data” and "data_dir" within the respective files Lasso_regression.R, Figures.R and Plot_composites.py. 
Lasso_regression.R calculates the Lasso logistic regression model. 
Figures.R creates the figures 1, 3, 5, 6, 7, 8, A3, A4 and the supplementary gifs of the article. 
Plot_composites.py creates the figures 2, A1 and A2 of the article.
Data_processing.R, correlation_computation.R, and additional_functions.R are subroutines required for Lasso_regression.R and Figures.R.
Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pix.RData contains the results of the Lasso logistic regression model using a lambda with standard error for the 995 grid points. It was created using Lasso_regression.R.

The corresponding climate and crop simulations to run the code are available from Tianyi Zhang (zhangty@post.iap.ac.cn) and Karin van der Wiel (wiel@knmi.nl) on request, respectively. Using the results of the Lasso logistic regression model (Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pix.RData, see description above) it is possible to create the figures 5, 7, 8, A3, A4 and the supplementary gifs with Figures.R without requiring the climate and crop simulation data.
