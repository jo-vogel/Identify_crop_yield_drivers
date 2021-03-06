# Identify_crop_yield_drivers
This repository provides the code of the article "Identifying meteorological drivers of extreme impacts: an application to simulated crop yields".
The code was created with contributions from Johannes Vogel, Pauline Rivoire, Cristina Deidda, Christoph A. Sauter and Elisabeth Tschumi.

This readme gives an overview of the structure of the code related to the submitted manuscript of “Identifying meteorological drivers of extreme impacts: an application to simulated crop yields”.
The codes were developed using R version 3.6 and Python version 3.7.

To run the code you need to assign the path where the code and the input data are stored on your computer to the variables “path_code” ,"path_data” and "data_dir" within the respective files Lasso_regression.R, Figures.R and Plot_composites.py. 
Lasso_regression.R calculates the Lasso logistic regression model. 
Figures.R creates the figures 1, 3, 5, 6, 7, 8, B4, B5 and the supplementary gifs of the article. 
Plot_composites.py creates the figures 2, B2 and B3 of the article.
Data_processing.R, correlation_computation.R, and additional_functions.R are subroutines required for Lasso_regression.R and Figures.R.
Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pix.RData contains the results of the Lasso logistic regression model using a lambda with standard error for the 995 grid points. It was created using Lasso_regression.R.
Required_variables.RData contains a list of data required for plotting Fig 3, 5, 8, B4 and the supplementary GIFs when running the file Figures.R without the crop and climate data. They were extracted as stated in the last line in Figures.R.

The corresponding climate and crop simulations to run the code are available from Tianyi Zhang (zhangty@post.iap.ac.cn) and Karin van der Wiel (wiel@knmi.nl) on request, respectively. The authors can then provide the pre processed data (extremeindices_and_monthlymeteovar_rescaled_995pix.Rdata) to the reader. Using the results of the Lasso logistic regression model (Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pix.RData, see description above) it is possible to create the figures 5, 7, 8, B4, B5 and the supplementary gifs with Figures.R without requiring the climate and crop simulation data.

Note that the terminology is slightly differs in the code compared to the article. In the code, "0" refers to a bad year and "1" refers to a normal year, while in the article it is vice versa.
