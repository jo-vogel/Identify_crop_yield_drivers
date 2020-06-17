# Identify_crop_yield_drivers
This repository provides the code of the article "Identifying meteorological drivers of extreme impacts: an application to simulated crop yields"

This readme gives an overview of the structure of the code related to the submitted manuscript of “Identifying meteorological drivers of extreme impacts: an application to simulated crop yields”.
The codes were developed using R version 3.6 and Python version 3.7.

To run the code you need to assign the path where the code and the input data are stored on your computer to the variables “path_code” ,"path_data” and "data_dir" within the respective files Lasso_regression.R, Figures.R and and Plot_composites.py. 
Lasso_regression.R calculates the Lasso logistic regression model. 
Figures.R creates the figures 1, 5, 6, 7, 8, A3, A4 and the supplementary gifs of the article. 
Plot_composites.py creates the figures 2, A1 and A2 of the article.
Data_processing.R and additional_functions.R are subroutines required for Lasso_regression.R and Figures.R.

The corresponding data to this repository can be found here: https://drive.google.com/drive/folders/1dVlsKwll1wtBGHcxjMzwgCl7u-oidw__
