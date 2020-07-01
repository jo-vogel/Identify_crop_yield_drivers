# Creation of the Figures 1, 5, 6, 7, 8, A3, A4 and the supplementary gifs of the article
# Authors: Pauline Rivoire, Johannes Vogel, Cristina Deidda, Elisabeth Tschumi

# This code was run under version 3.6 of R

#' Structure of the code
#' a) Load data
#' b) Process data
#' c) Create Fig. 1
#' d) Create Fig. 3
#' e) Create Fig. 5
#' f) Create Fig. 6
#' g) Create Fig. 7
#' h) Create Fig. 8
#' i) Create Fig. A3
#' j) Create Fig. A4
#' k) Create GIFs


##### Load data ####
####################

library(glmnet);library(InformationValue);library(ROCR);library(ggpubr);library(abind);library(maps);library(corrplot)
library(oce);library(stringr);library(ggplot2);library(viridis);library(raster);library(rgdal);library(pbapply)

path_data <- message("insert data directory here") # Insert path of the input data here
path_code <- message("insert code directory here") # Insert path of the code here
seed=1994 # random seed
train_size <- 70 # Percentage of data assigned to the training data set

message("Precise here if climate and crop data are available. The corresponding climate and crop simulations to run the code are available 
from Tianyi Zhang (zhangty@post.iap.ac.cn) and Karin van der Wiel (wiel@knmi.nl) on request, respectively. Using the results of the Lasso 
logistic regression model (Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed1994_train70_995pix.RData) it is 
possible to create the figures 5, 7, 8, A3, A4 and the supplementary gifs with Figures.R without requiring the climate and crop simulation data.")
climate_crop_provided <- FALSE


if (climate_crop_provided) {
  # Load the preprocessed standardized climate variables and crop yield for all 995 grid points in the northern hemisphere
  load(paste0(path_data,"/extremeindices_and_monthlymeteovar_rescaled_995pix.Rdata")) 
  # Mean yield and yield standard deviation for all 995 grid points
  load(paste0(path_data,"/RawMeanYield_995pix.Rdata"))
  load(paste0(path_data,"/RawSdYield_995pix.Rdata"))
} #end if crop yield provided


# Final selection of grid points
load(paste0(path_data,"/final_889pix_coords.Rdata"))
# Load matrix with all coordinates required for Fig. 8
coord_all <- read.csv2(paste0(path_data,"/coord_all.csv"), row.names=1)
# Shapefile of borders of the continents
continents <- readOGR(paste0(path_data,"/continent.shp")) 
message("This shapefile can be downloaded from https://www.arcgis.com/home/item.html?id=5cf4f223c4a642eb9aa7ae1216a04372.")

# The statistical model is calculated using  Lasso_regression.R
load(paste0(path_data,"/Lasso_lambda1se_month_xtrm_LASSO_threshbadyield005_seed",seed, "_train", train_size,"_995pix.Rdata"))
Model_chosen <- lasso_model_lambda1se # Lasso regression using lambda 1 standard error

# Load additional required functions
source(paste0(path_code,"/additional_functions.R"))


if (climate_crop_provided) {
  # Preprocess the data
  source(paste0(path_code,"/Data_processing.R"))
} else { # Load the preprocessed data
  load(paste0(path_data,"/Required_variables.RData"))
}



##### Adjust cutoff level #####

if (climate_crop_provided) {
  y1_train_list_simple_lasso <- y1_train_list
  x1_train_list_simple_lasso <- x1_train_list
  Model_chosen_889 <- list()
  y1_train_list_simple_lasso <- list()
  x1_train_list_simple_lasso <- list()
  work_pix_tmp <- numeric()
  final_pix_num <- length(final_pixels_coord$latitude) # see section 2.2 of the article
  
  for (pixel in 1:final_pix_num) {
    pix_in_995 <- final_pixels_coord$ref_in_995[pixel]
    y1_train_list_simple_lasso[[pixel]] <- y1_train_list[[pix_in_995]]
    x1_train_list_simple_lasso[[pixel]] <- x1_train_list[[pix_in_995]]
    Model_chosen_889[[pixel]] <- Model_chosen[[pix_in_995]]
    if(is.character(Model_chosen[[pix_in_995]])){work_pix_tmp[pixel]<-0} else {work_pix_tmp[pixel]<-1}
  }#end for pixel
  cost_fp_simple_lasso <- 100 # Misses: this should be associated with a higher cost, as it is more detrimental
  cost_fn_simple_lasso <- 100 # False alarms
  work_pix <- which(work_pix_tmp==1)
  
  # return the mean value, over all pixels, of the adjusted cutoff named segregation threshold
  segreg_th <- adjust_cutoff(model_vector = Model_chosen_889,x1_train_list = x1_train_list_simple_lasso, y1_train_list = y1_train_list_simple_lasso,
                             work_pix = work_pix, cost_fp = cost_fp_simple_lasso, cost_fn= cost_fn_simple_lasso)
} #end if crop yield provided


# General figure variables 
world <- map_data("world")
coord_subset <- cbind(final_pixels_coord$longitude, final_pixels_coord$latitude)
final_pix_num <- dim(coord_subset)[1]
ewbrks <- seq(-100,100,50)
nsbrks <- seq(10,50,10)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(abs(x), "°W"), ifelse(x > 0, paste(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(abs(x), "°S"), ifelse(x > 0, paste(x, "°N"),x))))


# Figure 1: Mean annual yield ####
##################################

if (climate_crop_provided) {
  DF_meanY <- data.frame(lon=Raw_mean_yield[,"longitudes"], lat = Raw_mean_yield[,"latitudes"],
                         meany = Raw_mean_yield[,"mean_yield"]/1000) # data frame containing mean annual yield (transferred from kg to tonnes) and associated coordinates
  pixels_excluded <- as.logical(1-(1:pix_num %in% final_pixels_coord$ref_in_995)) # Excluded grid points according to section 2.2 of the article
  DF_excluded_pix <- data.frame(lon = Raw_mean_yield[pixels_excluded,"longitudes"],
                                lat = Raw_mean_yield[pixels_excluded,"latitudes"])
  
  ggplot(data = DF_meanY, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3)+  geom_tile(aes(fill=DF_meanY$meany)) +
    scale_fill_gradient2(midpoint = max(DF_meanY$meany, na.rm = T)/2,
                         limits=c(0,max(DF_meanY$meany)),
                         low = "#f7fcb9", mid = "#addd8e", high = "#31a354") +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15))+
    # ylab(expression("Lat " ( degree*N))) +
    # xlab(expression("Lon " ( degree*E))) +
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
                ylim = c(min(DF_meanY$lat), max(DF_meanY$lat)),
                ratio = 1)+
    labs(fill=expression(paste("Mean yield [t ", ha^{-1},"]")))+
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14),
          axis.title.x=element_blank(),axis.title.y=element_blank())+
    geom_point(data = DF_excluded_pix, aes(x = DF_excluded_pix$lon, y = DF_excluded_pix$lat),
               color = "black", size = 0.89, pch=4)
  ggsave(filename = "Raw_mean_yield.png", width = 20, height = 4)
  
}#end if crop yield provided



# Figure 3: Linear correlation plot ####
########################################
if (climate_crop_provided) {
  # Compute correlation /!\ might take  ~10min
  source(paste0(path_code,"/correlation_computation.R"))
} else {
  load(paste0(path_data,"/correlation_vectors.Rdata"))
  france_meteovar_correlations <- list_correlations$france_meteovar_correlations
  france_xtrm_correlations <- list_correlations$france_xtrm_correlations
  global_meteovar_correlations <- list_correlations$global_meteovar_correlations
  global_xtrm_correlations <- list_correlations$global_xtrm_correlations
} #end if else crop yield provided



# Set p lot layout
pdf(file="correlation_plot.pdf", width = 12, height = 4)
layout(matrix(c(1, 2, 3, 4,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 5, 5, 5), ncol=4, byrow=TRUE),
       heights = c(1,1,1,0.5),    
       widths = c(8, 2.3, 9.5, 2.3))     

#define margins and colorpalette
par(oma=c(0,0,2,0),mai=c(0.2,0,0.2,0),cex=1, xpd=NA)
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))
#france meteovars
corrplot(france_meteovar_correlations[,c(4:14,18:22)], method = "circle", tl.col="black", xpd=NA, cl.pos = "n", cl.ratio=0.3, col = col2(50))
segments(11.5,0.5,x1=11.5,y1=12.5)
segments(15.5,0.5,x1=15.5,y1=12.5)
mtext("(a)",line=1,at=-0.9, cex=1.1)
#France extremal indicators
colnames(france_xtrm_correlations)<-c("GS                 ")
corrplot(france_xtrm_correlations, method = "circle", tl.col="black", col = col2(50), cl.pos = "n", cl.ratio=1)
mtext("(b)",line=1,at=-1.4, cex=1.1)
segments(2,0,2,11,lwd=3)
#global meteovars
corrplot(global_meteovar_correlations[,c(2:16,18:22)], method = "circle", tl.col="black", col = col2(50), cl.ratio=1, cl.pos="n")
mtext("(c)",line=1,at=-0.9, cex=1.1)
segments(15.5,0.5,x1=15.5,y1=12.5)
segments(19.5,0.5,x1=19.5,y1=12.5)
#global extremal indicators
colnames(global_xtrm_correlations)<-c("GS                 ")
corrplot(global_xtrm_correlations, method = "circle", tl.col="black", col = col2(50), cl.pos = "n", cl.ratio=1)
mtext("(d)",line=1,at=-1.4, cex=1.1)


par(mai=c(0,0.7,0,0.1))
drawPalette(zlim=c(-1,1), col = col2(50), pos=1, las=1)
dev.off()



# Figure 5: Critical success index (CSI) ####
#############################################

coeff  <-list() # List of all predictors of a given pixel
csi <- rep(NA, final_pix_num) # Critical success index
mypred <- vector("list", final_pix_num) # predicted crop yield
fitted_bad_yield <- vector("list", final_pix_num) # assign predicted crop yield to either bad or normal year crop yield

for (pixel in 1:final_pix_num) { #extract coefficients
  pix_in_995 <- final_pixels_coord$ref_in_995[pixel]
  coeff[[pixel]] <- coefficients(Model_chosen[[pix_in_995]])
}#end for pixel


if (climate_crop_provided) {
  for (pixel in 1:final_pix_num) { #extract predicted yield and csi
    pix_in_995 <- final_pixels_coord$ref_in_995[pixel]
    mypred[[pixel]] <- predict(Model_chosen[[pix_in_995]], as.matrix(x1_test_list[[pix_in_995]]),type="response")
    fitted_bad_yield[[pixel]] <- ifelse(mypred[[pixel]] > segreg_th,1,0)
    con_tab <- InformationValue::confusionMatrix(actuals = as.matrix(y1_test_list[[pix_in_995]]),
                                                 predictedScores = fitted_bad_yield[[pixel]], threshold = segreg_th)
    csi[pixel] <- con_tab["0","0"]/(con_tab["0","0"] + con_tab["1","0"] + con_tab["0","1"])
    if(is.na(con_tab["0","0"])){
      csi[pixel] <- 0
    }
  }#end for pixel
} else {
  load(paste0(path_data,"/csi_vector.Rdata"))
} #end if else crop yield provided

DF_csi <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], csi = csi)

ggplot(data = DF_csi, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=csi)) +
  scale_fill_gradient2(midpoint = max(csi, na.rm = T)/2,
                       low = "black", mid = "red3", high = "yellow") +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="CSI"
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave(filename = "CSImap_Lasso_lambda1se_adjcutoff_seed1994_training70_889GP.png", width = 20, height = 3.5)



# Figure 6: Correlation between Critical Success Index (CSI) and annual crop yield mean and variability ####
############################################################################################################

if (climate_crop_provided) {
  MeanY_CSI<-data.frame(cbind(Raw_mean_yield[final_pixels_coord$ref_in_995,"mean_yield"]/1000,csi))
  colnames(MeanY_CSI)<-c("mean_yield","csi") # data frame containing mean annual yield (transferred from kg to tonnes) and csi
  
  SDY_CSI<-data.frame(cbind(Raw_sd_yield[final_pixels_coord$ref_in_995,"sd_yield"]/1000,csi))
  colnames(SDY_CSI)<-c("sd_yield","csi") # data frame containing mean yield standard deviation (transferred from kg to tonnes) and csi
  
  p1<-ggplot(MeanY_CSI, aes(x=mean_yield, y=csi)) + geom_point()+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 17)) +
    labs( x = expression(paste("Mean yield [t ", ha^{-1},"]")))  + ylab("CSI") + theme(panel.grid.major = element_blank(),
                                                                                       plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"),
                                                                                       panel.grid.minor = element_blank(), 
                                                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  p2<-ggplot(SDY_CSI, aes(x=sd_yield, y=csi)) + geom_point()+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 17)) +
    labs( x = expression(paste("Mean yield [t ", ha^{-1},"]"))) + ylab("CSI")+ theme(panel.grid.major = element_blank(),
                                                                                     plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"), panel.grid.minor = element_blank(), 
                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggarrange(p1, p2, nrow = 1,ncol=2,labels = c("(a)", "(b)" ), font.label = list(size = 21, face="plain"))
} #end if crop yield provided

ggsave(filename = "Scatterplot_CSIvsYield_Lasso_lambda1se_adjcutoff_seed1994_training70_889GP.png", width = 15, height = 7)                      


                        
# Figure 7: Maps illustrating the selected predictors by the Lasso logistical regression ####
#############################################################################################

message("run Fig. 5 for this section first")
# Plot number of selected variables and climatic extreme indicators (Fig 7a and b) ####


nb_extr_kept <- numeric()
nb_coeff_kept <- numeric()
for (pixel in 1:final_pix_num) {
  nb_extr_kept[pixel] <- extreme_in_coeff(coeff[[pixel]])
  nb_coeff_kept[pixel] <- number_coeff_kept(coeff[[pixel]])
}#end grid point

levels_nb_var <- cut(nb_coeff_kept, breaks = c(0,5,10,15,20,25,30), right = F)
levels_nb_var <- gsub(","," - ",levels_nb_var,fixed=TRUE)


# Combination of met. variables (Fig. 7c) ####
met_type <- numeric(length=length(coeff)) # Associated numbers of all possible combinations of the 3 monthly mean predictors VPD, Tmax and Pr
meteo_type <- as.character(met_type) # Names of all possible combinations of the 3 monthly mean predictors VPD, Tmax and Pr
met_strings <- vector("list",length=length(coeff)) # Character string of predictor names
met_vpd <- numeric(length=length(coeff)) # Assign all cases of VPD selection as predictor
met_pr <- numeric(length=length(coeff)) # Assign all cases of Pr selection as predictor
met_temp <- numeric(length=length(coeff)) # Assign all cases of Tmax selection as predictor

for (pix in 1:length(coeff)) {
  if((length(which(coeff[[pix]]!=0))-1)<1){
    met_type[pix] <- 8
    meteo_type[pix] <- "None"
  } else {
    coeff_kept <- get_firstcoeffs(coeff = coeff[[pix]],
                                  nb_of_coeff = (length(which(coeff[[pix]]!=0))-1))
    met_string <- (c("vpd","pr","tmax") %in% unique(coeff_kept[,1])) 
    met_strings[[pix]] <- met_string
    if (identical(met_string, c(T,F,F))) {
      met_type[pix] <- 1
      meteo_type[pix] <- "VPD"
    } else if (identical(met_string, c(F,T,F))) {
      met_type[pix] <- 2
      meteo_type[pix] <- "Pr"
    } else if (identical(met_string, c(F,F,T))) {
      met_type[pix] <- 3
      meteo_type[pix] <- "Tmax"
    } else if (identical(met_string, c(T,T,F))) {
      met_type[pix] <- 4
      meteo_type[pix] <- "VPD & Pr"
    } else if (identical(met_string, c(T,F,T))) {
      met_type[pix] <- 5
      meteo_type[pix] <- "VPD & Tmax"
    } else if (identical(met_string, c(F,T,T))) {
      met_type[pix] <- 6
      meteo_type[pix] <- "Pr & Tmax"
    } else if (identical(met_string, c(T,T,T))) {
      met_type[pix] <- 7
      meteo_type[pix] <- "All"
    } else if (identical(met_string, c(F,F,F)) | is.null(met_string)) {
      met_type[pix] <- 8
      meteo_type[pix] <- "None"
    }
    if ("vpd" %in% unique(coeff_kept[,1])) {met_vpd[pix] <- 1} 
    if ("pr" %in% unique(coeff_kept[,1])) {met_temp[pix] <- 2} 
    if ("tmax" %in% unique(coeff_kept[,1])) {met_pr[pix] <- 3}
  }#end ifelse
}#end for pix


# Count the number of selected seasons for each grid point (Fig. 7d) ####

nb_of_seas <- numeric()
for (pix in 1:final_pix_num) {
  nb_of_seas[pix] <- count_seasons(coeff[[pix]])
}

# Fig. 7a: Total number of selected variables
DF_numbcoeff <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = levels_nb_var)
DF_numbcoeff$levels_nb_var <- gsub("\\[|\\)","",levels_nb_var)

P1 <- ggplot(data = DF_numbcoeff, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_numbcoeff$levels_nb_var)) +
  scale_fill_manual(values=c("0 - 5"="#f1eef6", "5 - 10"="#d4b9da", "10 - 15"="#c994c7",
                             "15 - 20"="#df65b0", "20 - 25"="#dd1c77", "25 - 30"="#980043"),
                    breaks=c("0 - 5", "5 - 10", "10 - 15", "15 - 20", "20 - 25", "25 - 30"),
                    label=c("0 - 4", "5 - 9", "10 - 14", "15 - 10", "20 - 24", "25 - 29"))+
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Nb of var."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))


# Fig. 7b: Plot number of selected extreme indices
DF_numbextr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], coeff_kep = nb_extr_kept)
DF_numbextr$coeff_kep <- as.factor(DF_numbextr$coeff_kep)

mycolors <- rev(hcl.colors(n=8,palette="Viridis"))

P2 <- ggplot(data = DF_numbextr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=coeff_kep)) +
  scale_fill_manual(values = mycolors, ) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Nb extr.\n ind."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))


# Fig. 7c: Combinations of meteorological variables
DF_meteo_type <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], met_type = meteo_type)

cols <- c("VPD" = "#7FC97F", "Pr" = "cadetblue2", "Tmax" = "#386CB0", "VPD & Pr" = "#824D99",
          "VPD & Tmax" = "#F0027F", "Pr & Tmax" = "darkred" , "All" = "#FDC086", "None" = "#FFFF99")

P3 <- ggplot(data = DF_meteo_type, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=met_type)) +
  scale_fill_manual(values = cols,breaks=c("VPD","Pr","Tmax","VPD & Pr","VPD & Tmax","Pr & Tmax", "All","None")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Combination\nof met.\nvariables"
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))


DF_nbseason <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], nb_season = nb_of_seas)
DF_nbseason$nb_season <- as.factor(DF_nbseason$nb_season)


# Fig. 7d: Plot number of selected seasons
P4 <- ggplot(data = DF_nbseason, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_nbseason$nb_season)) +
  scale_fill_manual(values = c("0"=viridis(6)[1], "1"=viridis(6)[3], "2"=viridis(6)[4],
                               "3"=viridis(6)[5], "4"=viridis(6)[6], "No met. var"="pink")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Nb of seas."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))

# Combine 4 subplots to one plot
L1 <- get_legend(P1+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L2 <- get_legend(P2+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L3 <- get_legend(P3+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L4 <- get_legend(P4+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))

ggarrange(P1 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L1,
          P2 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L2,
          P3 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L3,
          P4 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L4,
          nrow = 4, ncol=2,labels = c("(a)", "", "(b)", "",
                                      "(c)", "", "(d)", ""),
          label.x = -0.015
          ,widths=c(8,1), heights=c(1,1,1), font.label = list(size = 14, face = "plain", color ="black")
)
ggsave(filename = "4maps_nbvar_nbextr_combin-meteovar_nbseas_Lasso_lambda1se_adjcutoff_seed1994_training70_889GP.png", width = 30, height = 20)



# Figure 8: Percentage of grid points, where predictor is selected ####
#######################################################################

# Find indices of the coordinates
coord_subset_temp <- cbind(coord_subset,paste(coord_subset[,1],coord_subset[,2]))
coord_all_temp <- cbind(coord_all,paste(coord_all[,1],coord_all[,2]))
loc_pix <- which(coord_all_temp[,3] %in% coord_subset_temp [,3]) # locations of our grid points in the whole coordinate set

# Create continent polygons
africa <- subset(continents,subset=continents@data[["CONTINENT"]]=="Africa")
europe <- subset(continents,subset=continents@data[["CONTINENT"]]=="Europe")
no_am <- subset(continents,subset=continents@data[["CONTINENT"]]=="North America")
asia <- subset(continents,subset=continents@data[["CONTINENT"]]=="Asia")

# Connect final grid points to their coordinates
num_lon <- 320 # number of longitudes
num_lat <- 76 # number of latitudes
coord_assigned <- cbind(coord_all,rep(NA,num_lon*num_lat))
for (i in seq_along(1:final_pix_num)){
  coord_assigned[loc_pix[i],3] <- i
}

# Extract grid points by continent
loc_mat <- matrix(as.numeric(coord_assigned[,3]),nrow=num_lon,ncol=num_lat)
loc_ras <- raster(t(loc_mat[,num_lat:1]), xmn=min(coord_all$long_all), xmx=max(coord_all$long_all), ymn=min(coord_all$lati_all), ymx=max(coord_all$lati_all), crs=CRS(projection(continents)))
loc_afr_pixels <- extract(loc_ras,africa)
loc_eur_pixels <- extract(loc_ras,europe)
loc_no_am_pixels <- extract(loc_ras,no_am)
loc_asia_pixels <- extract(loc_ras,asia)

loc_afr_pixels_num <- loc_afr_pixels[[1]][!is.na(loc_afr_pixels[[1]])] 
loc_eur_pixels_num <- loc_eur_pixels[[1]][!is.na(loc_eur_pixels[[1]])]
loc_no_am_pixels_num <- loc_no_am_pixels[[1]][!is.na(loc_no_am_pixels[[1]])]
loc_asia_pixels_num <- loc_asia_pixels[[1]][!is.na(loc_asia_pixels[[1]])]

sum(!is.na(loc_eur_pixels[[1]])) #  233 grid points in Europe
sum(!is.na(loc_no_am_pixels[[1]])) #  419 grid points in North America
sum(!is.na(loc_afr_pixels[[1]])) # 22 grid points in Africa
sum(!is.na(loc_asia_pixels[[1]])) # 210 grid points in Asia
# 22+233+419+210=884; 889-884: 5 grid point are missing
# Add missing points
loc_eur_pixels_num <- c(loc_eur_pixels_num,450)
loc_no_am_pixels_num <- c(loc_no_am_pixels_num,391)
loc_asia_pixels_num <- c(loc_asia_pixels_num,3, 294, 563)

# Get data into right format for the plot
if (climate_crop_provided) {
  var_num <- apply(non_na_col,1,sum) # Number of predictor variables for each grid point
  numLevels_list <- sapply(1:pix_num, function(x){ rep(1,times=var_num[x])})
  for (i in 1:pix_num){
    names(numLevels_list[[i]]) <-  colnames(x1_test_list[[i]])
  }
}

coefs_seas <- sapply(1:length(coeff), function(x) names(numLevels_list[[final_pixels_coord$ref_in_995[x]]])[coeff[[x]][-1]!=0])
coefs_seas_vec <- unlist(coefs_seas)
coefs_seas_afr <- sapply(loc_afr_pixels_num, function(x) names(numLevels_list[[final_pixels_coord$ref_in_995[x]]])[coeff[[x]][-1]!=0])
coefs_seas_afr_vec <- unlist(coefs_seas_afr )
coefs_seas_eur <- sapply(loc_eur_pixels_num, function(x) names(numLevels_list[[final_pixels_coord$ref_in_995[x]]])[coeff[[x]][-1]!=0])
coefs_seas_eur_vec <- unlist(coefs_seas_eur )
coefs_seas_no_am <- sapply(loc_no_am_pixels_num, function(x) names(numLevels_list[[final_pixels_coord$ref_in_995[x]]])[coeff[[x]][-1]!=0])
coefs_seas_no_am_vec <- unlist(coefs_seas_no_am )
coefs_seas_asia <- sapply(loc_asia_pixels_num, function(x) names(numLevels_list[[final_pixels_coord$ref_in_995[x]]])[coeff[[x]][-1]!=0])
coefs_seas_asia_vec <- unlist(coefs_seas_asia )

coefs_seas_afr_tab <- as.data.frame(table(coefs_seas_afr_vec))
coefs_seas_eur_tab <-  as.data.frame(table(coefs_seas_eur_vec))
coefs_seas_no_am_tab <- as.data.frame(table(coefs_seas_no_am_vec))
coefs_seas_asia_tab <- as.data.frame(table(coefs_seas_asia_vec))
coefs_seas_vec_tab <- as.data.frame(table(coefs_seas_vec))

# correct variable names
list_coefs <- list(coefs_seas_afr_tab,coefs_seas_eur_tab,coefs_seas_no_am_tab,coefs_seas_asia_tab,coefs_seas_vec_tab)
for (i in 1:length(list_coefs)){
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="vpd", replacement = "VPD")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="pr", replacement = "Pr")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="APr", replacement = "Apr") # recorrect April
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tmax", replacement = "Tmax")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="txx", replacement = "TXx")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tnn", replacement = "TNn")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="rx5", replacement = "Rx5day")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tx90p", replacement = "TX90p")
  list_coefs[[i]][,1] <- gsub(x=list_coefs[[i]][,1], pattern="tn10p", replacement = "TN10p")
}
coefs_seas_afr_tab <- list_coefs[[1]];coefs_seas_eur_tab <- list_coefs[[2]];coefs_seas_no_am_tab <- list_coefs[[3]];coefs_seas_asia_tab <- list_coefs[[4]];coefs_seas_vec_tab <- list_coefs[[5]]
colnames(coefs_seas_afr_tab) <- c("Variables","Freq_Afr")
colnames(coefs_seas_eur_tab) <- c("Variables","Freq_Eur")
colnames(coefs_seas_no_am_tab) <- c("Variables","Freq_No_Am")
colnames(coefs_seas_asia_tab) <- c("Variables","Freq_Asia")
colnames(coefs_seas_vec_tab) <- c("Variables","Freq_All")

# Merge to one file
coefs_all_cont <- merge(coefs_seas_afr_tab,coefs_seas_eur_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_no_am_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_asia_tab,all=T,by="Variables")
coefs_all_cont <- merge(coefs_all_cont,coefs_seas_vec_tab,all=T,by="Variables")
coefs_all_cont <- coefs_all_cont[order(coefs_all_cont$Freq_All,decreasing=F),]

colName <- colnames(coefs_all_cont)
coefs_all_cont_889 <- cbind(coefs_all_cont[,1],coefs_all_cont[,2:6]/pixel*100) # calculate percentage of all grid points
coefs_all_cont <- data.frame(coefs_all_cont[,1],coefs_all_cont[,2]/sum(!is.na(loc_afr_pixels[[1]]))*100,coefs_all_cont[,3]/sum(!is.na(loc_eur_pixels[[1]]))*100,
                             coefs_all_cont[,4]/sum(!is.na(loc_no_am_pixels[[1]]))*100,coefs_all_cont[,5]/sum(!is.na(loc_asia_pixels[[1]]))*100, coefs_all_cont[,6]/pixel*100) # calculate percentage of all grid points according to the respective number of grid points at each continent
colnames(coefs_all_cont) <- colName

coefs_all_cont_mat <- as.matrix(coefs_all_cont_889)
row.names(coefs_all_cont_mat) <- coefs_all_cont_mat[,1]
coefs_all_cont_mat[which(is.na(coefs_all_cont_mat))] <- 0
coefs_all_cont_mat <- coefs_all_cont_mat[,-c(1,6)]

# 4 Plots in one line (all, N.America, Europe, Asia)
line2user <- function(line, side) { # from https://stackoverflow.com/questions/14660372/common-main-title-of-a-figure-panel-compiled-with-parmfrow
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

pdf(file="barplot_variables_lambda1se_separate_continents.pdf", width=16,height=12)
par(mfrow=c(1,4),mar=c(c(5, 4, 1, 0.5)),oma=c(0,7.5,2,0))
barplot(t(coefs_all_cont_mat),horiz=T,las=1,col=c("brown3","DarkOrange2","goldenrod3","burlywood1"),
        xlab="",cex.names=1.8,font=1,cex=1.7,
        legend.text=c("Africa","Europe","North America","Asia"),args.legend =list(x= "bottomright",cex=1.5,text.font=1),main="")
mtext("All continents",side=3,font=2,line=-1.4, cex=1.4)
mtext("(a)",side=3,font=1,line=0.6, adj=0.01,cex=1.6)
barplot(coefs_all_cont$Freq_No_Am,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6,font=1,cex=1.7)
mtext("North America",side=3,font=2,line=-1.4, cex=1.4)
mtext("(b)",side=3,font=1,line=0.6, adj=0.01,cex=1.6)
barplot(coefs_all_cont$Freq_Eur,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6,font=1,cex=1.7)
mtext("Europe",side=3,font=2,line=-1.4, cex=1.4)
mtext("(c)",side=3,font=1,line=0.6, adj=0.01,cex=1.6)
text(line2user(line=mean(par('mar')[c(2, 3)]), side=2), 
     line2user(line=3.5, side=1), 'Percentage of grid points, for which the predictor is included in the regression model', xpd=NA, cex=2.2)
barplot(coefs_all_cont$Freq_Asia,horiz=T,las=1,col="LightCyan3",main="",
        xlab="",cex.names=0.6,font=1,cex=1.7)
mtext("Asia",side=3,font=2,line=-1.4, cex=1.4)
mtext("(d)",side=3,font=1,line=0.6, adj=0.01,cex=1.6)
dev.off()



# Figure A3: Number of months in the growing season ####
########################################################

nb_month_GS <- integer(length = final_pix_num) # Number of months in the growing season

if (climate_crop_provided) {
  for (pixel in 1:final_pix_num) {
    pix_in_995 <- final_pixels_coord$ref_in_995[pixel]
    nb_month_GS[pixel] <- sum(substr(colnames(x1_train_list[[pix_in_995]]), start = 1, stop = 3)=="pr_")
  }#end for pixel
} else {
  load(paste0(path_data,"/nb_months_GS_vector.Rdata"))
} #end if else crop yield provided

levels_nb_month <- cut(nb_month_GS, breaks = c(5,8,11,14), right = F)
levels_nb_month <- gsub(","," - ",levels_nb_month,fixed=TRUE)
DF_GSmonth <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], levels_nb_month = levels_nb_month)
DF_GSmonth$levels_nb_month <- gsub("\\[|\\)","",levels_nb_month)

ggplot(data = DF_GSmonth, aes(x=DF_GSmonth$lon, y=DF_GSmonth$lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF_GSmonth$levels_nb_month)) +
  scale_fill_manual(values=c("5 - 8" = "#edf8b1", "8 - 11" = "#41b6c4", "11 - 14" = "#2c7fb8" ),
                    breaks=c("5 - 8", "8 - 11", "11 - 14"),
                    label=c("5 - 8"="5 - 7", "8 - 11"="8 - 10", "11 - 14"="11 - 13")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(DF_GSmonth$lat)-1, max(DF_GSmonth$lat+1)),
              ratio = 1)+
  labs(fill="Growing season\nlength (months)"
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave(filename = "nb_months_studied_3classes.png", width = 20, height = 3.5)



# Figure A4: Selected climate extreme indicators ####
#####################################################

# 4 climate extreme indicators dtr, frs, TX90p and Rx5day
coeff_dtr <- numeric()
coeff_frs <- numeric()
coeff_tx90p <- numeric()
coeff_rx5 <- numeric()

for (pix in 1:length(coeff)) {
  coeff_dtr[pix] <- coeff[[pix]]["dtr",]
  coeff_frs[pix] <- coeff[[pix]]["frs",]
  coeff_tx90p[pix] <- coeff[[pix]]["tx90p",]
  coeff_rx5[pix] <- coeff[[pix]]["rx5",]
}#end for pix


# Plot grid points where coefficient dtr, frs, TX90p and Rx5day are selected
# dtr
DF1_dtr <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], selected_coeff = (coeff_dtr!=0))
DF1_dtr$selected_coeff <- as.factor(DF1_dtr$selected_coeff)

P1 <- ggplot(data = DF1_dtr, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF1_dtr$selected_coeff)) +
  scale_fill_manual(values = c("FALSE"="gray", "TRUE"="#984ea3"),
                    labels=c("FALSE"="No", "TRUE"="Yes"),
                    breaks=c("TRUE", "FALSE")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Selection of\ndtr coeff."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))

# frs
DF2_frs <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], selected_coeff = (coeff_frs!=0))
DF2_frs$selected_coeff <- as.factor(DF2_frs$selected_coeff)

P2 <- ggplot(data = DF2_frs, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF2_frs$selected_coeff)) +
  scale_fill_manual(values = c("FALSE"="gray", "TRUE"="#984ea3"),
                    labels=c("FALSE"="No", "TRUE"="Yes"),
                    breaks=c("TRUE", "FALSE")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Selection of\nfrs coeff."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))

# rx5
DF3_rx5 <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], selected_coeff = (coeff_rx5!=0))
DF3_rx5$selected_coeff <- as.factor(DF3_rx5$selected_coeff)

P3 <- ggplot(data = DF3_rx5, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF3_rx5$selected_coeff)) +
  scale_fill_manual(values = c("FALSE"="gray", "TRUE"="#984ea3"),
                    labels=c("FALSE"="No", "TRUE"="Yes"),
                    breaks=c("TRUE", "FALSE")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Selection of\nRx5day coeff."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))

# tx90p
DF4_tx90p <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], selected_coeff = (coeff_tx90p!=0))
DF4_tx90p$selected_coeff <- as.factor(DF4_tx90p$selected_coeff)

P4 <- ggplot(data = DF4_tx90p, aes(x=lon, y=lat)) +
  geom_polygon(data = world, aes(long, lat, group=group),
               fill="white", color="black", size=0.3) +
  geom_tile(aes(fill=DF4_tx90p$selected_coeff)) +
  scale_fill_manual(values = c("FALSE"="gray", "TRUE"="#984ea3"),
                    labels=c("FALSE"="No", "TRUE"="Yes"),
                    breaks=c("TRUE", "FALSE")) +
  theme(panel.ontop = F, panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
  # ylab(expression("Lat " ( degree*N))) +
  # xlab(expression("Lon " ( degree*E))) +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
              ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
              ratio = 1)+
  labs(fill="Selection of\nTX90p coeff."
  )+
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 14))

# Combine 4 subplots to one plot
L1 <- get_legend(P1+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L2 <- get_legend(P2+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L3 <- get_legend(P3+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))
L4 <- get_legend(P4+ theme(legend.title = element_text(size=25),
                           legend.text = element_text(size=25),
                           legend.key.size = unit(1.8,"line")))

ggarrange(P1 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L1,
          P2 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L2,
          P3 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L3,
          P4 + theme(legend.position = "none",
                     axis.title.x=element_blank(),axis.title.y=element_blank(),
                     axis.text.x = element_text(size = 20),
                     axis.text.y = element_text(size = 20)),
          L4,
          nrow = 4, ncol=2,labels = c("(a)", "", "(b)", "",
                                      "(c)", "", "(d)", ""),
          label.x = -0.015
          ,widths=c(8,1), heights=c(1,1,1), font.label = list(size = 14, face = "plain", color ="black")
)

ggsave(filename = "selection_of_coeff-dtr-frs-rx5-tx90.png", width = 30, height = 20)



# GIFs: Grid points with inclusion of monthly predictors ####
#############################################################

message("run Fig. 5 for this section first")
allvariables_adj <- allvariables # Adjust names
allvariables_adj <- gsub(x=allvariables_adj, pattern="vpd", replacement = "VPD")
allvariables_adj <- gsub(x=allvariables_adj, pattern="tmax", replacement = "Tmax")
allvariables_adj <- gsub(x=allvariables_adj, pattern="pr", replacement = "Pr")
allvariables_adj <- gsub(x=allvariables_adj, pattern="APr", replacement = "Apr") # recorrect April

for (varia in 1:length(allvariables)) {
  varia_name <- allvariables[varia]
  varia_in_pix <- numeric()
  plots <- vector("list",length=(length(allvariables)))
  for (pix in 1:final_pix_num) {
    varia_in_pix[pix] <- (varia_name %in% row.names(coeff[[pix]])[which(coeff[[pix]]!=0)])
  }
  varia_name <- allvariables_adj[varia]
  
  DF_var <- data.frame(lon=coord_subset[,1], lat = coord_subset[,2], var_in = varia_in_pix)
  DF_var$var_in <- as.factor(DF_var$var_in)
  
  ggplot(data = DF_var, aes(x=lon, y=lat)) +
    geom_polygon(data = world, aes(long, lat, group=group),
                 fill="white", color="black", size=0.3) +
    geom_tile(aes(fill=DF_var$var_in)) +
    scale_fill_manual(drop=F,values = c("1"="#984ea3", "0"="gray"),
                      label= c("1"="Yes", "0"="No"),
                      breaks=c("1","0"), limits=c(0,1)) +
    theme(panel.ontop = F, panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15))+
    # ylab(expression("Lat " ( degree*N))) +
    # xlab(expression("Lon " ( degree*E))) +
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    coord_fixed(xlim = c(min(coord_subset[,1])-1, max(coord_subset[,1]+1)),
                ylim = c(min(coord_subset[,2])-1, max(coord_subset[,2]+1)),
                ratio = 1)+
    labs(fill="Selection",title = varia_name)+
    theme(plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15),
          legend.text = element_text(size = 14))
  ggsave(filename=paste0(varia_name,".jpg"), width = 20, height = 6)
}

# save the variables required for running this file witout crop and climate data
# save(global_meteovar_correlations, global_xtrm_correlations, france_meteovar_correlations, france_xtrm_correlations, csi, nb_month_GS, numLevels_list, allvariables, file="Required_variables.RData") 
