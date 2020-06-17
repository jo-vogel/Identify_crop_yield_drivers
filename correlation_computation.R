# Compute correlation between crop yield and meteorological variables (monthly data), and extreme indicators
# For one pixel (France) and all pixels
# Authors : Elisabeth Tschumi, Pauline Rivoire

# This code was run under version 3.6 of R

library("sp");library("ncdf4");library("corrplot");library('oce')

monthlymeteovar_xtrm <- readRDS(paste0(path_data,"/extremeindices_and_monthlymeteovar_995pix.Rdata"))


# remove gridpoints that weren't used ####
pixelstokeep <- final_pixels_coord$ref_in_995

unused <- setdiff(1:995,pixelstokeep)

#put data in arrays
monthlymeteovar_array <- array(NA,c(7,995,17,1600)) # extract the meteorological variables
for (varia in 1:7){
  monthlymeteovar_array[varia,,,] <- monthlymeteovar_xtrm[[varia+3]]
}

yield_xtrm_array <- array(NA,c(8,995,1600)) #extract the extreme indicators
for (varia in 1:7){
  yield_xtrm_array[varia+1,,] <- monthlymeteovar_xtrm[[varia+10]]
}
yield_xtrm_array[1,,] <- monthlymeteovar_xtrm[[3]]

#remove unused pixels
for (pix in 1:length(unused)){
  monthlymeteovar_array[,unused[pix],,]<-NA
}

for (pix in 1:length(unused)){
  yield_xtrm_array[,unused[pix],]<-NA
}


# Correlation at every pixel ####

message('gives warnings due to NAs')

global_correlations_meteovars <- array(NA,c(7,995,17))
for (varia in 1:7){ #loop on meteo variabels
  for (pix in 1:995){ #loop on pixels
    for (month in 1:17){ #loop on months
      global_correlations_meteovars[varia,pix,month] <- cor(monthlymeteovar_array[varia,pix,month,], yield_xtrm_array[1,pix,])
    } # end for month
  } #end for pix
} #end for varia
global_correlations_xtrm <- array(NA,c(7,995))
for (varia in 2:8){ #loop on extreme indicators
  for (pix in 1:995){ #loop on pixels
    global_correlations_xtrm[varia-1,pix] <- cor(yield_xtrm_array[varia,pix,], yield_xtrm_array[1,pix,])
  } #end for pix
} #end for varia


# global mean of correlation per season ####

global_seasonal_correlations_meteovars <- array(NA,c(7,995,4))
winter <- apply(monthlymeteovar_array[,,c(5:7,17),],c(1,2,4),mean, na.rm=TRUE)
spring <- apply(monthlymeteovar_array[,,c(8:10),],c(1,2,4),mean, na.rm=TRUE)
summer <- apply(monthlymeteovar_array[,,c(1,11:13),],c(1,2,4),mean, na.rm=TRUE)
autumn <- apply(monthlymeteovar_array[,,c(2:4,14:16),],c(1,2,4),mean, na.rm=TRUE)
for (varia in 1:7){ #loop on variables
  for (pix in 1:995){ #loop on pixels
    for (season in 1:4){ #loop on seasons
      if (season==1){ #winter
        global_seasonal_correlations_meteovars[varia,pix,season] <- cor(winter[varia,pix,], yield_xtrm_array[1,pix,])
      }else if (season==2){ #spring
        global_seasonal_correlations_meteovars[varia,pix,season] <- cor(spring[varia,pix,], yield_xtrm_array[1,pix,])
      }else if (season==3){ #summer
        global_seasonal_correlations_meteovars[varia,pix,season] <- cor(summer[varia,pix,], yield_xtrm_array[1,pix,])
      }else if (season==4){ #autumn
        global_seasonal_correlations_meteovars[varia,pix,season] <- cor(autumn[varia,pix,], yield_xtrm_array[1,pix,])
      }
    } # end for season
  } # end for pix
} #end for varia


# mean per growing season (GS) ####

global_GS_correlations_meteovars <- array(NA,c(7,995))
GS <- apply(monthlymeteovar_array[,,,],c(1,2,4),mean, na.rm=TRUE)
for (varia in 1:7){ # loop on meteo variables
  for (pix in 1:995){ # loop on pixels
    global_GS_correlations_meteovars[varia,pix] <- cor(GS[varia,pix,], yield_xtrm_array[1,pix,])
  } #end for ipix
} #end for varia


# make global means ####
global_means_correlations_meteovars <- apply(global_correlations_meteovars, c(1,3), mean, na.rm=TRUE)
global_means_correlations_xtrm <- apply(global_correlations_xtrm, c(1), mean, na.rm=TRUE)

global_seasonal_means_correlations_meteovars <- apply(global_seasonal_correlations_meteovars, c(1,3), mean, na.rm=TRUE)
global_GS_means_correlations_meteovars <- apply(global_GS_correlations_meteovars, c(1), mean, na.rm=TRUE)


# extract France pixel ####
x <- which.min((monthlymeteovar_xtrm$longitudes-1.1)^2 + (monthlymeteovar_xtrm$latitudes-47.7)^2)

france_correlations_meteovars <- global_correlations_meteovars[,x,]
france_correlations_xtrm <- global_correlations_xtrm[,x]
france_seasonal_correlations_meteovars <- global_seasonal_correlations_meteovars[,x,]
france_GS_correlations_meteovars <- global_GS_correlations_meteovars[,x]


# put all correlation data together ####
global_meteovar_correlations <- array(NA,c(7,22))
global_meteovar_correlations[,1:17]<- global_means_correlations_meteovars
global_meteovar_correlations[,18:21]<- global_seasonal_means_correlations_meteovars
global_meteovar_correlations[,22]<- global_GS_means_correlations_meteovars
rownames(global_meteovar_correlations) <- c(":T[d]","Pr","Wind","Rad",":T[max]",":T[min]","VPD")
colnames(global_meteovar_correlations) <- c("Aug1", "Sep1", "Oct1","Nov1","Dec1","Jan2","Feb2","Mar2","Apr2","May2","Jun2","Jul2","Aug2","Sep2","Oct2","Nov2","Dec2","winter (DJF)","spring (MAM)","summer (JJA)","autumn (SON)","GS    ")

global_xtrm_correlations <- array(NA,c(7,1))
global_xtrm_correlations[,]<- global_means_correlations_xtrm
rownames(global_xtrm_correlations) <- c("dtr","frs","TXx","TNn","Rx5day","TX90p","TN10p")
colnames(global_xtrm_correlations) <- c("GS    ")

france_meteovar_correlations <- array(NA,c(7,22))
france_meteovar_correlations[,1:17]<- france_correlations_meteovars
france_meteovar_correlations[,18:21]<- france_seasonal_correlations_meteovars
france_meteovar_correlations[,22]<- france_GS_correlations_meteovars
rownames(france_meteovar_correlations) <- c(":T[d]","Pr","Wind","Rad",":T[max]",":T[min]","VPD")
colnames(france_meteovar_correlations) <- c("Aug1", "Sep1", "Oct1","Nov1","Dec1","Jan2","Feb2","Mar2","Apr2","May2","Jun2","Jul2","Aug2","Sep2","Oct2","Nov2","Dec2","winter (DJF)","spring (MAM)","summer (JJA)","autumn (SON)","GS    ")

france_xtrm_correlations <- array(NA,c(7,1))
france_xtrm_correlations[,]<- france_correlations_xtrm
rownames(france_xtrm_correlations) <- c("dtr","frs","TXx","TNn","Rx5day","TX90p","TN10p")
colnames(france_xtrm_correlations) <- c("GS    ")