# MAF-NMFI
Models for estimating malaria-attributable fever and non-malarial febrile illness from household survey data

### R script to run NMFI model with treatment propensity

setwd("E:/Users/zool1232.NDPH/Work/NMFI/")
rm(list=ls())

nmfi.data <- read.csv("Data/cellnum.data_1108_all_withcovariates.csv")
ncoarse <- nmfi.data$cell500

long.lat <- cbind(nmfi.data$long,nmfi.data$lat)
year <- nmfi.data$year_start
inholdout <- nmfi.data$in.holdout

##### Fit model to partial dataset (i.e., removing holdouts from likelihood function) to quantify predictive capacity

# Build INLA mesh
library(raster)
Pflimits <- raster("Z:\\cubes\\5km\\Pf_limits\\Pf_limits.tif")  # load spatial limits raster
Pflimits <- rasterToPoints(Pflimits) # point values for raster instead of cellnumber values
Pflimits <- Pflimits[Pflimits[,3] %in% c(0,2),1:2] # Setting boundaries for mesh
Pflimits.mainland <- Pflimits[Pflimits[,2]>Pflimits[,1]*2-102,]
Pflimits.Mad <- Pflimits[Pflimits[,2]<Pflimits[,1]*2-102,]
library(INLA)
library(splancs)
bdy.mainland <- inla.nonconvex.hull(as.matrix(Pflimits.mainland),convex=-0.05,resolution=100)
bdy.Mad <- inla.nonconvex.hull(as.matrix(Pflimits.Mad),convex=-0.05,resolution=100)
Africa.mesh <- inla.mesh.2d(rbind(long.lat,bdy.mainland$loc,bdy.Mad$loc),max.edge=4,cut=1.5,offset=c(5,5))

nmfi.data$pbg.predicted <- nmfi.data$fever.onemodel.yearly
nmfi.data$pbg.predicted[is.na(nmfi.data$pbg.predicted)] <- mean(nmfi.data$pbg.predicted,na.rm=T)
