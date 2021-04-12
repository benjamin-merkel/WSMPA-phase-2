#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages

pckg <- c('raster','fasterize','ncdf4','sf','sp','rgdal','rgeos',
          'nngeo','stringr','readxl','classInt','smoothr','usdm',
          'biomod2','ecospat',"e1071",
          'viridis','RColorBrewer',"scales","remote","angstroms","data.table",
          'adehabitatHR','concaveman','corrplot','RStoolbox',"PresenceAbsence")

for(i in 1:length(pckg)) do.call("library", list(pckg[i]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define directories

maud_directory       <- getwd()
gbm_output_directory <- paste0(maud_directory, "/data/gbm")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define projections 

#EPSG:4326 lat lon unprojected
lonlat.proj <- CRS('+proj=longlat +datum=WGS84 +no_defs') 
# EPSG:102021 south pole stereographic -> conformal map projection which has its line of true latitude at the South Pole
south_pole_stereo.proj      <- CRS("+proj=stere +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 
# EPSG:102020 South Pole Lambert Azimuthal Equal Area
south_pole_equal_area.proj  <- CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load common map files 

grat.30     <- st_read("map data/ne_50m_graticules_30.shp")
grat.10     <- st_read("map data/ne_50m_graticules_10.shp")
grat.5      <- st_read("map data/ne_50m_graticules_5.shp")
land        <- st_read("map data/ne_10m_land.shp")
ice_shelf   <- st_read("map data/Ice_shelf.shp")
ice_shelf   <- st_buffer(ice_shelf,dist = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WSMPAP2 custom extents 

dml.ext            <- extent(0, 30, -80, -60)
dml.wider.ext      <- extent(-10, 40, -80, -55)
extd         <- extent(-186577.5,1629845,1575388,3326885) #extent(-186577.5,1599845,1575388,3326885)

pts = matrix(c(0  ,seq(0,30,0.1),30 ,0,
               -80,rep(-60,301) ,-80,-80),
             ncol=2, byrow=F)
wsmpap2.box = st_as_sf(spPolygons((pts)))
st_crs(wsmpap2.box) <- lonlat.proj
wsmpap2.box <- st_difference(wsmpap2.box,land)
wsmpap2.box <- nngeo::st_remove_holes(wsmpap2.box)

wsmpap2.border <- st_cast(wsmpap2.box, "LINESTRING")
wsmpap2.border <- st_transform(wsmpap2.border, st_crs(ice_shelf))
wsmpap2.border <- st_difference(wsmpap2.border, ice_shelf)
wsmpap2.border <- st_cast(wsmpap2.border, "LINESTRING")
wsmpap2.border <- wsmpap2.border[st_length(wsmpap2.border) == max(st_length(wsmpap2.border)),]


pts = matrix(c(-10  ,seq(-10,40,0.1),40 ,-10,
               -80,rep(-50,501) ,-80,-80),
             ncol=2, byrow=F)
wider.wsmpap2 = st_as_sf(spPolygons((pts)))
st_crs(wider.wsmpap2) <- lonlat.proj
wider.wsmpap2 <- st_difference(wider.wsmpap2,land)
wider.wsmpap2 <- nngeo::st_remove_holes(wider.wsmpap2)

wider.wsmpap2.border <- st_cast(wider.wsmpap2, "POLYGON")
wider.wsmpap2.border <- st_cast(wider.wsmpap2.border[st_area(wider.wsmpap2.border) == max(st_area(wider.wsmpap2.border)),], "LINESTRING")
wider.wsmpap2.border <- st_transform(wider.wsmpap2.border, st_crs(ice_shelf))
wider.wsmpap2.border <- st_difference(wider.wsmpap2.border, ice_shelf)
wider.wsmpap2.border <- st_cast(wider.wsmpap2.border, "LINESTRING")
wider.wsmpap2.border <- wider.wsmpap2.border[st_length(wider.wsmpap2.border) == max(st_length(wider.wsmpap2.border)),]


pts = matrix(c(-90  ,seq(-90,90,0.1),90 ,-90,
               -89,rep(-50,1801) ,-89,-89),
             ncol=2, byrow=F)
ext90E90W50S = st_as_sf(spPolygons((pts)))
st_crs(ext90E90W50S) <- lonlat.proj
# ext90E90W50S <- st_difference(ext90E90W50S,land)
# ext90E90W50S <- nngeo::st_remove_holes(ext90E90W50S)

ext90E90W50S.border <- st_cast(ext90E90W50S, "LINESTRING")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reproject 

wsmpap2.box   <- st_transform(wsmpap2.box,south_pole_equal_area.proj)
wsmpap2.border <- st_transform(wsmpap2.border,south_pole_equal_area.proj)
wider.wsmpap2 <- st_transform(wider.wsmpap2,south_pole_equal_area.proj)
wider.wsmpap2.border <- st_transform(wider.wsmpap2.border,south_pole_equal_area.proj)
ext90E90W50S  <- st_transform(ext90E90W50S,south_pole_equal_area.proj)
ext90E90W50S.border  <- st_transform(ext90E90W50S.border,south_pole_equal_area.proj)
land50        <- st_transform(st_crop(land, extent(-180,180,-90,-50)),south_pole_equal_area.proj)
land          <- st_transform(land,south_pole_equal_area.proj)

ice_shelf   <- st_transform(ice_shelf,south_pole_equal_area.proj)
grat.30     <- st_transform(grat.30,south_pole_equal_area.proj)
grat.10     <- st_transform(st_crop(grat.10, extent(-180,180,-85,0)),south_pole_equal_area.proj)
grat.5      <- st_transform(grat.5, south_pole_equal_area.proj)



# plot model domains
ratio <- (extent(ext90E90W50S)[2]-extent(ext90E90W50S)[1])/
  (extent(ext90E90W50S)[4]-extent(ext90E90W50S)[3])
setwd(maud_directory)
png("figures/model domains.png", res = 800, width=20*ratio, height = 20, units="cm")
par(mar=c(1,1,1,1))
plot(st_geometry(ext90E90W50S),border="transparent",col=brewer.pal(3,"Set3")[1])
plot(st_geometry(ext90E90W50S.border),add=T,border="transparent",col=brewer.pal(3,"Set3")[1])
plot(st_geometry(wider.wsmpap2),add=T,border="transparent",col=brewer.pal(3,"Set3")[2])
plot(st_geometry(wsmpap2.box),add=T,border="transparent",col=brewer.pal(3,"Set3")[3])
plot(st_geometry(wider.wsmpap2),lty=2,add=T,lwd=1.5)
plot(st_geometry(ice_shelf),add=T,border=grey(0.9),col=grey(0.9),lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(grat.10),lwd=0.5,col=grey(0.2),add=T)
# plot(st_geometry(ext90E90W50S),lty=2,add=T)
plot(st_geometry(wsmpap2.border),add=T,lwd=1.5)
text(c(-4250,-3190, -2100, -1000),rep(-50,4), cex=0.7,las = 2,
       labels = c(expression(paste("50",degree,"S")), 
                  expression(paste("60",degree,"S")), 
                  expression(paste("70",degree,"S")), 
                  expression(paste("80",degree,"S"))))
axis(3,at = c(-3800,-2650,-1650,-800, 0, 800, 1650, 2650, 3800),  
     tick = F, line = -1, cex.axis=0.7,
     labels = c(expression(paste("40",degree,"W")), 
                expression(paste("30",degree,"W")), 
                expression(paste("20",degree,"W")), 
                expression(paste("10",degree,"W")), 
                expression(paste("0",degree,"")), 
                expression(paste("10",degree,"E")),
                expression(paste("20",degree,"E")),
                expression(paste("30",degree,"E")),
                expression(paste("40",degree,"E"))))
box()
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create polygon for pseudo absences 

grat.50S             <- grat.5[grat.5$direction=="S" & grat.5$degrees==50,][,1]
grat.50S             <- rbind(grat.50S, grat.50S[1,])
grat.50S             <- concaveman(grat.50S) # works
grat.50S             <- st_buffer(grat.50S,dist = 0)
background_polygon <- st_difference((grat.50S),st_geometry(ice_shelf))
background_polygon <- st_difference(background_polygon,st_geometry(land50))
background_raster  <- raster(as_Spatial(background_polygon),res=10)
values(background_raster) <- 1:ncell(background_raster)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load environmental covariates

# Initial selection of relevant covariates
# World ocean atlas temperature (surface and 200m depth): WOA_temp_0, WOA_temp_200
# World ocean atlas salinity (surface and 200m depth): WOA_sal_0, WOA_sal_200
# World ocean atlas silicate and nitrate (surface): WOA_si_0, WOA_ni_0
# World ocean atlas dissolved oxygen concentration (surface and 200m depth): WOA_ox_0, WOA_ox_200
# Sea ice metrics (length of the sea ice season, timing of sea ice retreat, summer sea ice concentration): NSIDC_ice_duration, NSIDC_ice_retreat, NSIDC_ice_conc
# Summer bloom duration (bd), initiation timing (BI), timing of peak concentration (bt) and average Chl a concentration during the bloom (bm) derived from OC-CCI: OCCCI_bd, OCCCI_bi, OCCCI_bt, OCCCI_bm
# Bathymetry and its derivative distance to to the continental shelf (defined as 1000 m depth): bath, dis_1000
# Mixed layer depth from Pellichero et al. 2019: Pellichero_ml_depth

setwd(maud_directory)
env.init.selected.names <- names(readRDS("data/WSMPA P2 Initial environmental covariate raster stack.rds"))
env.init.selected <- stack("data/WSMPA P2 Initial environmental covariate raster stack.tif")
names(env.init.selected) <- env.init.selected.names
# env.init.selected <- readRDS("data/WSMPA P2 Initial environmental covariate raster stack.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define species

species     <- "P.antarctica"
species     <- "E.superba"
species     <- "E.crystallorophias"


if(species == "E.superba") {
  model.domain <- wider.wsmpap2
  model.domain.border <- wider.wsmpap2.border
  range <- "wider WSMPAP2"
} else {
  model.domain <- ext90E90W50S
  model.domain.border <- ext90E90W50S.border
  range <- "90W-90E"
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# crop covariates for model domain

env.cor.test <- crop(mask(env.init.selected,wider.wsmpap2),model.domain)

M  <- cor(as.data.frame(env.cor.test), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))
hc <- hclust(as.dist(1-abs(M)), method = "complete")

opar <- par(mfrow=c(1,1),mar=c(2,4,2,0))
plot(hc,hang=-1,las=1,main='',xlab="",ylab="pearson correlation")
rect.hclust(hc,h=0.3)
par(opar)


# test colinearity of covariates in model domain 
vif2.select <- vifcor(env.cor.test, th = 0.7)
vif2.select

# remove colinear parameters at cutoff th = 0.7
env.selected <- env.init.selected[[vif2.select@results$Variables]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load response data

setwd(maud_directory)
data <- readRDS(paste0("data/WSMPA P2_",species,"_response data.rds"))

# extract environmental covariates for each observation and remove locations with NA values
data <- cbind(data, extract(env.selected, data))
data <- data[complete.cases(data.frame(data[,names(env.selected)])),]

presence              <- data[data$presence_absence %in% 1,]
absence               <- data[data$presence_absence %in% 0,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate pseudo-absences using one class support vector machines (OSCVM)

# calculate how many pseudo-absence locations are needed
pseudo.size <- nrow(presence)-nrow(absence)
if(pseudo.size<0) pseudo.size <- 0

if(pseudo.size>0){
  pseudo.start                  <- data.frame(rasterToPoints(crop(background_raster,model.domain)))
  coordinates(pseudo.start)     <- cbind(pseudo.start$x,pseudo.start$y)
  proj4string(pseudo.start)     <- projection(background_raster)
  pseudo.start                  <- SpatialPointsDataFrame(pseudo.start, data.frame(background.bin=extract(background_raster, pseudo.start)))
  pseudo.start                  <- cbind(pseudo.start,extract(env.init.selected, pseudo.start))
  pseudo.start                  <- pseudo.start[!duplicated(pseudo.start$background.bin),]
  pseudo.start$presence_absence <- rep(NA,nrow(pseudo.start))
  
  
  # create training dataset for OSCVM
  selected.variables <- names(env.selected)
  train <- data.frame(rbind(presence2,absence2))
  train$presence <- TRUE
  train$presence[train$presence_absence==0] <- FALSE
  train <- train[train$presence==T,]
  trainLabels<-train[,"presence"]
  trainpredictors<-train[,selected.variables]
  
  # determine nu and gamma to be used
  tuned <- tune.svm(x = trainpredictors, 
                    y = trainLabels,
                    nu =  seq(1,50,1)/10000,
                    gamma = 0.01,#10^(-2:0), 
                    type='one-classification'                 
  )
  tuned
  
  # run OCSVM with calibrated nu and gamma
  svm.model <-svm(x = trainpredictors, 
                  y = trainLabels,
                  type='one-classification', 
                  nu = tuned$best.parameters$nu, 
                  gamma = tuned$best.parameters$gamma
  )
  
  # determine which background locations are outside the hypervolume created
  ps2 <- pseudo.start[,selected.variables]
  ps2 <- ps2[complete.cases(data.frame(ps2)),]
  testpredictors  <- data.frame(ps2)[,selected.variables]
  
  svm.pred <- predict(object  = svm.model, 
                      newdata = testpredictors)
  ps2$svm.pred <- svm.pred
  
  ps3 <- spTransform(ps2[ps2$svm.pred==F,],proj4string(presence))
  ps3$presence_absence <- NA
  ps3 <- ps3[,"presence_absence"]
  ps3 <- cbind(ps3, extract(env.selected, ps3))
  ps3 <- ps3[complete.cases(data.frame(ps3[,names(env.selected)])),]
  ps3 <- as_Spatial(st_intersection(st_as_sf(ps3), model.domain))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODELLING

# define model settings
myBiomodOption <- BIOMOD_ModelingOptions(GBM=list(n.trees = 10000, 
                                                  cv.folds=10,
                                                  n.cores=3,
                                                  interaction.depth = 5,
                                                  shrinkage = 0.01))

# run procedure for 10 randomly determined sets of pseudo-absence locations based on OCSVM
biomod_gbm <- gbm_Eval <- dat.used <- NULL

nruns = 10
random.runs <- 10
if(pseudo.size==0) random.runs <- 1

for(randomrun in 1:random.runs){
  
  if(random.runs>1) {
    dat.used[[randomrun]] <- rbind(presence2[,"presence_absence"], 
                    ps3[sample(1:nrow(ps3),pseudo.size),"presence_absence"],
                    absence2[,"presence_absence"])
  } else {
    dat.used[[randomrun]] <- rbind(presence2[,"presence_absence"], 
                                   absence2[,"presence_absence"])
  }
  names(dat.used[[randomrun]]) <- paste(substr(species,1,3),randomrun,sep="_")
  
  biomod_data <- BIOMOD_FormatingData(resp.var = dat.used[[randomrun]], 
                                      resp.name= paste(substr(species,1,3),randomrun,sep="_"),
                                      expl.var = env.selected)
  
  setwd(gbm_output_directory)
  biomod_gbm[[randomrun]] <- BIOMOD_Modeling(data              = biomod_data,
                                             models            = "GBM",
                                             models.options    = myBiomodOption,
                                             NbRunEval         = nruns,
                                             DataSplit         = 70,
                                             Prevalence        = NULL,
                                             Yweights          = NULL,
                                             VarImport         = 10,
                                             models.eval.meth  = c('TSS','ROC',"ACCURACY","KAPPA"),
                                             SaveObj           = TRUE,
                                             rescal.all.models = FALSE,
                                             do.full.models    = T,
                                             modeling.id       = paste('GBM',range,randomrun)) 
  
  
  gbm_Eval[[randomrun]] <- get_evaluations(biomod_gbm[[randomrun]])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load models

biomod_gbm <- NULL
setwd(gbm_output_directory)
temp.path <- list.files(pattern=".models.out", recursive = TRUE)
temp.path <- temp.path[grepl(substr(species,1,3), temp.path)][1:10]
if(exists("myBiomodOuttmp")) { rm(myBiomodOuttmp) }
for(i in 1:random.runs)  biomod_gbm[[i]] <- get(load(temp.path[i]))

gbm_Eval <- NULL
for(randomrun in 1:random.runs){
  # biomod_gbm[[randomrun]] <- BIOMOD_LoadModels(biomod_gbm[[randomrun]])
  gbm_Eval[[randomrun]] <- get_evaluations(biomod_gbm[[randomrun]])
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check model evaluations

for(randomrun in 1:random.runs){
  tss1 <- gbm_Eval[[randomrun]]["TSS","Testing.data",,1:nruns,]
  if(randomrun==1) tss2 <- tss1 else tss2 <- c(tss2, tss1)
  roc1 <- gbm_Eval[[randomrun]]["ROC","Testing.data",,1:nruns,]
  if(randomrun==1) roc2 <- roc1 else roc2 <- c(roc2, roc1)
}
summary(roc2)
summary(tss2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check marginal response plots for each covariate across models

setwd(gbm_output_directory)
for(randomrun in 1:random.runs){
  resp2 <- response.plot2(
    models = BIOMOD_LoadModels(biomod_gbm[[randomrun]]),
    Data = get_formal_data(biomod_gbm[[randomrun]], 'expl.var'),
    show.variables = get_formal_data(biomod_gbm[[randomrun]],'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var.metric = 'mean',
    plot=F,
    legend = F,
    data_species = get_formal_data(biomod_gbm[[randomrun]], 'resp.var')
  )
  if(randomrun==1) Biomodresponse <- resp2 else Biomodresponse <- rbind(Biomodresponse, resp2) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check variable importance across all models

for(randomrun in 1:random.runs){
  var_imp <- get_variables_importance(biomod_gbm[[randomrun]])
  # var_imp[,,,"AllData"]
  vi3b <- t(var_imp[,,1:nruns,"AllData"])
  if(randomrun==1) vi3 <- vi3b else vi3 <- rbind(vi3b, vi3)
}
vi4 <- apply(vi3,2,mean)


setwd(maud_directory)
png(paste0("figures/",species ," variable importance.png"), res = 800, width=17, height = 17, units="cm")
opar <- par(mfrow=c(1,1),mar=c(4,10,1,1))
boxplot(vi3[,order(vi4)], horizontal=T,las=1,
        col = grey(1), border = "skyblue3", lty = 1, lwd = 2,
        ylab = "", xlab = "variable importance")
points(vi4[order(vi4)],1:length(vi4),cex=1.5,pch=19,col="darkred")
par(opar)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict model output over model domain

setwd(gbm_output_directory)
new.env <- stack(crop(env.selected,model.domain))
for(randomrun in 1:random.runs){
  proj1 <- BIOMOD_Projection(modeling.output = biomod_gbm[[randomrun]],
                             new.env = new.env,
                             proj.name ='current',
                             selected.models = "all",
                             binary.meth =NULL,
                             compress =F,
                             clamping.mask = F,
                             output.format ='.grd')
  if(randomrun==1) biomod_gbm_Proj <- proj1@proj@val else biomod_gbm_Proj <- stack(biomod_gbm_Proj, proj1@proj@val)
}

# mean model prediction based on 10 models using all available data
gbm_ras_full <- mean(raster::subset(biomod_gbm_Proj, grep('Full', names(biomod_gbm_Proj), value = T)))

# mean model prediction based on 100 models using 70% of available data
gbm_ras_mean <- mean(raster::subset(biomod_gbm_Proj, grep('RUN', names(biomod_gbm_Proj), value = T)))

# standard deviation of model prediction based on 100 models using 70% of available data
gbm_ras_sd   <- calc(raster::subset(biomod_gbm_Proj, grep('RUN', names(biomod_gbm_Proj), value = T)), sd)

# plot model predictions
opar <- par(mfrow=c(2,2))
plot(mask(gbm_ras_mean,model.domain), main="Full model", 
     col = rev(brewer.pal(10,"Spectral")), breaks=seq(0,1000,100))
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(wsmpap2.box),add=T)

plot(mask(gbm_ras_full,model.domain), main="Mean model", 
     col = rev(brewer.pal(10,"Spectral")), breaks=seq(0,1000,100))
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(wsmpap2.box),add=T)

plot(mask(gbm_ras_sd,model.domain),   main="model SD",   
     col = rev(plasma(10)), breaks=seq(0,300,30))
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(wsmpap2.box),add=T)
par(opar)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE

setwd(maud_directory)
writeRaster(gbm_ras_mean, filename=paste0("data/",species,"_",range,"_GBM_MEAN_100_ensemble_with OCSVM.tif"), format="GTiff", overwrite=TRUE)
writeRaster(gbm_ras_full, filename=paste0("data/",species,"_",range,"_GBM_FULL_MEAN_100_ensemble_with OCSVM.tif"), format="GTiff", overwrite=TRUE)
writeRaster(gbm_ras_sd  , filename=paste0("data/",species,"_",range,"_GBM_SD_100_ensemble_with OCSVM.tif"), format="GTiff", overwrite=TRUE)

saveRDS(dat.used         , paste0("data/",species,"_",range,"_response_data.rds"))
saveRDS(Biomodresponse   , paste0("data/",species,"_",range,"_GBM_response_curve_data.rds"))
saveRDS(vi3              , paste0("data/",species,"_",range,"_GBM_estimated variable importance.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting

# load rasters
gbm_ras_mean <- raster(paste0("data/",species,"_",range,"_GBM_MEAN_100_ensemble_with OCSVM.tif"))
gbm_ras_sd   <- raster(paste0("data/",species,"_",range,"_GBM_SD_100_ensemble_with OCSVM.tif"))
dat.used     <- readRDS(paste0("data/",species,"_",range,"_response_data.rds"))


# plot available data distribution
ratio <- (extent(model.domain)[2]-extent(model.domain)[1])/
  (extent(model.domain)[4]-extent(model.domain)[3])
setwd(maud_directory)
png(paste0("figures/",species ," available data.png"), res = 800, 
    width=20*ratio, height = 20, units="cm")
par(mar=c(1,1,1,1))
plot(st_geometry(model.domain),lty=2,lwd=1,border="transparent")
plot(st_geometry(model.domain.border),lty=2,add=T,lwd=2)
plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
plot(st_geometry(grat.10),lwd=0.5,col=grey(0.8),add=T)
plot(st_geometry(wsmpap2.border),add=T,lwd=1.5)
if(pseudo.size>0) {for(i in 1:10) plot(dat.used[[i]],add=T,pch=21,bg=5,col="white",lwd=0.5,cex=0.5)}
plot(data,add=T,pch=21,bg=4,col="white",lwd=0.5)
plot(data[data$presence_absence==1,],add=T,pch=21,bg=2,col="white",lwd=0.5)
if(species=='E.superba') {
  axis(2,at = c(4280, 3160, 2020), tick = F, line = -1, cex.axis=0.7,
     labels = c(expression(paste("50",degree,"S")), 
                expression(paste("60",degree,"S")), 
                expression(paste("70",degree,"S"))))
} else {
  text(c(-4250,-3190, -2100, -1000),rep(-50,4), cex=0.7,las = 2,
       labels = c(expression(paste("50",degree,"S")), 
                  expression(paste("60",degree,"S")), 
                  expression(paste("70",degree,"S")), 
                  expression(paste("80",degree,"S"))))
}
axis(3,at = c(-3800,-2650,-1650,-800, 0, 800, 1650, 2650, 3800),  
     tick = F, line = -1, cex.axis=0.7,
     labels = c(expression(paste("40",degree,"W")), 
                expression(paste("30",degree,"W")), 
                expression(paste("20",degree,"W")), 
                expression(paste("10",degree,"W")), 
                expression(paste("0",degree,"")), 
                expression(paste("10",degree,"E")),
                expression(paste("20",degree,"E")),
                expression(paste("30",degree,"E")),
                expression(paste("40",degree,"E"))))
box()
dev.off()


# plot mean and sd of model output over model domain
plot.sd=T
for(plot.sd in c(T, F)){

  if(plot.sd==T) layer1 <- gbm_ras_sd else  layer1 <- gbm_ras_mean 
  if(plot.sd==T) brks = seq(0,ceiling(max(values(gbm_ras_sd),na.rm=T)),ceiling(max(values(gbm_ras_sd),na.rm=T))/10) else brks = seq(0,1000,100) 
  if(plot.sd==T) legend.axis.labels = c("less","more") else legend.axis.labels = c("low","high")
  if(plot.sd==T) legend.labels = "uncertainty" else legend.labels = "probability"
  vid <- brewer.pal(11,"Spectral")
  cols_used <- colorRampPalette(vid[c(11,6,1)])(10)
  if(plot.sd==T) cols_used <- plasma(10)
  ratio <- (extent(wider.wsmpap2)[2]-extent(wider.wsmpap2)[1])/
    (extent(wider.wsmpap2)[4]-extent(wider.wsmpap2)[3])
  setwd(maud_directory)
  if(plot.sd==F) png(paste0("figures/",species ," model prediction mean map.png"), res = 800,   width=20*ratio, height = 20, units="cm")
  if(plot.sd==T)  png(paste0("figures/",species ," model prediction SD map.png"), res = 800,   width=20*ratio, height = 20, units="cm")
  par(mar=c(1,1,1,1))
  plot(st_geometry(wider.wsmpap2),lty=2)  
  plot(mask(layer1, wider.wsmpap2),  col = cols_used, breaks=brks, add=T, legend =F)
  plot(st_geometry(ice_shelf),add=T,border="grey",col="grey",lwd=0.1)
  plot(st_geometry(land),add=T,border=grey(0.4),col=grey(0.4),lwd=.1)
  plot(st_geometry(grat.10),lwd=0.5,col=grey(0.8),add=T)
  plot(st_geometry(wider.wsmpap2.border),lty=2,add=T,lwd=1.5)
  plot(st_geometry(wsmpap2.border),add=T,lwd=1.5)
  axis(2,at = c(4280, 3160, 2020), tick = F, line = -1, cex.axis=0.7,
       labels = c(expression(paste("50",degree,"S")), 
                  expression(paste("60",degree,"S")), 
                  expression(paste("70",degree,"S"))))
  axis(3,at = c(-800, 0, 800, 1650, 2600),  tick = F, line = -1, cex.axis=0.7,
       labels = c(expression(paste("10",degree,"W")), 
                  expression(paste("0",degree,"")), 
                  expression(paste("10",degree,"E")),
                  expression(paste("20",degree,"E")),
                  expression(paste("30",degree,"E"))))
  box()
  plot(layer1,legend.only = T, col=cols_used, breaks=brks, 
       axis.args=list(at=c(min(brks)+0.04*max(brks),max(brks)-0.04*max(brks)), #brks,
                      labels=legend.axis.labels,#brks/10 ,
                      cex.axis=1,tick=F,line=-0.2),
       legend.args=list(text=legend.labels, side=4, font=2, line=2.5, cex=1),
       legend.width=0.2, legend.shrink=0.75,
       smallplot=c(0.85,0.87, 0.07,0.32))
  dev.off()
}


# plot marginal response curves for each covariate
ncols <- ceiling((nlayers(env.selected)+1)/4)
setwd(maud_directory)
png(paste0("figures/",species ," marginal response curves.png"),
    units="cm",res=500,width=4*ncols, height=15)
par(mfrow = c(4,ncols),mar=c(2,2,2,0.1))
for(i in unique(Biomodresponse$expl.name)){
  xx <- Biomodresponse[Biomodresponse$expl.name==i,]
  xx <- xx[complete.cases(xx),]
  xx <- xx[!grepl("Full", xx$pred.name),]
  # xx$expl.val <- round(xx$expl.val,0)
  response <- data.frame(mean   = tapply(xx$pred.val, xx$expl.val, mean),
                         median = tapply(xx$pred.val, xx$expl.val, median),
                         sd     = tapply(xx$pred.val, xx$expl.val, sd),
                         se     = tapply(xx$pred.val, xx$expl.val, sd)/sqrt(length(unique(xx$pred.name))),
                         q0.025 = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.025)),
                         q0.25  = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.25)),
                         q0.75  = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.75)),
                         q0.975 = tapply(xx$pred.val, xx$expl.val, FUN = function(x) quantile(x,0.975)))
  response$x <- as.numeric(as.character(rownames(response)))
  response$lower.ci <- response$mean - 2 * response$se
  response$upper.ci <- response$mean + 2 * response$se
  
  plot(response$x, response$mean,type="l",lwd=2,ylim=c(0,1),main=i,las=1,ylab="",xlab="",cex.axis=0.8,col="white")
  for(j in unique(Biomodresponse$pred.name)){
    xx2 <- xx[xx$pred.name==j,]
    lines(xx2$expl.val, xx2$pred.val, lwd=1, col=rgb(1,0,0,0.2))
  }
  lines(response$x, frollmean(response$mean, 10), lwd=3, col=1)
  
  rug(biomod_data@data.env.var[,i])
}
plot(1,1,col="white",ann=F,axes=F)
legend("bottomright", legend = c("ensemble mean","single model"),
       lwd=c(4,1),lty=1,col=c(1,2),cex=0.9)

par(opar)
dev.off()



