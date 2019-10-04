# Jonathan B. Koch 
# Orginator: Jesse Tabor (UH-Hilo)
# BiodMod2: X. sonorina

## set where the program saves and looks for data ##
setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW")
# www.unil.ch/hsdm
# list files
list.files()

#test

## download the "packages" that have the required functions needed for the study ##
# install.packages("rgbif")
# install.packages("biomod2")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("sf")
# install.packages("rgdal")
# install.packages("ade4")
# install.packages ("maptools")
# instal.packages("raster")
# install.packages("rasterVis")
# install.packages("latticeExtra")
# install.packages("lattice")
# install.packages("sp")
#test

# load libraries
library(rgbif)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(sf)
library(ade4)
library(rgdal)
library(raster)
library(maptools)
library(rasterVis)
library(latticeExtra)
library(lattice)
library(ggmap)
library(maps)
library(mapdata)
library(sp)

## load species data into R
df <- read.table("vosnesenskii_unique_occurrences_filterd_gbif.csv",
                    header=TRUE, 
                    sep=",",
                    row.names=NULL)

# read in shapefile
world <- map_data("world") # from maps
north_america <- subset(world, region %in% c("Canada", "Mexico", "USA"))

# create boundary of data visualization
max_y = max(df$decimalLatitude+1)
min_y = min(df$decimalLatitude-1)
max_x = max(df$decimalLongitude+1)
min_x = min(df$decimalLongitude-1)

# check out data
ggplot() + 
  geom_polygon(data = north_america, aes(x=long, y = lat, group = group)) + 
  coord_sf(xlim = c(min_x, max_x), ylim = c(min_y, max_y), expand = FALSE) + #xlim = longitude, ylim = latitude
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "black", size = 5)

# B. vosnesnskii, remove occurence that don't match Koch et al. 2012 or Williams et al. 2014
nrow(df)
df <- subset(df, decimalLongitude <= -115 )

# check out data
ggplot() + 
  geom_polygon(data = north_america, aes(x=long, y = lat, group = group)) + 
  coord_sf(xlim = c(min_x, max_x), ylim = c(min_y, max_y), expand = FALSE) + #xlim = longitude, ylim = latitude
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "black", size = 5)

# create boundary of data visualization
nrow(df)
max_y = max(df$decimalLatitude+1)
min_y = min(df$decimalLatitude-1)
max_x = max(df$decimalLongitude+1)
min_x = min(df$decimalLongitude-1)

# check out data
ggplot() + 
  geom_polygon(data = north_america, aes(x=long, y = lat, group = group)) + 
  coord_sf(xlim = c(min_x, max_x), ylim = c(min_y, max_y), expand = FALSE) + #xlim = longitude, ylim = latitude
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "black", size = 5)


# create mask from max/min latitude (y) and longitude (x)
# create spatial polygons
# https://rstudio-pubs-static.s3.amazonaws.com/202536_7a122ff56e9f4062b6b012d9921afd80.html

# make bounding box
# to clip bioclimatic rasters
# no need to this is if you did it the first time
# x_coord <- c(max_x, max_x, min_x, min_x)
# y_coord <- c(min_y, max_y, max_y, min_y)
# xym <- cbind(x_coord, y_coord)
# xym

# polygon library
# don't need to make polygon after you
# you clipped the B. vosnesenskii extent
# for bioclimatic variables
# p = Polygon(xym)
# ps = Polygons(list(p),1)
# sps = SpatialPolygons(list(ps))
# plot(sps)

# Enviromental variables selection
# stack of all your climate variables together

# bio <- stack(list.files("bio",pattern = "bio",
#                        full.names = T),
#             RAT = FALSE)

# bio_1 <- bio[[c(1)]]
# bio_1 <- mask(bio_1,sps)
# bio_1 <-crop(bio_1,sps)

# bio_10 <- bio[[c(2)]]
# bio_10 <- mask(bio_10,sps)
# bio_10 <-crop(bio_10,sps)

# bio_11 <- bio[[c(3)]]
# bio_11 <- mask(bio_11,sps)
# bio_11 <-crop(bio_11,sps)

# bio_12 <- bio[[c(4)]]
# bio_12 <- mask(bio_12,sps)
# bio_12 <-crop(bio_12,sps)

# bio_13 <- bio[[c(5)]]
# bio_13 <- mask(bio_13,sps)
# bio_13 <-crop(bio_13,sps)

# bio_14 <- bio[[c(6)]]
# bio_14 <- mask(bio_14,sps)
# bio_14 <-crop(bio_14,sps)

# bio_15 <- bio[[c(7)]]
# bio_15 <- mask(bio_15,sps)
# bio_15 <-crop(bio_15,sps)

# bio_16 <- bio[[c(8)]]
# bio_16 <- mask(bio_16,sps)
# bio_16 <-crop(bio_16,sps)

# bio_17 <- bio[[c(9)]]
# bio_17 <- mask(bio_17,sps)
# bio_17 <-crop(bio_17,sps)

# bio_18 <- bio[[c(10)]]
# bio_18 <- mask(bio_18,sps)
# bio_18 <-crop(bio_18,sps)

# bio_19 <- bio[[c(11)]]
# bio_19 <- mask(bio_19,sps)
# bio_19 <-crop(bio_19,sps)

# bio_2 <- bio[[c(12)]]
# bio_2 <- mask(bio_2,sps)
# bio_2 <-crop(bio_2,sps)

# bio_3 <- bio[[c(13)]]
# bio_3 <- mask(bio_3,sps)
# bio_3 <-crop(bio_3,sps)

# bio_4 <- bio[[c(14)]]
# bio_4 <- mask(bio_4,sps)
# bio_4 <-crop(bio_4,sps)

# bio_5 <- bio[[c(15)]]
# bio_5 <- mask(bio_5,sps)
# bio_5 <-crop(bio_5,sps)

# bio_6 <- bio[[c(16)]]
# bio_6 <- mask(bio_6,sps)
# bio_6 <-crop(bio_6,sps)

# bio_7 <- bio[[c(17)]]
# bio_7 <- mask(bio_7,sps)
# bio_7 <-crop(bio_7,sps)

# bio_8 <- bio[[c(18)]]
# bio_8 <- mask(bio_8,sps)
# bio_8 <-crop(bio_8,sps)

# bio_9 <- bio[[c(19)]]
# bio_9 <- mask(bio_9,sps)
# bio_9 <-crop(bio_9,sps)

# bioclim <- stack(c(bio_1, bio_10, bio_11,
#                    bio_12, bio_13, bio_14,
#                    bio_15, bio_16, bio_17,
#                    bio_18, bio_19, bio_2,
#                    bio_3, bio_4, bio_5,
#                    bio_6, bio_7, bio_8,
#                    bio_9))

# create directory to hold clipped and masked bioclim variables
# only create this directory once 
# dir.create("vosnesenskii")

# setwd()
setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW/vosnesenskii")

# write the rasters
# rasters are already made so only do this once
writeRaster(bioclim, filename=names(bioclim), bylayer=TRUE, format="raster", overwrite = TRUE)

# load the rasters back into the session
bio_1 <-raster("bio_1.grd")
bio_2 <-raster("bio_2.grd")
bio_3 <-raster("bio_3.grd")
bio_4 <-raster("bio_4.grd")
bio_5 <-raster("bio_5.grd")
bio_6 <-raster("bio_6.grd")
bio_7 <-raster("bio_7.grd")
bio_8 <-raster("bio_8.grd")
bio_9 <-raster("bio_9.grd")
bio_10 <-raster("bio_10.grd")
bio_11 <-raster("bio_11.grd")
bio_12 <-raster("bio_12.grd")
bio_13 <-raster("bio_13.grd")
bio_14 <-raster("bio_14.grd")
bio_15 <-raster("bio_15.grd")
bio_16 <-raster("bio_16.grd")
bio_17 <-raster("bio_17.grd")
bio_18 <-raster("bio_18.grd")
bio_19 <-raster("bio_19.grd")

# turn them into a raster brick
bioclim <- stack(c(bio_1, bio_10, bio_11,
                      bio_12, bio_13, bio_14,
                      bio_15, bio_16, bio_17,
                      bio_18, bio_19, bio_2,
                      bio_3, bio_4, bio_5,
                      bio_6, bio_7, bio_8,
                      bio_9))


##obtain identifiers of the cells where species occures (lat/long).
nrow(df)
points_df <- data.frame(df)[1:nrow(df), c("decimalLongitude","decimalLatitude")]

##extract the value from each cell and put it in one centeral location
df_cell_id <- cellFromXY(subset(bioclim,1),points_df)

##principal component analysis, convert raster object into a data frame to run PCA and
##remove non-defined area from the dataset, gives values to every point 
bioclim_df <- na.omit(as.data.frame(bioclim))

##check and see what it did
head(bioclim_df)
nrow(bioclim_df)

##pca study over the whole study area
##dudi.pca performs a principal component analysis of a data frame and 
##returns the results as objects of class pca and dudi.

##PCA scores on the first two axes
pca.df <- dudi.pca(bioclim_df, scannf = F, nf = 2)

##plot(pca_ZA$li[, 1:2])
par(mfrow=c(1, 1))
plot(pca.df$li[, 1:2])

##tail of distributions
## allows me to get rid of abherrant data
sort(pca.df$li[, 1], decreasing = TRUE)[1:7]


##IDs of points to remove
to_remove <- which(pca.df$li[, 1] > 11)
to_remove

## remove points and recompute PCA
nrow(bioclim_df)
if(length(to_remove)){ ## remove outliers
  bioclim_df <- bioclim_df[ - to_remove,]
  pca.df <- dudi.pca(bioclim_df,scannf = F, nf = 2)  
}
nrow(bioclim_df)

# plot the PCA again
pca.df <- dudi.pca(bioclim_df,scannf = F, nf = 2)
plot(pca.df$li[, 1:2])


# save picture
tiff("Fig1_bioclim_variable_selection.tiff",
     width = 6, height = 4, units = 'in', res = 300)
par(mfrow=c(1, 2))

##Discrminiate species prescences from the entire hawaiian enviromental space. (scatterplot)
### need to figure out what df = 5 means
s.class(pca.df$li[, 1:2],fac = factor(rownames(bioclim_df)%in% df_cell_id,
                                      levels = c("FALSE","TRUE"),
                                      labels = c("background","df")),
        col = c("light gray","black"),csta = 0,
        cellipse = 2,
        cpoint = .3,
        pch = 16)

##create 2 sides of the plot a and b
mtext("(a)",side = 3,line = 3,adj = 0)

##second side of the plot panel b to help decide which variables to keep
s.corcircle(pca.df$co,clabel = .5)
mtext("(b)",side = 3, line = 3, adj = 0)

dev.off()
##subselect four variables from the full enviromental set
bioclim_sub <- stack(subset(bioclim, c("bio_2","bio_9","bio_15","bio_18")))

##Biomod2 modeling procedure

##next we put the data into the right format by using the BIOMOD_formattingData function
df_data <- BIOMOD_FormatingData(resp.var = rep(1,nrow(df)),
                                   expl.var = bioclim_sub,
                                   resp.xy = df[,c('decimalLongitude','decimalLatitude')],
                                   resp.name = "B. vosnesenskii",
                                   PA.nb.rep = 3,PA.nb.absences = 10000,
                                   PA.strategy = 'random')

##check to see if it worked
df_data

#plot of selected pseudo-absences
##displays the location of the psuedo absences in the three datasets compared to the species occurences
plot(df_data)

##this line picks the modeling options
df_opt <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic',
                                               interaction.level = 1),
                                    GBM = list(n.trees = 1000),
                                    GAM = list(algo ='GAM_mgcv'))


##next line we run the model with 4 diffrent models
df_models <- BIOMOD_Modeling(data = df_data,
                             models = c("GLM","GBM","RF","GAM",
                                        "CTA","ANN","SRE","FDA",
                                        "MARS"),
                             models.options = df_opt, 
                             NbRunEval = 4, DataSplit = 80,
                             VarImport = 3, do.full.models = F,
                             modeling.id = "ex.2")


##get model evaluation scores 
xsono_models_scores <- get_evaluations(df_models)

#xsono_models_scores is a 5 deminsional array containing the scores for the models

##dim: Get or set the number of rows, columns, and layers of a Raster* object.
##You cannot use this function to set the dimensions of a RasterStack object.
dim(xsono_models_scores)

##Retrieve or set the dimnames of an object.
dimnames(xsono_models_scores)

##model_scores_graph: This function is a graphic tool to represent evaluation scores of models
##produced with biomod2 according to 2 different evaluation methods. 

##Models can be grouped in several ways (by algo, by CV run, ...) 
##to highlight potential differences in models quality due to chosen models,
##cross validation sampling bias. 
##Each point represents the average evaluation score across each group. 
##Lines represents standard deviation of evaluation scores of the group.

##GAM, GBM, GLM, RF scores
models_scores_graph(df_models, by = "models",metrics = c("ROC","TSS"),
                    xlim = c(0.5,1), ylim = c(0.5,1))

##points represent mean of evaluation score for a given condition 
##and the lines represent the associated standard deviations

## get variable importance in the models
(xsono_models_var_import <- get_variables_importance(df_models))

##RUN1, RUN2, RUN3, RUN4 scores
models_scores_graph(df_models, by = "cv_run",metrics = c("ROC","TSS"),
                    xlim = c(0.5,1), ylim = c(0.5,1))

##PA1, PA2, PA3 scores
models_scores_graph(df_models, by = "data_set",metrics = c("ROC","TSS"),
                    xlim = c(0.5,1), ylim = c(0.5,1))

#calculate the MEAN of variable importance by algorighm
apply(xsono_models_var_import, c(1,2), mean)

xsono_models <- df_models

#to do this we first have to load the produced models
xsono_glm <- BIOMOD_LoadModels(xsono_models, models = 'GLM')
xsono_gbm <- BIOMOD_LoadModels(xsono_models, models = 'GBM')
xsono_rf <- BIOMOD_LoadModels(xsono_models, models = 'RF')
xsono_cta <- BIOMOD_LoadModels(xsono_models, models = 'CTA')
xsono_ann <- BIOMOD_LoadModels(xsono_models, models = 'ANN')
xsono_sre <- BIOMOD_LoadModels(xsono_models, models = 'SRE')
xsono_fda <- BIOMOD_LoadModels(xsono_models, models = 'FDA')
xsono_mars <- BIOMOD_LoadModels(xsono_models, models = 'MARS')
xsono_gam <- BIOMOD_LoadModels(xsono_models, models = 'GAM')
# xsono_phi <- BIOMOD_LoadModels(xsono_models, models = 'MAXENT.Phillips')
# xsono_tsu <- BIOMOD_LoadModels(xsono_models, models = 'MAXENT.Tsuruoka')

##these are graphical visualizations of the response curve of each varaible
##how does each enviromental variable influence probability of prescence?
##each line cooresponds to a diffrent model
glm_eval_strip <- biomod2::response.plot2(models = xsono_glm, Data = get_formal_data(xsono_models,'expl.var'),
                                          show.variables = get_formal_data(xsono_models,'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
gbm_eval_strip <- biomod2::response.plot2(models = xsono_gbm, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
rf_eval_strip <- biomod2::response.plot2(models = xsono_rf, Data = get_formal_data(xsono_models, 'expl.var'),
                                         show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                         do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                         display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
gam_eval_strip <- biomod2::response.plot2(models = xsono_gam, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
cta_eval_strip <- biomod2::response.plot2(models = xsono_cta, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
ann_eval_strip <- biomod2::response.plot2(models = xsono_ann, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
sre_eval_strip <- biomod2::response.plot2(models = xsono_sre, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
fda_eval_strip <- biomod2::response.plot2(models = xsono_fda, Data = get_formal_data(xsono_models, 'expl.var'),
                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
mars_eval_strip <- biomod2::response.plot2(models = xsono_mars, Data = get_formal_data(xsono_models, 'expl.var'),
                                           show.variables = get_formal_data(xsono_models, 'expl.var.names'),
                                           do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
                                           display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
# phi_eval_strip <- biomod2::response.plot2(models = xsono_phi, Data = get_formal_data(xsono_models, 'expl.var'),
#                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
#                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
#                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))
# tsu_eval_strip <- biomod2::response.plot2(models = xsono_tsu, Data = get_formal_data(xsono_models, 'expl.var'),
#                                          show.variables = get_formal_data(xsono_models, 'expl.var.names'),
#                                          do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
#                                          display_title = FALSE, data_species = get_formal_data(xsono_models, 'resp.var'))

##Ensembal modeling
##to reduce the number of outputs we only concider two "ensembaling" options:
##committee averaging and weighted mean. we also produce coefficent of variation 
##that tell us the extent the models agree or disagree.

xsono_ensemble_models <- BIOMOD_EnsembleModeling(modeling.output = xsono_models,
                                                 em.by = 'all',
                                                 eval.metric = 'TSS',
                                                 eval.metric.quality.threshold = 0.6,
                                                 models.eval.meth = c('KAPPA','TSS','ROC'),
                                                 prob.mean = FALSE,
                                                 prob.cv = TRUE, 
                                                 committee.averaging = TRUE,
                                                 prob.mean.weight = TRUE,
                                                 VarImport = 0)

##now check the scores for the ensembal models
(xsono_ensemble_models_scores <- get_evaluations(xsono_ensemble_models))

##Current projections##
xsono_models_proj_current <- BIOMOD_Projection(modeling.output = xsono_models,
                                               new.env = bioclim_sub,
                                               proj.name = "current", 
                                               binary.meth = 'TSS',
                                               output.format = ".img", 
                                               do.stack = FALSE)

xsono_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(EM.output = xsono_ensemble_models,
                                                                 projection.output = xsono_models_proj_current,
                                                                 binary.meth = "TSS",
                                                                 output.format = "img",
                                                                 do.stack = FALSE)

#Load 2050 bioclim variables
bioclim_world_2050_BC45 <- stack(c(bio_2 = "WorldClim_data/2050/BC_45/bc45bi502.tif",
                                   bio_9 = "WorldClim_data/2050/BC_45/bc45bi509.tif",
                                   bio_15 = "WorldClim_data/2050/BC_45/bc45bi5015.tif",
                                   bio_18 = "WorldClim_data/2050/BC_45/bc45bi5018.tif"))

#crop of our area
bioclim_HI_2050_BC45 <- crop(bioclim_world_2050_BC45, hawaii)
bioclim_HI_2050_BC45 <- mask(bioclim_HI_2050_BC45, hawaii)
bioclim_HI_2050_BC45 <- stack(bioclim_HI_2050_BC45)

##Save this raster stack on the hard drive if needed
xsono_models_proj_2050_BC45 <- BIOMOD_Projection(modeling.output = xsono_models,
                                                 new.env = bioclim_HI_2050_BC45,
                                                 proj.name = "2050_BC45", 
                                                 binary.meth = "TSS",
                                                 output.format = '.img',
                                                 do.stack = FALSE)

xsono_ensemble_models_proj_2050_BC45 <- BIOMOD_EnsembleForecasting(EM.output = xsono_ensemble_models,
                                                                   projection.output = xsono_models_proj_2050_BC45, 
                                                                   binary.meth = "TSS",
                                                                   output.format = ".img",
                                                                   do.stack = FALSE)

#Load 2070 bioclim variables
bioclim_world_2070_BC45 <- stack(c(bio_2 = "WorldClim_data/2070/BC_45/bc45bi702.tif",
                                   bio_9 = "WorldClim_data/2070/BC_45/bc45bi709.tif",
                                   bio_15 = "WorldClim_data/2070/BC_45/bc45bi7015.tif",
                                   bio_18 = "WorldClim_data/2070/BC_45/bc45bi7018.tif"))

#crop of our area
bioclim_HI_2070_BC45 <- crop(bioclim_world_2070_BC45, hawaii)
bioclim_HI_2070_BC45 <- mask(bioclim_HI_2070_BC45, hawaii)
bioclim_HI_2070_BC45 <- stack(bioclim_HI_2070_BC45)

#you may save these raster on the hard drive
xsono_models_proj_2070_BC45 <- BIOMOD_Projection(modeling.output = xsono_models,
                                                 new.env = bioclim_HI_2070_BC45,
                                                 proj.name = "2070_BC45",
                                                 binary.meth = "TSS",
                                                 output.format = ".img",
                                                 do.stack = FALSE)

xsono_ensemble_models_proj_2070_BC45 <- BIOMOD_EnsembleForecasting(EM.output = xsono_ensemble_models,
                                                                   projection.output = xsono_models_proj_2070_BC45,
                                                                   binary.meth = "TSS",
                                                                   output.format = ".img",
                                                                   do.stack = FALSE)

#get the ensembal models projection stack
stk_xsono_ef_2070_BC45 <- get_predictions(xsono_ensemble_models_proj_2070_BC45)

#keep committee averaging and weighted mean projections only ensamble models 
stk_xsono_ef_2070_BC45 <- subset(stk_xsono_ef_2070_BC45, grep("EMca|EMwmean", names(stk_xsono_ef_2070_BC45)))

#simplify the layer names for plotting conveniences
names(stk_xsono_ef_2070_BC45) <-sapply(strsplit(names(stk_xsono_ef_2070_BC45), "_"), getElement, 2)

## make sure library(rasterVis) and library(lattice) are on
levelplot(stk_xsono_ef_2070_BC45, main = "X. Sonorina ensemble projections\nin 2070 with BC45",
          col.regions = colorRampPalette(c("grey90","yellow4","green4"))(100))

##species range change
#load binary projections
xsono_bin_proj_current <- stack(c(ca = "x.sonorina/proj_current/individual_projections/x.sonorina_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbinimg.grd",
                                  wm = "x.sonorina/proj_current/individual_projections/x.sonorina_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbinimg.grd"))
xsono_bin_proj_2050_BC45 <- stack(c(ca = "x.sonorina/proj_2050_BC45/individual_projections/x.sonorina_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                                    wm = "x.sonorina/proj_2050_BC45/individual_projections/x.sonorina_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"))
xsono_bin_proj_2070_BC45 <- stack(c(ca = "x.sonorina/proj_2070_BC45/individual_projections/x.sonorina_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                                    wm = "x.sonorina/proj_2070_BC45/individual_projections/x.sonorina_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"))

###SRC means species range change
##calculate the range change from current to 2050 and 2070
#SRC current -> 2050
SRC_current_2050_BC45 <- BIOMOD_RangeSize(xsono_bin_proj_current, xsono_bin_proj_2050_BC45)
SRC_current_2050_BC45$Compt.By.Models

#SRC current -> 2070
SRC_current_2070_BC45 <- BIOMOD_RangeSize(xsono_bin_proj_current, xsono_bin_proj_2070_BC45)
SRC_current_2070_BC45$Compt.By.Models

### plot the predicted changes
xsono_src_map <- stack(SRC_current_2050_BC45$Diff.By.Pixel, SRC_current_2070_BC45$Diff.By.Pixel)
names(xsono_src_map) <- c("ca cur-2050", "wm cur-2050", "ca cur-2070", "wm cur-2070")

my.at <- seq(-2.5,1.5,1)
myColorkey <- list(at = my.at, labels = list(c("lost", "pres", "abs", "gain"),at = my.at[-1]-0.5))
rasterVis::levelplot(xsono_src_map, main = "X. sonorina range change", colorkey = myColorkey, layout = c(2,2))
## should see some plots

##impact of model scenario/time slice
##choose a refrence for which all other model predictions will be compared
ref <- subset(xsono_bin_proj_current, "ca")

#define the facets(variables) we want to study
mods <- c("GLM", "GBM", "RF", "GAM","CTA","ANN","SRE","FDA","MARS","MAXENT.Phillips","MAXENT.Tsuruoka","caByTSS","wmeanByTSS")
data_set <- c("PA1", "PA2", "PA3", "mergedData")
cv_run <- c("RUN1", "RUN2", "RUN3", "RUN4", "mergedRun")

#construct combination of all facets
groups <- as.matrix(expand.grid(models = mods, data_set = data_set, cv_run = cv_run, stringsAsFactors = FALSE))

#load all projections we have produced
all_bins_proj_files <- list.files(path = "x.sonorina", pattern = "_TSSbin.img$", full.names = TRUE, recursive = TRUE)

#we want to focus on current vs 2070. We thus removed the projections of 2050
current_and_2070_proj_files <- grep(all_bins_proj_files, pattern = "2070",value=T)

#only keep projections that match our selected facets groups pg. 384
selected_bin_proj_files <- apply(groups,1,
                                 function(x){
                                   proj_file<-NA 
                                   match_tab <- sapply(x, grepl, current_and_2070_proj_files) 
                                   match_id <- which(apply(match_tab,1,all))
                                   if(length(match_id)) proj_file <- current_and_2070_proj_files[match_id]
                                   return(proj_file)
                                 })

#remove non-matching groups
to_remove <- which(is.na(selected_bin_proj_files))
if(length(to_remove)){
  groups <- groups[-to_remove, ]
  selected_bin_proj_files <- selected_bin_proj_files[-to_remove]}

#build a stack of selected projections
proj_groups <- stack(selected_bin_proj_files)

ProbDensFunc(initial = ref, projections = proj_groups, groups = t(groups),plothist = FALSE,csvn = FALSE)