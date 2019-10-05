# Jonathan B. Koch 
# Script: 01_PacNW_appositus_occurrences_crop_bioclim_vars.R
# Paper: Climate Change Pacific Northwest Bumble Bees
# Date: 04 October 2019
# Website for help with biomod2 scripts: # www.unil.ch/hsdm

## set where the program saves and looks for data ##
setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW/occurences_filtered")

# list files
list.files()

# load libraries
library(biomod2)
library(ggplot2)
library(ade4)
library(raster)
library(maptools)
library(maps)
library(mapdata)
library(sp) # can't be installed with HPC version of R


## load species data into R
df <- read.table("vosnesenskii_unique_occurrences_filterd_gbif.csv",
                 header=TRUE, 
                 sep=",",
                 row.names=NULL)

# read in shapefile
world <- map_data("world") # from ggplot2
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
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "yellow", size = 2)

# B. vosnesnskii, remove occurence that don't match Koch et al. 2012 or Williams et al. 2014
nrow(df)
df <- subset(df, decimalLongitude <= -115 )

# check out data
ggplot() + 
  geom_polygon(data = north_america, aes(x=long, y = lat, group = group)) + 
  coord_sf(xlim = c(min_x, max_x), ylim = c(min_y, max_y), expand = FALSE) + #xlim = longitude, ylim = latitude
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "yellow", size = 2)

# create boundary of data visualization
nrow(df)
max_y = max(df$decimalLatitude+0.5)
min_y = min(df$decimalLatitude-0.5)
max_x = max(df$decimalLongitude+0.5)
min_x = min(df$decimalLongitude-0.5)

# check out data
ggplot() + 
  geom_polygon(data = north_america, aes(x=long, y = lat, group = group)) + 
  coord_sf(xlim = c(min_x, max_x), ylim = c(min_y, max_y), expand = FALSE) + #xlim = longitude, ylim = latitude
  geom_point(data = df, aes(x = df$decimalLongitude, y = decimalLatitude), color = "yellow", size = 2)


# create mask from max/min latitude (y) and longitude (x)
# create spatial polygons
# https://rstudio-pubs-static.s3.amazonaws.com/202536_7a122ff56e9f4062b6b012d9921afd80.html

# make bounding box
# to clip bioclimatic rasters
# no need to this is if you did it the first time
x_coord <- c(max_x, max_x, min_x, min_x)
y_coord <- c(min_y, max_y, max_y, min_y)
xym <- cbind(x_coord, y_coord)
xym

# polygon library library(sp)
# don't need to make polygon after you
# you clipped the B. vosnesenskii extent
# for bioclimatic variables
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

# Enviromental variables selection
# stack of all your climate variables together
# change directory to location of bio

setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW")
bio <- stack(list.files("bio",pattern = "bio",
                        full.names = T), RAT = FALSE)

bio_1 <- bio[[c(1)]]
bio_1 <- mask(bio_1,sps)
bio_1 <-crop(bio_1,sps)

bio_10 <- bio[[c(2)]]
bio_10 <- mask(bio_10,sps)
bio_10 <-crop(bio_10,sps)

bio_11 <- bio[[c(3)]]
bio_11 <- mask(bio_11,sps)
bio_11 <-crop(bio_11,sps)

bio_12 <- bio[[c(4)]]
bio_12 <- mask(bio_12,sps)
bio_12 <-crop(bio_12,sps)

bio_13 <- bio[[c(5)]]
bio_13 <- mask(bio_13,sps)
bio_13 <-crop(bio_13,sps)

bio_14 <- bio[[c(6)]]
bio_14 <- mask(bio_14,sps)
bio_14 <-crop(bio_14,sps)

bio_15 <- bio[[c(7)]]
bio_15 <- mask(bio_15,sps)
bio_15 <-crop(bio_15,sps)

bio_16 <- bio[[c(8)]]
bio_16 <- mask(bio_16,sps)
bio_16 <-crop(bio_16,sps)

bio_17 <- bio[[c(9)]]
bio_17 <- mask(bio_17,sps)
bio_17 <-crop(bio_17,sps)

bio_18 <- bio[[c(10)]]
bio_18 <- mask(bio_18,sps)
bio_18 <-crop(bio_18,sps)

bio_19 <- bio[[c(11)]]
bio_19 <- mask(bio_19,sps)
bio_19 <-crop(bio_19,sps)

bio_2 <- bio[[c(12)]]
bio_2 <- mask(bio_2,sps)
bio_2 <-crop(bio_2,sps)

bio_3 <- bio[[c(13)]]
bio_3 <- mask(bio_3,sps)
bio_3 <-crop(bio_3,sps)

bio_4 <- bio[[c(14)]]
bio_4 <- mask(bio_4,sps)
bio_4 <-crop(bio_4,sps)

bio_5 <- bio[[c(15)]]
bio_5 <- mask(bio_5,sps)
bio_5 <-crop(bio_5,sps)

bio_6 <- bio[[c(16)]]
bio_6 <- mask(bio_6,sps)
bio_6 <-crop(bio_6,sps)

bio_7 <- bio[[c(17)]]
bio_7 <- mask(bio_7,sps)
bio_7 <-crop(bio_7,sps)

bio_8 <- bio[[c(18)]]
bio_8 <- mask(bio_8,sps)
bio_8 <-crop(bio_8,sps)

bio_9 <- bio[[c(19)]]
bio_9 <- mask(bio_9,sps)
bio_9 <-crop(bio_9,sps)

bioclim <- stack(c(bio_1, bio_10, bio_11,
                   bio_12, bio_13, bio_14,
                   bio_15, bio_16, bio_17,
                   bio_18, bio_19, bio_2,
                   bio_3, bio_4, bio_5,
                   bio_6, bio_7, bio_8,
                   bio_9))


# create directory to hold clipped and masked bioclim variables
# only create this directory once 
# dir.create("vosnesenskii")

# setwd()
setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW/bioclim_cropped/vosnesenskii")

# WRITE the RASTERS
# rasters are already made so only do this once
writeRaster(bioclim, filename=names(bioclim), bylayer=TRUE, format="raster", overwrite = TRUE)

# load the rasters back into the session
list.files()
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

# setwd
setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW/occurences_pca_filtered")

# write df to .csv so it can be save for the final analysis
write.csv(df, "vosnesenskii_unique_occurrences_filterd_PCA_gbif.csv",
          row.names = FALSE)