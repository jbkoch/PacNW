# Run R script for slurm
# Jonathan B. Koch 
# Orginator: Jesse Tabor (UH-Hilo)
# BiodMod2: X. sonorina
# 4 October 2019

# load libraries
library(biomod2)
library(raster)
library(ade4)

# pathway on kochj on HPC @UH
# "/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/.."

# Get masked bioclim grids

# get data
# setwd()
# setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW") #Jon Mac
setwd("/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/occurences")
df <- read.csv("vosnesenskii_unique_occurrences_filterd_PCA_gbif.csv", 
               sep = ",",
               header = TRUE)

# raster data
# setwd("/Users/jonathankoch/Google Drive/git_myrepo/PacNW/vosnesenskii")
setwd("/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/vosnesenskii")
list.files()

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
bioclim

## subselect variables from the full enviromental set based on
## PCA results of target species
bioclim_sub <- stack(subset(bioclim, c("bio_4","bio_7","bio_2",
                                       "bio_14", "bio_18","bio_5",
                                       "bio_10", "bio_1", "bio_3",
                                       "bio_11", "bio_6", "bio_15",
                                       "bio_19", "bio_12", "bio_9")))

bioclim_sub

##Biomod2 modeling procedure

##next we put the data into the right format by using the BIOMOD_formattingData function
##default pseudo-absence is 10,000
df_data <- BIOMOD_FormatingData(resp.var = rep(1,nrow(df)),
                                expl.var = bioclim_sub,
                                resp.xy = df[,c('decimalLongitude','decimalLatitude')],
                                resp.name = "Bombus vosnesenskii",
                                PA.nb.rep = 3,PA.nb.absences = 10000,
                                PA.strategy = 'random')

df_data
#plot of selected pseudo-absences
##displays the location of the psuedo absences in the three datasets compared to the species occurences

# commit files outside of git 
# setwd("/Users/jonathankoch/Documents/PacNW")
setwd("/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/Figures/vosnesenskii")
pdf("Fig2_inputdata_PA.pdf")
par(mfrow=c(1, 1))
plot(df_data)
dev.off()

# change directory for the files of modeling
setwd("/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/SDMresults/vosnesenskii/output")

##this line picks the modeling options
df_opt <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = "/mnt/lts/nfs_fs02/hilobio_all/kochj/biomod_inR/PacNW-master/maxent_files"))

##next line we run the model with 4 diffrent models
df_models <- BIOMOD_Modeling(data = df_data,
                             models = c("MAXENT.Phillips", "GLM"),
                             models.options = df_opt, 
                             NbRunEval = 4, DataSplit = 80,
                             VarImport = 3, do.full.models = F,
                             modeling.id = "ex.2")

##get model evaluation scores 
xsono_models_scores <- get_evaluations(df_models)
xsono_models_scores

# setwd
write.txt(xsono_models_scores, "xsono_models_scores.txt")
