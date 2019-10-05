#!/bin/bash
#SBATCH -J 01_R_BIOMOD.job # Name for your job
#SBATCH -c 10 # Number of cores requested, a node has 20 cores in the community.q partition
#SBATCH -n 1 # Number of tasks when using mpi.
#SBATCH --mem-per-cpu 6400 #the minimum amount of memory per core to request
#SBATCH -N 1 # Number of nodes, Current Maximum is 20 nodes in the community.q partition -
#SBATCH -t 4200 # Runtime in minutes. The Maximum runtime currently is 72 hours, 4320 minutes -
#SBATCH -p shared #community.q is Partition to submit to t
#SBATCH -o job.01_R_BIOMOD-%A.out # Standard out goes to this file
#SBATCH -e job.01_R_BIOMOD-%A.err # Standard err goes to this file

# useful website for r scripts on HPC
# https://www.hawaii.edu/its/ci/hpc-tutor/using-r/
# script to do R in cluster
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