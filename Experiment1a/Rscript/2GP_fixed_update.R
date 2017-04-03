###################################################################
##                    Gaussian Update                            ##
## This script performs a Gaussian Update of the GIA processes   ##
## Given                                                         ##
## 1. fixed parameter values                                     ##
## 2. subset of GPS obs                                          ##
## 3. GIA forward model prior                                    ##
###################################################################
library(rgl)
library(GEOmap)
library(INLA)
library(MVST)
load("o:/glbm/Experiment1a/Mesh_GIA.RData")