#### This file load the GIA prior and GPS data and prepare the data for the BHM model

#### 1 Load GIA prior
if(Sys.info()["sysname"] == "Windows"){
  ice6g <- read.table("Z:/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  ice6g <- read.table("~/globalmass/WP2-SolidEarth/GIAforwardModels/textfiles/GIA_Pel-6-VM5.txt", header = T)
}else{
  ice6g <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}

#### 1 Load GPS data

if(Sys.info()["sysname"] == "Windows"){
  GPSV4b <- read.table("Z:/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  GPSV4b <- read.table("~/globalmass/WP2-SolidEarth/GPS/NGL/BHMinputFiles/GPS_v04b.txt", header = T)
}else{
  GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}

