#!/usr/bin/env Rscript


##The intent of this script is to 'tune' crown delineation parameters
## all parameters will be varied randomly withing set range
## and run 10000 times
## Author:Jack
## Start date: 2/21/19

#All methods are based on a seed/marker generated from  lidR::tree_detection()

## All methods are from the lidR package: 
# https://cran.r-project.org/web/packages/lidR/lidR.pdf
# https://github.com/Jean-Romain/lidR/wiki

## set seed so that the random parameter tuning is reproducible 
# set.seed(1)



##Load in required packages
library(lidR)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(sp)
library(maptools)
library(dplyr)
library(rgeos)
library(snow)
library(cleangeo)
max_iteration = 300 ## how many loops?

#load in CHM, LAS, and plots
chm <- raster("CHM_trim.tif")
##chm <- raster::focal(chm, w=matrix(1,3,3), fun=mean) 
las <- readLAS("las_trim.las")
plots <- readOGR("20m_plots.shp"); names(plots) <- c("plotID")


#_____MANUAL CROWN PREP
manual <- readOGR("manual_crown_delineation_fixed_geom.shp")
man_geom_report <- clgeo_CollectionReport(manual)
if(length(unique(man_geom_report$valid == 2))){
  manual <- clgeo_Clean(manual)
}

# clean it up a bit
manual@data$crownID <-NULL

manual$man.ID <- 1:length(manual) # give unique ID
manual$man.area <- area(manual) # calculate area
manual <- manual[,c(2,1,3)] # reorder my columns
# -- calculate manual centroid, and remove any polygons where the centroid falls outside plot bounds
man_cent <- gCentroid(manual, byid = T)
manIN <-over(man_cent, plots) # calculate if centroids fall WITHIN plot boundaries; this will also assign plot IDs
manIN$man.ID <-1:length(manIN$plotID) # make identifer column, with same name as manual to join on 
manual@data <- dplyr::left_join(manual@data, manIN) ## joined 
## now we need to remove the polygons outside the plot
manual <- manual[!is.na(manual@data$plotID),]
# remove manual excess -- left with just manual
rm(man_cent, manIN, man_geom_report); gc()
#
#

#
#
accuracy_assessment <- function(auto, manual, plots){
  ## Auto crowns will be spit out with just one column: Tree ID
  auto$aut.ID <- 1:length(auto@data$treeID) # assign polygon ID
  auto <- auto[,-1]# remove treeID -- useless?
  
  if(identicalCRS(auto, manual) == FALSE){ # I think CRS will always be wrong -- this checks and corrects
    auto <- spTransform(auto, crs(manual))
  }
  plot_buf <- gBuffer(plots, width = 3, byid = T)
  aut_cent <- gCentroid(auto, byid = T)# calculate centroid
  autIN <- over(aut_cent, plot_buf) #find centroids in polygons -- this also assigns plot ID
  autIN$aut.ID <- 1:length(autIN$plotID) # assign polygon ID to be joined on.
  auto@data <- dplyr::left_join(auto@data, autIN) # join centroid 'over' data with auto 
  auto <- auto[!is.na(auto@data$plotID),] # remove crowns that are not in plot
  auto$aut.area <- area(auto) # calculate area of autos
  row.names(auto@data)<-NULL
  auto<- gBuffer(auto,byid = T,  width = 0 )
  
  aut_geom_report <- clgeo_CollectionReport(auto)
  if (length(unique(aut_geom_report$valid == 2))){
    print("Fixing broken geometry")
    auto <- clgeo_Clean(auto)
  }
  
  
  #intersection between auto and manual
  inter <- raster::intersect(auto, manual)
  inter <- inter[,c(4,1,2,6,3,5)] # rearrange, and drop duplicate colums
  colnames(inter@data)[colnames(inter@data)=="plotID.1"] <- "plot.ID" # clean up colnames
  # calculate intersection area
  inter$int.area <- area(inter)
  for (i in 1:length(inter@polygons)){
    inter@polygons[[i]]@ID <- as.character(i)
    
  }
  inter <- gBuffer(inter, byid=TRUE, width=0)
  
  
  
  inter$int.man <- inter$int.area/inter$man.area # ratio of intersected area to manual area
  inter$int.aut <- inter$int.area/inter$aut.area # ratio of intersected area to auto area
  #determine number of true positives
  #currently this will just classify an intersection as a true positive (1) or a miss (0)
  inter$TP <- ifelse(inter$int.man >= 0.5 & inter$int.aut >= 0.5, 1, 0)
  inter@data$TP <- as.factor(inter@data$TP)
  
  #compile data frame with slots for TP count
  
  tpX <- as.data.frame(inter@data %>%
                         group_by(plot.ID, TP)%>%
                         tally())
  tpX$TP <- as.numeric(as.character(tpX$TP))
  #this unfortunate chunk of code is for the circumstances when a plot as no TP
  for(i in 1:NROW(tpX)){
    if(length(tpX[tpX$TP[[i]] == 1])== 0){
      tpX[i, c(2:3)] <- c(1,0)}
  }
  TP_agg <- as.data.frame(tpX %>%
                            group_by(plot.ID, TP, n) %>%
                            tally())
  TP_agg <- TP_agg[order(TP_agg$plot.ID, TP_agg$n, decreasing = T),]
  TP_agg <- TP_agg[!duplicated(TP_agg$plot.ID),c(1:3)]
  TP_agg <- TP_agg[,c(1,3)]
  names(TP_agg) <- c("plotID", "TP")
  
  #compile data.frame with counts of manual crowns per plot
  man_agg <- as.data.frame(count(manual@data, plotID))
  names(man_agg) <- c("plotID", "MC")
  
  #compile data.frames with counts of auto crowns per plot
  aut_agg <- as.data.frame(count(auto@data, plotID))
  names(aut_agg) <- c("plotID", "AC")
  
  
  # merge to one data.frame
  CrownCount <- merge(TP_agg, man_agg, by = "plotID")
  CrownCount <- merge(CrownCount, aut_agg, by = "plotID")
  #add new column for "true positive accuracy
  # simply the ratio of true positives to manually delineated crowns 
  CrownCount$accuracy <- CrownCount$TP/CrownCount$MC
  CrownCount$plotID <- as.integer(as.character(CrownCount$plotID))
  CrownCount <- CrownCount[order(CrownCount$plotID), ]
  #Overall_Accuracy <- (sum(CrownCount$TP))/(sum(CrownCount$CC)) # this should likely just be left as a post function analysis
  
  
  return(CrownCount) ## Can easily change this to return a table to look at individual plots
}
## delineation method: simple watershed
# Cite:


## build empty lists before loop -- these will be populated with iteration, parameters, and accuracy score
iterate <- c() #necessary? counts rep
TH_Tree <- c()
EXT <- c()
TOL <- c()
#-error
nMC <-c() #number of manual crowns
nAC <-c() # number of automatic crowns
Acc <-c() #overall accuracy
#plot 1-15 accuracies
P1 <-c() 
P2 <-c()
P3 <-c()
P4 <-c()
P5 <-c()
P6 <-c()
P7 <-c()
P8 <-c()
P9 <-c()
P10 <-c()
P11 <-c()
P12 <-c()
P13 <-c()
P14 <-c()
P15 <-c()

for (i in 1:max_iteration){
  print(i)
  ##define parameter variablity
  thTree <-signif(runif(1,2,10),4) 
  ext <-sample(1:3,1)
  tol <-signif(runif(1,0.001,1),4)
  #________________________________________________________________
  
  ws_crowns <- lastrees(las, watershed(chm = chm, th_tree = thTree, tol = tol, ext = ext))
  auto <- tree_hulls(ws_crowns, type = "concave")
  
  #_________________________________________________________________
   ##Error
   CrownCount<- accuracy_assessment(auto, manual,plots); print (CrownCount)
   Overall_Accuracy <- (sum(CrownCount$TP))/(sum(CrownCount$MC))
   Nmc <- sum(CrownCount$MC)
   Nac <- sum(CrownCount$AC)
   rownames(CrownCount)<- paste0("P",CrownCount$plotID)
   CCtrans <- as.data.frame(t(as.matrix(CrownCount)))
   
   Acc[[i]] <- Overall_Accuracy
   nMC[[i]] <- Nmc
   nAC[[i]] <- Nac
   if((length(CCtrans$P1) == 0)){P1[[i]]<-0}else {P1[[i]]<-CCtrans$P1[[5]]}  
   if((length(CCtrans$P2) == 0)){P2[[i]]<-0}else {P2[[i]]<-CCtrans$P2[[5]]}  
   if((length(CCtrans$P3) == 0)){P3[[i]]<-0}else {P3[[i]]<-CCtrans$P3[[5]]}  
   if((length(CCtrans$P4) == 0)){P4[[i]]<-0}else {P4[[i]]<-CCtrans$P4[[5]]}  
   if((length(CCtrans$P5) == 0)){P5[[i]]<-0}else {P5[[i]]<-CCtrans$P5[[5]]}  
   if((length(CCtrans$P6) == 0)){P6[[i]]<-0}else {P6[[i]]<-CCtrans$P6[[5]]}  
   if((length(CCtrans$P7) == 0)){P7[[i]]<-0}else {P7[[i]]<-CCtrans$P7[[5]]}  
   if((length(CCtrans$P8) == 0)){P8[[i]]<-0}else {P8[[i]]<-CCtrans$P8[[5]]}  
   if((length(CCtrans$P9) == 0)){P9[[i]]<-0}else {P9[[i]]<-CCtrans$P9[[5]]}  
   if((length(CCtrans$P10) == 0)){P10[[i]]<-0}else {P10[[i]]<-CCtrans$P10[[5]]}  
   if((length(CCtrans$P11) == 0)){P11[[i]]<-0}else {P11[[i]]<-CCtrans$P11[[5]]}  
   if((length(CCtrans$P12) == 0)){P12[[i]]<-0}else {P12[[i]]<-CCtrans$P12[[5]]}  
   if((length(CCtrans$P13) == 0)){P13[[i]]<-0}else {P13[[i]]<-CCtrans$P13[[5]]}  
   if((length(CCtrans$P14) == 0)){P14[[i]]<-0}else {P14[[i]]<-CCtrans$P14[[5]]}  
   if((length(CCtrans$P15) == 0)){P15[[i]]<-0}else {P15[[i]]<-CCtrans$P15[[5]]}  
   
  ## compile parameters to predefined lists
   iterate[[i]] <- i
   TH_Tree[[i]] <- thTree
   EXT[[i]] <- ext
   TOL[[i]] <- tol
 
}
ws_tune <- as.data.frame(cbind(iterate, TH_Tree, EXT, TOL,Acc, nMC, nAC,
                                     P1,P2,P3,P4,P5,P6,P7,P8,P9,
                                     P10,P11,P12,P13,P14,P15))
                                    

write.csv(ws_tune, "simpleWS_tune_300buf.csv")


