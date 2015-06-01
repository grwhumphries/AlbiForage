

source("C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/R/Resampling_WAAL_data_April27/Predict_Functions.R")


Predictors<-c("SWH","Bathy","Dist_to_coast","SST","WindDifferential","WindDirToFlight","Period","T_before","T_after",
              "Sp_before","WD_before","WDTF_before","HeadWind","ResT_before","ResT_after","ResT_during","Approach","Depart","NumTurnsBefore",
              "NumTurnsAfter","MajTurnsBefore","SPslope","Fnspd","Sp_around_event","torSeg1","AccelSeg1","Pathchange1","PathChange_headwind1",
              "WinDifSeg1","torSeg2","AccelSeg2","Pathchange2","PathChange_headwind2","WinDifSeg2")


WS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/Bird_Tracks_Resampled_April27_2015/Summarized_data/"
ModWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/Bird_Tracks_Resampled_April27_2015/ModelData/"
PlotWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/Bird_Tracks_Resampled_April27_2015/TEMP/"
Species<-"WAAL"

Nmodels<-5
samples<-200
Nclus<-3
RF.mtry<-3
RF.ntree<-1000
range<-10




Create.Train(WS,ModWS,Species,Nmodels,samples,Nclus,RF.mtry,RF.ntree,Predictors,PLOT=T,OUTPUT=T,verbose=TRUE,threed=TRUE)

filenames=list.files(path=WS, pattern="*.txt",full.names=TRUE)                                 
Tracks = lapply(filenames, function(x){read.table(file=x,sep=",",header=T)})

for(Track in Tracks){

  Id<-as.character(Track$BirdName[1])
  df<-cluslab.to.track(Track,Id,ModWS,verbose=TRUE)
  dfprd<-RF.preds(df,Id,ModWS,RF.mtry,RF.ntree,Predictors,CV=FALSE,verbose=TRUE)  
  prdcor<-Pred.Correct(dfprd,write.out=TRUE,Id,PlotWS,verbose=TRUE)
  SummaryTable<-Summary.table(prdcor, range)
  Assessment<-assess.model(SummaryTable, Id)
  print(Assessment)

}


WS<-PlotWS
Plot.track(plot.all=FALSE,PlotWS,interactive=TRUE,Sit.Fly=FALSE,verbose=TRUE)







Track<-read.table(file="WAAL_12503_Resampled_120_Output.txt",sep=",",header=T)
ModWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/Bird_Tracks_Resampled_April27_2015/ModelData/"
setwd(WS)
Id<-"Bird_12503"

#### Add the Importance and variable trimming function
#### Incubating v Brooding birds









