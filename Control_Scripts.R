## Test2

source("C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/AlbiForage/Predict_Functions.R")


Predictors<-c("SWH","Bathy","Dist_to_coast","SST","WindDifferential","WindDirToFlight","Period","T_before","T_after",
              "Sp_before","WD_before","WDTF_before","HeadWind","ResT_before","ResT_after","ResT_during","Approach","Depart","NumTurnsBefore",
              "NumTurnsAfter","MajTurnsBefore","SPslope","Fnspd","Sp_around_event","torSeg1","AccelSeg1","Pathchange1","PathChange_headwind1",
              "WinDifSeg1","torSeg2","AccelSeg2","Pathchange2","PathChange_headwind2","WinDifSeg2")


#Predictors<-c("ResT_during","torSeg2","ResT_before","Fnspd")



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

filenames=list.files(path=WS, pattern="*.txt",full.names=TRUE)                                 
Tracks = lapply(filenames, function(x){read.table(file=x,sep=",",header=T)})


zz<-1
while(length(Predictors)>3){
  print(zz)
  
  Create.Train(WS,ModWS,Species,Nmodels,samples,Nclus,RF.mtry,RF.ntree,Predictors,PLOT=T,OUTPUT=T,verbose=TRUE,threed=FALSE)
  
  v.Imps<-VarImps(ModWS,RF.mtry,RF.ntree,Predictors,verbose=TRUE)
  
  assess.matrix<-matrix(ncol=5,nrow=0)
  colnames(assess.matrix)<-c("Prop.corr","Total.Over","Beh.Type","Bird","Tdiff")

  for(Track in Tracks){
    

    Id<-as.character(Track$BirdName[1])
    df<-cluslab.to.track(Track,Id,ModWS,verbose=TRUE)
    dfprd<-RF.preds(df,Id,ModWS,RF.mtry,RF.ntree,Predictors,CV=FALSE,verbose=F)  
    prdcor<-Pred.Correct(dfprd,write.out=F,Id,PlotWS,verbose=FALSE)
    
    SummaryTable<-Summary.table(prdcor, range)
    print(SummaryTable)
    
    Assessment<-assess.model(SummaryTable, Id)
    assess.matrix<-rbind(assess.matrix,Assessment)
  }

  matrix.name<-paste(PlotWS,"Output_assessment/assess_",zz,".csv",sep="")
  write.csv(assess.matrix,matrix.name,row.names=F)

  pred.name<-paste(PlotWS,"Output_assessment/predictors_",zz,".txt",sep="")
  write(Predictors,pred.name)
  
  Predictors<-as.vector(v.Imps$Predictors)[1:(nrow(v.Imps)-2)]
  zz<-zz+1
}






WS<-PlotWS
Plot.track(plot.all=FALSE,PlotWS,interactive=TRUE,Sit.Fly=FALSE,verbose=TRUE)



Track<-read.table(file="WAAL_1403_Resampled_120_Output.txt",sep=",",header=T)
ModWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/Bird_Tracks_Resampled_April27_2015/ModelData/"
setwd(WS)
Id<-"Bird_4040403A"


#### Incubating v Brooding birds









