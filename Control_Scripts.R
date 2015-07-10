## Test 4

source("C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/AlbiForage/Predict_Functions.R")


Predictors<-c("Period","T_before","T_after",
              "Sp_before","WD_before","WDTF_before","HeadWind","ResT_before","ResT_after","ResT_during","Approach","Depart","NumTurnsBefore",
              "NumTurnsAfter","MajTurnsBefore","SPslope","Fnspd","Sp_around_event","torSeg1","AccelSeg1","Pathchange1",
              "torSeg2","AccelSeg2","Pathchange2")

#"SWH","Bathy","Dist_to_coast","SST",
#Predictors<-c("ResT_during","torSeg2","ResT_before","Fnspd")
#"WinDifSeg1",,"WinDifSeg2","PathChange_headwind2""PathChange_headwind1""WindDifferential","WindDirToFlight"


WS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/analysis/resampled_120_sec/summarized/"
ModWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/analysis/resampled_120_sec/model_data/"
PlotWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/analysis/resampled_120_sec/plots/"
Species<-"WAAL"

Nmodels<-30
samples<-60
Nclus<-2
RF.mtry<-5
RF.ntree<-6000
range<-10

filenames=list.files(path=WS, pattern="*.txt",full.names=TRUE)                                 
Tracks = lapply(filenames, function(x){read.table(file=x,sep=",",header=T)})

Create.Train(WS,ModWS,Species,Nmodels,samples,Nclus,RF.mtry,RF.ntree,Predictors,PLOT=T,OUTPUT=T,verbose=TRUE,threed=TRUE)

v.Imps<-VarImps(ModWS,RF.mtry,RF.ntree,Predictors,verbose=TRUE)

assess.matrix<-matrix(ncol=5,nrow=0)
colnames(assess.matrix)<-c("Prop.corr","Total.Over","Beh.Type","Bird","Tdiff")

for(Track in Tracks){
  try({
    Id<-as.character(Track$BirdName[1])
    print(Id)
    
    df<-cluslab.to.track(Track,Id,ModWS,verbose=TRUE)
    dfprd<-RF.preds(df,Id,ModWS,RF.mtry,RF.ntree,Predictors,CV=FALSE,verbose=TRUE)  
    prdcor<-Pred.Correct(dfprd,write.out=T,Id,PlotWS,verbose=FALSE)
    
    SummaryTable<-Summary.table(prdcor, range)
    print(SummaryTable)
    
    Assessment<-assess.model(SummaryTable, Id)
    assess.matrix<-rbind(assess.matrix,Assessment)
  })
}


Plot.track(plot.all=TRUE,PlotWS,interactive=FALSE,Sit.Fly=FALSE,verbose=TRUE,shape.out=TRUE)


#### Incubating v Brooding birds




WS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/analysis/resampled_120_sec/summarized/"
ModWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/analysis/resampled_120_sec/model_data/"
PlotWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/analysis/resampled_120_sec/plots/"
Species<-"BBAL"

Nmodels<-30
samples<-60
Nclus<-2
RF.mtry<-5
RF.ntree<-6000
range<-10

filenames=list.files(path=WS, pattern="*.txt",full.names=TRUE)                                 
Tracks = lapply(filenames, function(x){read.table(file=x,sep=",",header=T)})



assess.matrix<-matrix(ncol=5,nrow=0)
colnames(assess.matrix)<-c("Prop.corr","Total.Over","Beh.Type","Bird","Tdiff")

for(Track in Tracks){
  try({
    Id<-as.character(Track$BirdName[1])
    print(Id)
    
    df<-cluslab.to.track(Track,Id,ModWS,verbose=TRUE)
    dfprd<-RF.preds(df,Id,ModWS,RF.mtry,RF.ntree,Predictors,CV=FALSE,verbose=TRUE)  
    prdcor<-Pred.Correct(dfprd,write.out=T,Id,PlotWS,verbose=FALSE)
    
    SummaryTable<-Summary.table(prdcor, range)
    print(SummaryTable)
    
    Assessment<-assess.model(SummaryTable, Id)
    assess.matrix<-rbind(assess.matrix,Assessment)
  })
}
Plot.track(plot.all=TRUE,PlotWS,interactive=FALSE,Sit.Fly=FALSE,verbose=TRUE,shape.out=TRUE)









############# THIS BIT IS FOR DOING THE PREDICTOR CLIPPING FOR VARIABLE IMPORTANCE #############
zz<-1
while(length(Predictors)>3){
  print(zz)
  
  Create.Train(WS,ModWS,Species,Nmodels,samples,Nclus,RF.mtry,RF.ntree,Predictors,PLOT=T,OUTPUT=T,verbose=TRUE,threed=TRUE)
  
  v.Imps<-VarImps(ModWS,RF.mtry,RF.ntree,Predictors,verbose=TRUE)
  
  assess.matrix<-matrix(ncol=5,nrow=0)
  colnames(assess.matrix)<-c("Prop.corr","Total.Over","Beh.Type","Bird","Tdiff")
  
  for(Track in Tracks){
    try({
      Id<-as.character(Track$BirdName[1])
      print(Id)
      
      df<-cluslab.to.track(Track,Id,ModWS,verbose=TRUE)
      dfprd<-RF.preds(df,Id,ModWS,RF.mtry,RF.ntree,Predictors,CV=FALSE,verbose=TRUE)  
      prdcor<-Pred.Correct(dfprd,write.out=T,Id,PlotWS,verbose=FALSE)
      
      SummaryTable<-Summary.table(prdcor, range)
      print(SummaryTable)
      
      Assessment<-assess.model(SummaryTable, Id)
      assess.matrix<-rbind(assess.matrix,Assessment)
    })
  }
  
  matrix.name<-paste(PlotWS,"Output_assessment/assess_",zz,".csv",sep="")
  write.csv(assess.matrix,matrix.name,row.names=F)
  
  pred.name<-paste(PlotWS,"Output_assessment/predictors_",zz,".txt",sep="")
  write(Predictors,pred.name)
  
  Predictors<-as.vector(v.Imps$Predictors)[1:(nrow(v.Imps)-2)]
  zz<-zz+1
}






