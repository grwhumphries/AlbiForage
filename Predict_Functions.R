####################################################################################
#### This file contains all the functions used for making predictions on tracks ####
####################################################################################

## These functions were written by Grant Humphries.  May, 2015
## Department of Neurobiology, Physiology and Behavior
## Lab of Dr. Gabrielle Nevitt
## For the prediction of foraging events using random forests clustering 


#############################################################################################
## Load libraries ##

library(randomForest)
library(iplots)
library(ggplot2)
library(ggvis)
library(dplyr)
library(rgl)
source("C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/RF_Clustering.R")
##############################################################################################



################################################################
#### FUNCTION FOR CREATION OF TRAINING DATA FOR ANY SPECIES ####
################################################################

Create.Train<-function(WS,OutWS,Species,Nmodels,samples,Nclus,RF.mtry,RF.ntree,Predictors,PLOT=F,OUTPUT=T,verbose=TRUE,threed=F){

#########################################################################################################################
# WS          = Workspace where tracks are stored                       (Character string)
# OutWS       = Workspace where output training data will be stored     (Character string)
# Species     = Species name                                            (Character string)
# Nmodels     = Number of training data dataframes to create            (Integer)
# samples     = Number of "non events" to sample at each iteration      (Integer)
# Nclus       = Number of clusters to identify using pamNew             (Integer)
# RF.mtry     = Random Forests MTRY setting                             (Integer)
# RF.ntree    = Random Forests maximum trees setting                    (Integer)
# Predictors  = A vector of predictors to use in the model              (Vector of strings)
# PLOT        = A boolean: output plot = T, no plot = F                 (boolean)
# OUTPUT      = A boolean: output training data to file?                (boolean)
# threed      = A boolean: do we output a 3d scatter plot               (boolean)
                
#########################################
# Variable definitions in function
########################################

#### filenames = List of all filenames with .txt at the end
#### datalist = All the files in filenames opened as a single dataframe
#### A = The foraging events selected from the entire datalist
#### B = The dataset with only the predictor variables from "Predictors" vector plus identification variables
#### C = The dataset, B, with only complete cases selected
#### D = Dataframe with identification information to be appended after
#### E = Dataframe D bound to the RF cluster labels assigned by the RFdist function, and the Row index for merging together
#### f = Merged dataset E to A - puts clusterlabels in the correct positions to create training data
#### G = List of numbers from 1 to Nmodels: represents the number of times to randomly sample for training data
#### H = All "non-events" from the original datalist
#### J = Randomly sampled "non-events" from dataframe H
#### K = Non-events and Events combined into single dataframe
#### L = Output file name for the modeling data

##########################################################################################################################

  if(verbose) cat("Creating training data","\n","\n")

  
  filenames=list.files(path=WS, pattern="*.txt",full.names=TRUE)                                 ## read in the file names from the path as long as they are a .txt ##
  datalist = do.call(rbind,lapply(filenames, function(x){read.table(file=x,sep=",",header=T)}))  ## Rbind all the files together as a large data table ##
  
  A<-tbl_df(datalist[which(datalist$X.event==1),])                                               ## Select the foraging events and define a row index to match back the RF cluster labels
  A$RowIndex<-c(1:nrow(A))                                                                       
    
  B<-tbl_df(data.frame(A[,Predictors],RowIndex=A$RowIndex))                                      ## Select only the predictors from the list above and find all the complete cases (RF won't run with missing data)
  C<-B[complete.cases(B),]                                                                       ## We append the index/ID values 
  
  D<-data.frame(Species=rep(Species,nrow(C)),RowIndex=C$RowIndex)                                ## Create a dataframe that has species name and Indices that can be appended 
  
  C<-dplyr::select(C,-RowIndex)                                                                  ## Remove the RowIndex for modeling
  
  ##############################################################################################################
  
  if(verbose) cat("Running the RF clustering algorithm, this may take a moment...","\n","\n")    #### This runs the clustering algorithm 
  distRF<-RFdist(C, mtry1=RF.mtry, RF.ntree, 50, addcl1=T,addcl2=F,imp=T, oob.prox1=T)
  
  cmd1<-cmdscale(as.dist(distRF$cl1),3)                                                          ### Setting the value of 3 means we can get three axes for a 3d figure if we want
  RFclusterLabel = pamNew(cmd1, Nclus)
  if(verbose) cat("Done...","\n","\n")

  E<-data.frame(D,clusterlab=RFclusterLabel)                                                     ### Create the data frame that has the events classified, then merge it
  f<-merge(A,E,by="RowIndex")
  
  f<-dplyr::select(f,-RowIndex)                                                                  ### Remove the "RowIndex" label
  
  ##############################################################################################################
  ##### Do we create output files?
  if(OUTPUT==T){
    
    if(verbose) cat("Output training data being created, located in: \n",OutWS,"\n","\n")
        
    G<-c(1:Nmodels)                                                                              ### Now randomly select non-events, but there are a lot of non events, so we do this "Nmodels" number of times
    
    for(i in G){    
      
      if(verbose) cat("Creating training data: ",i,"\n")

      H<-tbl_df(datalist[which(datalist$X.event==0),])                                           ## select all non-events
      J<-sample_n(H,samples,replace=F)                                                           ## randomly sample non-events
      J$Species<-Species                                                                         ## Assign species name column
      J$clusterlab<-as.factor(0)                                                                 ## Assign a factor value of 0 here to delineate a "non event"
      
      K<-rbind(f,J)
      L<-paste(OutWS,"/ModelData_",Species,"_",i,".txt",sep="")
      
      write.table(K,L,sep=",",row.names=F)    
      
      
    }
  }else{if(verbose) cat("\n","warning: OUTPUT set to FALSE, no training data being created","/n")}
  
  ###### Do we plot the three dimensional output?
  if(threed==T){
    plot3d(cmd1[,1:3],col=RFclusterLabel,type="s",size=1,box=F,xlab="Scaling Dimension 1",ylab="Scaling Dimension 2", zlab = "Scaling Dimension 3")
  }


  ######## Do we plot the output of the clusters? 
  if(PLOT==T){
    
    jpegnm<-paste(OutWS,"/clusters.jpg",sep="")
    
    if(verbose) cat("\n","Writing plot: \n",jpegnm,"\n")
        
    jpeg(filename = jpegnm,width = 600,height=600,units="px")
    plot(cmd1,type="n",xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
    text(cmd1, label=RFclusterLabel,col=RFclusterLabel)
    dev.off()
    
    
        
  }else{if(verbose) cat("Plotting set to FALSE, no output plot created","/n")}

if(verbose) cat("complete...")

}



#######################################################################
#### GET MEAN VARIABLE IMPORTANCES FROM RF MODELS FOR VAR TRIMMING ####
#######################################################################

### Will want to put plotting functions in here to get out variable plots ###

VarImps<-function(WS,RF.mtry,RF.ntree,Predictors,verbose=TRUE){
  #########################################################################################################################
  # WS          = Workspace where Model training data are located                 (Character string)
  # Predictors  = Vector of predictors to be used for modeling                    (Vector)
  # CV          = If TRUE, will do predictions using independent subsets          (Boolean)
  # verbose     = If TRUE, will output iteratively                                (Boolean)
  # RF.mtry     = Random Forests MTRY setting                                     (Integer)
  # RF.ntree    = Random Forests maximum trees setting                            (Integer)
  
  
  #####################################################
  # Variable definitions in function
  #####################################################
  #### Modlist    = List of full path names from the training data files
  #### OUTPUT     = matrix that will be returned in this function
  #### Mod        = Iteratively drawn Model from the ModList
  #### colname    = We use this variable to define the column names where predictions are stored
  #### A          = The filename "Mod" opened as a data frame
  #### B          = Complete cases of data frame A
  #### C          = Data frame with ID variables for re-merging
  #### rf1        = random forests model object
  #### OUT        = Track output with associated predictions
  
  if(verbose){cat("--------------------------------------------------------------","\n")
              cat("Calculating variable importances across training data samples", "\n")           
  }
  
  
  set.seed(248)
  
  Modlist<-list.files(path=WS,pattern="*.txt",full.names=TRUE)
  rfimps<-matrix(nrow=length(Predictors),ncol=0)
  
  
  colname<-1
  for(Mod in Modlist){
    
    if(verbose) cat("running model using ",Mod,"\n")
    
    A<-tbl_df(read.table(Mod,sep=",",header=T))                                                             # Read the table Mod  
    A<-tbl_df(data.frame(A[,Predictors],clusterlab=A$clusterlab,Eindex=A$Eindex,BirdName=A$BirdName))       # Select the "Predictors" columns
    A$RowIndices<-c(1:nrow(A))                                                                              # Set row-indices so the data can be matched back to the tracks
    
    if(length(which(names(A)=="torSeg1"))==1){                                                              #### The variables torSeg1 and torSeg2 sometimes have problems with INF values...
      A<-A[is.finite(A$torSeg1),]
    }
    
    if(length(which(names(A)=="torSeg2"))==1){
      A<-A[is.finite(A$torSeg2),]
    }
    
    B<-A[complete.cases(A),]
    C<-data.frame(clusterlab=B$clusterlab,Eindex=B$Eindex,RowIndices=B$RowIndices,BirdName=B$BirdName)      # Create a dataframe with the indices, cluster labels, and row indices
    B<-data.frame(dplyr::select(B,-Eindex,-RowIndices,-BirdName))
    
    rf1 <- randomForest(as.factor(clusterlab)~.,mtry=RF.mtry,ntree=RF.ntree,data = B,importance=T)          #### B is the training data for the model
    rfIp<-data.frame(importance(rf1)[,"MeanDecreaseGini"])                                                  # Gets variable importances for trimming, etc.
    rfimps<-cbind(rfimps,rfIp[,1])
    
    
  }
  
  rfimp.mean<-rowMeans(rfimps)                                                                              # Calculate the mean importance
  rfimp.frame<-tbl_df(data.frame(Predictors=Predictors,MeanDecreaseGini=rfimp.mean))
  
  rfImp<-arrange(rfimp.frame,-MeanDecreaseGini)                                                             # Sorts the dataframe from most important to least
  
  
  return(rfImp)
}



#######################################################################################
#### THIS WILL MATCH CLUSTER LABELS FROM MODEL DATA (create.train()) TO THE TRACKS ####
#######################################################################################


cluslab.to.track<-function(Track,Id,WS,verbose=TRUE){
#########################################################################################################################
# Track       = Dataframe of track that will have predictions                   (Data frame)
# Id          = Id of the bird being predicted                                  (Character string)
# WS          = This is the workspace of the ModelData (training data)          (Character string)

#####################################################
# Variable definitions in function
#####################################################
#### Modlist    = List of full path names from the training data files
#### A          = First training data file - we are only interested in where events occur (not the non-events)
#### B          = The subsetted cluster labels > 0, and ID characteristics
#### C          = The subset of the cluster labels, which belong to bird Id

  if(verbose) cat("Matching cluster labels with Track for ",Id, "\n")  

  Modlist<-list.files(path=WS,pattern="*.txt",full.names=TRUE)
  ### We read the very first Modlist file because the clusterlab data are all the same (only non-events) vary between each training file
  A<-tbl_df(read.table(Modlist[[1]],sep=",",header=T))  
  ### Now subset out cluster labels > 0
  A<-filter(A,clusterlab>0)
  B<-tbl_df(data.frame(Index=A$Eindex,BirdName=A$BirdName,Type=A$clusterlab))
  
  ### Select values from B, where BirdName = the Id (bird of interest)
  C<-filter(B,BirdName==Id)

  Track$clusterlab<-"0"
  
  ### Go through each row of the events and append them in the correct location in "Track"

  for(jj in 1:nrow(C)){
    
    #C$Index[jj] is the location (row) of the track, where that event occurred for row jj of the dataframe "C"
    Track$clusterlab[C$Index[jj]]<-C$Type[jj]
    
  }
  
  if(verbose) cat(Id,"complete","\n")

  return(Track)
  
}


###################################################################################
#### THIS FUNCTION CREATES PREDICTIONS ON THE TRACK USED IN CLUSLAB.TO.TRACK() ####
###################################################################################


RF.preds<-function(Track,Id,WS,RF.mtry,RF.ntree,Predictors,CV=TRUE,verbose=TRUE){
#########################################################################################################################
# Track       = This is the output from the function cluslab.to.track()         (Data frame)
# Id          = Id of the bird being predicted                                  (Character string)
# WS          = Workspace where Model training data are located                 (Character string)
# Predictors  = Vector of predictors to be used for modeling                    (Vector)
# CV          = If TRUE, will do predictions using independent subsets          (Boolean)
# verbose     = If TRUE, will output iteratively                                (Boolean)
# RF.mtry     = Random Forests MTRY setting                                     (Integer)
# RF.ntree    = Random Forests maximum trees setting                            (Integer)


#####################################################
# Variable definitions in function
#####################################################
#### Modlist    = List of full path names from the training data files
#### OUTPUT     = matrix that will be returned in this function
#### Mod        = Iteratively drawn Model from the ModList
#### colname    = We use this variable to define the column names where predictions are stored
#### A          = The filename "Mod" opened as a data frame
#### B          = Complete cases of data frame A
#### C          = Data frame with ID variables for re-merging
#### rf1        = random forests model object
#### OUT        = Track output with associated predictions

  if(verbose){cat("--------------------------------------------------------------","\n")
    cat("Creating predictions for bird ",Id,"\n")           
  }


  set.seed(248)
  
  Modlist<-list.files(path=WS,pattern="*.txt",full.names=TRUE)
  OUTPUT<-matrix(nrow=nrow(Track),ncol=0)

  colname<-1
  for(Mod in Modlist){
    
    if(verbose) cat("running model using ",Mod,"\n")
  
    A<-tbl_df(read.table(Mod,sep=",",header=T))                                                             # Read the table Mod
  
    A<-tbl_df(data.frame(A[,Predictors],clusterlab=A$clusterlab,Eindex=A$Eindex,BirdName=A$BirdName))       # Select the "Predictors" columns
    A$RowIndices<-c(1:nrow(A))                                                                              # Set row-indices so the data can be matched back to the tracks
    
      
    if(length(which(names(A)=="torSeg1"))==1){                                                              #### The variables torSeg1 and torSeg2 sometimes have problems with INF values...
      A<-A[is.finite(A$torSeg1),]
    }
    
    if(length(which(names(A)=="torSeg2"))==1){
      A<-A[is.finite(A$torSeg2),]
    }
    
    #### If we want to use cross validation (across individual birds), instead of predicting back to itself
    if(CV){
      
      B<-filter(A,BirdName!=Id)                                                                             # Selects all the data not associated with the bird we are predicting to
      B<-B[complete.cases(B),]
    }else{B<-A[complete.cases(A),]}
        
    C<-data.frame(clusterlab=B$clusterlab,Eindex=B$Eindex,RowIndices=B$RowIndices,BirdName=B$BirdName)      # Create a dataframe with the indices, cluster labels, and row indices
    B<-data.frame(dplyr::select(B,-Eindex,-RowIndices,-BirdName))
    
          
    rf1 <- randomForest(as.factor(clusterlab)~.,mtry=RF.mtry,ntree=RF.ntree,data = B,importance=T)          #### B is the training data for the model
            
    prds<-predict(rf1,Track)                                                                                #### Now we predict back to the Track 
    OUTPUT<-cbind(OUTPUT,prds)
    colnames(OUTPUT)[ncol(OUTPUT)]<-paste("prds",colname,sep="")
    colname<-colname+1
  }
  
  OUT<-cbind(Track,OUTPUT)
  
  if(verbose){cat("model predictions successful","\n")
              cat("--------------------------------------------------------------","\n")
  }

  return(OUT)
}


####################################################################################################
#### THIS FUNCTION WILL TAKE ALL THE PREDICTION COLUMNS AND THEN CORRECT FOR "OVER PREDICTIONS" ####
####################################################################################################


Pred.Correct<-function(Track,write.out=TRUE,Id=NA,WS=NA,verbose=TRUE){
#########################################################################################################################
# Track       = Dataframe of track and predictions: output from RF.preds()                   (Data frame)
# write.out   = If TRUE, the Track will be written out to a csv file with predictions        (Boolean)
                # This is required if plotting the outputs
# Id          = Id of the track, used to create file name for output                         (Character string)
# WS          = Output workspace where output is created                                     (Character string)

#####################################################
# Variable definitions in function
#####################################################
#### A        = All prediction columns
#### B        = Row means of all the predictions
#### C        = Column with predictions that have 100% certainty (same across all prd columns) 
#### D        = Run Length Encoding data - calculates where repeated sequences occur in a vector
#### E        = The data from D organized into a dataframe for access
#### f        = Empty vector of 0s that we use to assign predictions - this is the output of this function
#### G        = Output filepath for Track when write.out=T



  ### Selects all the columns with "prd" in the name (prediction columns), subtracts 1 (predictions are made from 1+, but need to start at 0)
  A<-Track[,grep("prd",names(Track))] -1
  
  B<-rowMeans(A,na.rm=TRUE)
  ### In case some NANs are formed - we convert those to 0s
  B[is.nan(B)]<-0
  
  ### This removes predicted points that have uncertainty associated with them
  ### essentially, the predictions must agree across all model runs
  ### otherwise it is given a score of 0
  C<-sapply(B,function(x){if(x-floor(x)==0){x=x}else{x=0}})
  
  ### Lets hard wire in that the final point has a score of 0 in order to avoid infinite loops!
  C[length(C)]<-0
  
  
  ### This next chunk of code will go through the corrected predictions and then find 
  ### Repeats that occur in sequence
  ### rle finds sequences of numbers that repeat in a vector
  D<-rle(C)
  
  ### bind the values and lengths together, then add the row names (which corresponds to index values in the Track dataframe (A))
  E<-data.frame(D$values,D$length)
  E$rows<-as.numeric(rownames(E))
  
  ### Now we create an empty vector with all "0"s that is the same length of C in which to place our predictions
  f<-rep(0,length(C))
  
  ### Now you can loop through the dataframe E and get the index values where repeated values occur and then place them 
  
  for(i in 1:nrow(E)){
    ## If D.values (i.e. the value of clusterlabel) is greater than 1, AND it repeats more than 2 times
    if(E$D.values[i]>0 & E$D.length[i]>2){
      ## the value of E$rows[i-1] is the index BEFORE the sequence starts repeating.. so we add 1 to get the starting point
      Start<-E$rows[i-1]+1
      ## The number of times it repeats is in the value E$D.length[i]
      Len<-E$D.length[i]
      
      ## We don't want to place the prediction right at the start, nor in the middle
      ## as the behavior usually takes place somewhere around the beginning of these sequences
      ## we will place them at the "quarter way" mark (i.e 1/4th of the way along the sequence)
      
      qrtpoint<-Start+(floor((Len)/4))
      
      ## We take the value of the cluster label and place it at the quarter way mark of the sequence
      f[qrtpoint]<-E$D.value[i]
    }
    
  }
  
  Track$predicted<-f
  
  ## If write.out == true, then write the output
  if(write.out==TRUE){
    
    if(verbose) cat("write.out set to TRUE. Writing output for",Id,"\n","\n")
    
    G<-paste(WS,"/",Id,".txt",sep="")
    write.table(Track,G,sep=",")
    
  }
  
return(Track)

}

###################################################################################################
#### THIS FUNCTION CREATES THE SUMMARY TABLE THAT WILL BE USED TO CALCULATE SUMMARY STATISTICS ####
###################################################################################################


Summary.table<-function(Track, range){
#########################################################################################################################
# Track       = Dataframe of track and predictions: output from Pred.Correct()                    (Data frame)
# range       = error range for predictions by number of points from the prediction               (Integer)  
  
####################################################
# Variable definitions in function
#####################################################
#### A        = All prediction columns
#### B        = Matrix that matches up cluster labels (observed events) with predicted events
#### C        = Second - First datetimes to get track resolution for calculating differences between obs and predicted
#### i        = index value in for loop for nrows of A
#### st       = starting index (based on the range variable) for which to search for a match with predicted value
#### end      = ending index (based on range variable) for which to search for a match with predicted value
#### Match    = The Observed values between st and end which match the Predicted value
#### tdiff    = If there is a Match, then the time between observed and predicted is calculated
#### Btyp     = Value of the Prediction (behavior type / cluster label)
#### T.evs    = total number of events in the whole dataset which match Btyp 
#### rw       = the row that will be combined into the matrix, B

  A<-tbl_df(data.frame(Observed=Track$clusterlab,Predicted=Track$predicted))
  
  ###### THIS bit of code will match up the predictions to observed events by the behavior types   
  B<-matrix(ncol=3,nrow=0)
  colnames(B)<-c("Btyp","total.evs","tdiff")
  
  # This loop will search through the dataframe "A" by the Predicted values
  # When it comes across a predicted value < 0 (i.e.  a foraging event)
  # It will then search through the Observed values within a range (see variable, range) of points/distance from the event
  # If the predicted matches the observed value (within the range), it calculates the distance (in time) between predicted and observed
  
  #### This is the track resolution, used to calculate time between events
  C<-as.POSIXct(Track$LocalDateTime[2])-as.POSIXct(Track$LocalDateTime[1])
  
  for(i in 1:nrow(A)){
    
    if(A$Predicted[i] > 0){
      end<- i+range
      if(i < range){st=1}else{st=i-range}                # If there is a prediction less than the range value, start search changes to row 1
      if(end > nrow(A)){end<-nrow(A)}                    # If there is a predction closer to the end than the range value , end search changes to last row
      
      Match<-which(A$Observed[st:end]==A$Predicted[i])   # Calculate which events match predictions +/- range index values
      
      tdiff<-abs(Match - range)*C                        # time difference between prediction and event
      if(length(tdiff)==0){tdiff=NA}                     # If tdiff cannot be calculated (no matching observed and predicted), set it to NA
      
      Btyp<-A$Predicted[i]
      t.evs<-length(which(A$Observed==A$Predicted[i]))
      rw<-data.frame(Btyp,t.evs,tdiff)
      B<-rbind(B,rw)
    } 
  }
  
  return(B)

}



assess.model<-function(SummaryTable, Id){
#########################################################################################################################
# SummaryTable      = Dataframe output from Summary.table()                           (Data frame)
# Id                = Id of the bird being predicted                                  (Character string)                 
  
####################################################
# Variable definitions in function
#####################################################
#### A        = Matrix where summary data will be stored
#### B        = Matrix where bird name (Id) will be stored 
#### j        = index j for for loop
#### C        = subset of SummaryTable for the unique value of "Btyp" - i.e. the cluster type of interest
#### D        = Total number of rows
#### E        = Total number of correct predictions for a cluster label (Btyp)
#### f        = Proportion of correctly classified cluster labels
#### G        = Total number of predictions of a specific cluster label overall
#### H        = Total number of over-predicted points
#### k        = Mean time differences

  A<-matrix(ncol=6,nrow=0)
  B<-matrix(ncol=1,nrow=0)
  
  for(j in unique(SummaryTable[,"Btyp"])){
    
    C<-subset(SummaryTable,SummaryTable[,"Btyp"]==j)       ## First, subset for each unique value of "Btyp" (cluster label)
    
    D<-C[1,"t.evs"]                                        ## We pull the first value for total number of events (i.e. total number of observed events in the Track)
      
    E<-length(which(!is.na(C[,"tdiff"])))                  ## This is the total number of times that the value Btyp was correctly predicted in the track
        
    f<-E/D * 100                                           ## Proportion of correctly classified events of specific cluster label
      
    G<-nrow(C)                                             ## Represents the total number of predictions made by the algorithm
      
    H<-G-E                                                 ## Total number of over-predictions from particular Btyp
    
    k<-mean(C[,"tdiff"],na.rm=TRUE)                        ## Mean time difference of predictions v observed
    
    rw<-data.frame(f,D,G,H,j,as.numeric(k))                ## put proportion correct, total number of over predictions and Btyp in a row then append to matrix A
  
    A<-rbind(A,rw)
    B<-rbind(B,Id)                                         ## in order to prevent formation of a character matrix, bind names to another one
  }
  
  A<-data.frame(A)
  names(A)<-c("prop.corr","Total.observed","Total.preds","Total.Over","Beh.Type","mean.tdiff")
  A$Bird<-B[,1]
  
  return(A)
  
}


##########################################################################################
#### THIS FUNCTION IS FOR PLOTTING THE TRACKS  - CAN BE DONE INDIVIDUALLY, OR IN BULK ####
#### THERE IS FUNCTIONALITY HERE FOR INTERACTIVE PLOTTING AS WELL TO VIEW THE TRACKS  ####
##########################################################################################


Plot.track<-function(plot.all=TRUE,WS,interactive=FALSE,Sit.Fly=FALSE,verbose=TRUE){
#########################################################################################################################
# plot.all        = If true, plots for all tracks in WS will be plotted                                           (boolean)
# WS              = WS where tracks are located for plotting                                                      (Character string)
# interactive     = If true, user will be able to select an individual bird and then plot it interactively        (boolean)
# Sit.Fly         = If true, the values of Sitting v Flying will be plotted on the map                            (boolean)

  
####################################################
# Variable definitions in function
#####################################################
#### Date         = Current date and time for creation of new output workspace
#### outputspace  = Newly created workspace
#### A            = List of file names in workspace where new files are output
#### B            = List of all files in A open as dataframes
#### C            = Index of B in for loop for plotting all birds (i.e., individual tracks)
#### C.ID         = Name of the track being plotted
#### D            = Filtering every dataframe for ease of plotting
#### E            = All rows where foraging events are predicted
#### f            = dataframe E bound to D with newly formed column which represents predictions and observed values
#### Gp           = The initial ggplot 
#### G            = The ggplot layers
#### H            = List of files in workspace
#### I            = All the tracks in list H as a single data frame
#### J            = Subsetted track as selected by user
#### K            = Subsetted predicted points for plotting
#### Z            = Control variable in while loop for interactive plots

### First we write an error catch in case plot.all=T and interactive=T as they cannot both coincide.

  if(plot.all==TRUE & interactive==TRUE) stop("plot.all and interactive arguments cannot both be set to TRUE")
  if(plot.all==FALSE & interactive==FALSE) stop("both plotting arguments are set to FALSE, one must be set to TRUE")

  if(plot.all==TRUE){
    ## This calculates the current date/time and then collapses it with no special characters for creating output workspace for plots
    Date<-gsub(":","",gsub(" ","",date()))
    outputspace<-paste(WS,"/",Date,"/",sep="")
    dir.create(outputspace)
    
    if(verbose){
      cat("Plotting all set to TRUE, plotting all tracks in directory ",WS, "\n")
      cat("creating new directory", outputspace,"\n")
      cat("--------------------------------------------------------------------------------","\n")
    }
    
    ### List all the files in the workspace
    A<-list.files(path=WS, pattern="*.txt",full.names=TRUE)
    ### Then we use lapply to open them all simultaneously so we don't have to open them one at a time
    B<-lapply(A, function(x){read.table(file=x,sep=",",header=T)})
    
    ### Now we are looping through all the dataframes in B
    for(C in B){
     
      ## this creates the ID of the track for plotting and saving
      C.ID<-as.character(C$BirdName[1])
            
      if(verbose) cat("Creating output plot for",C.ID,"\n")
      ## Select the data we need from C
      D<-data.frame(sitfly=as.factor(C$Sit.Fly),observed=C$clusterlab,predicted=C$predicted,x=C$X.Longitude,y=C$X.Latitude)
      D$Index<-D$observed
      
      ## This selects all the points where predicted >0 (i.e. all the predicted foraging events)
      E<-tbl_df(D) %>% filter(predicted>0)    #,Sit.Fly=="Fly"
      
      ## We add 10 so that we can place them in a dataframe for plotting - placing into the column "Index" 
      E$Index<-E$predicted+10
      
      ## We rbind the two frames here  and then convert the newly formed Index column to a factor
      f<-rbind(D,E)
      f$Index<-factor(f$Index)
      
      ###### Calling in the ggplot here
      Gp<-ggplot(data=f,aes(x=x,y=y))
      
      G<-Gp+ggtitle(C.ID)+
        
        geom_point(data=f,aes(x=x,y=y,color=Index,size=Index,shape=Index))+ 

        scale_color_manual(values=c("0"="black","1"="orange","2"="blue","3"="purple","4"= "red",
                                    "11" = "orange", "12" = "blue", "13" = "purple", "14" = "red"),
                           
                           labels=c("0"="Not foraging","1"="foraging type 1","2"="foraging type 2","3"="foraging type 3","4"= "foraging type 4",
                                    "11" = "predicted 1", "12" = "predicted 2", "13" = "predicted 3", "14" = "predicted 4"))+                   
                           
                      
        scale_size_manual(values=c("0"=0.5,"1"=4,"2"=4,"3"=4,"4"= 4,"11"=7 , "12"=7, "13"=7, "14"=7),
                          
                          labels=c("0"="Not foraging","1"="foraging type 1","2"="foraging type 2","3"="foraging type 3","4"= "foraging type 4",
                                   "11" = "predicted 1", "12" = "predicted 2", "13" = "predicted 3", "14" = "predicted 4"))+
        
        
        scale_shape_manual(values=c("0"=20,"1"=15,"2"=16,"3"=17,"4"= 18,"11"=0 , "12"=1, "13"=2, "14"=5),
                           
                           labels=c("0"="Not foraging","1"="foraging type 1","2"="foraging type 2","3"="foraging type 3","4"= "foraging type 4",
                                    "11" = "predicted 1", "12" = "predicted 2", "13" = "predicted 3", "14" = "predicted 4"))+     
              
        #geom_point(data=Trck.sit,aes(x=X.Longitude,y=X.Latitude),color="black",size=0.5)+
        theme(axis.line=element_blank(),
              axis.ticks=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              title = element_text(family="serif",color="black",size=12),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title=element_blank(),
              legend.background = element_rect(color="black",linetype=2),
              legend.key=element_blank())
      
      G
      
      filename<-paste(outputspace,C.ID,".tiff",sep="")
      ggsave(filename,width=140,height=100,units="mm",dpi=600)
      
      
    }
  }


  if(interactive==TRUE){
  
    ### This function will retrieve all the values at the key value selected in the ggvis plot
    all_values <- function(x) {
      if(is.null(x)) return(NULL)
      row<- K[K$Key == x$Key,]
      
      paste0(names(row), ": ",as.character(row), collapse = "<br />")
    }
    
    
    ### This function was taken from: http://q2a.science/add-a-plot-title-to-ggvis-i249368.htm by user tonytonov
    ### It is a function to add the title to the interactive plot, which is not supported in ggvis as of May 2015
    add_title <- function(vis, ..., x_lab = "X units", title = "Plot Title") 
    {
      add_axis(vis, "x", title = x_lab) %>% 
        add_axis("x", orient = "top", ticks = 0, title = title,
                 properties = axis_props(
                   axis = list(stroke = "white"),
                   labels = list(fontSize = 0)
                 ), ...)
    }
    
    
    ### List all the files in the workspace
    H<-list.files(path=WS, pattern="*.txt",full.names=TRUE)
    ### Then we use lapply to open them all simultaneously so we don't have to open them one at a time
    I<-do.call(rbind,lapply(H, function(x){read.table(file=x,sep=",",header=T)}))
    
    
    ### We list all the available bird names in the directory
    cat("Below are all available bird tracks in the directory: \n")
    cat("-------------------------------------------------------------- \n")
    cat(as.vector(unique(I$BirdName)),"\n")
    cat("-------------------------------------------------------------- \n")
    
    ### Set up a while control loop for the user to go through and select plots
    ### will use Z as the control variable here. When set to F, the loop will end. 
    Z<-TRUE
    while(Z==TRUE){
            
      I.ID<-readline("Please type in the name of the track you are selecting: ")
      
      ### This filters all the datapoints for the track selected by the user
      J<-tbl_df(I) %>% filter(BirdName==I.ID)
    
      if(nrow(J)==0){cat("error: the track you typed in does not exist, please try again")}else{Z<-FALSE}

    } #while Z==TRUE
    
    
    ### This filters all rows where predictions are > 0 so we can visualize predicted points
    ### The key value is for the "all_values" function, which is used for the tooltip
    K<-filter(J,predicted>0)
    K$Key<-1:nrow(K)
    
    ### We want to show where the cluster labels (observed values) exist as well
    K1<-filter(J,clusterlab>0)
    
    
    ggvis() %>%
      layer_points(~X.Longitude,~X.Latitude,data=J,size := 10)%>%
      layer_points(~X.Longitude,~X.Latitude,data=K1,size := 200, shape=~factor(clusterlab),fill:="red") %>%
      layer_points(~X.Longitude,~X.Latitude,data=K, size := 100,fill=~factor(predicted), key := ~Key) %>%
      add_tooltip(all_values, "hover") %>%
      add_title(title=I.ID) %>%
      hide_legend("fill")
    
    
    
  } #if Interactive==TRUE
} #Function 
  
  
  


