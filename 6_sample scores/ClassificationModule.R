.myClassificationBenchmarking=function(inputDataList,ncores=1,expClassName="proband",methodsList=list("reg","logReg","lda","qda","sqda","ridgeReg","lassoReg","ridgeLogReg","lassoLogReg","elasticNetLogReg","randomForest","boosting","bagging","neuralNet")){
  
  #inputDataList=resultsPosterior;ncores=1;expClassName="proband";methodsList=list("qda")
  x=unique(unlist(lapply(inputDataList,function(x){x$method})))
  PRresults=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(PRresults)=x
  colnames(PRresults)=methodsList
  PRresults=data.frame(PRresults)
  
  ROCresults=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(ROCresults)=x
  colnames(ROCresults)=methodsList
  ROCresults=data.frame(ROCresults)
  
  results9=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(results9)=x
  colnames(results9)=methodsList
  results9=data.frame(results9)
  
  results95=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(results95)=x
  colnames(results95)=methodsList
  results95=data.frame(results95)
  
  results85=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(results85)=x
  colnames(results85)=methodsList
  results85=data.frame(results85)
  
  perc9=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(perc9)=x
  colnames(perc9)=methodsList
  perc9=data.frame(perc9)
  
  perc95=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(perc95)=x
  colnames(perc95)=methodsList
  perc95=data.frame(perc95)
  
  perc85=matrix(NA,nrow=length(x),ncol=length(methodsList))
  row.names(perc85)=x
  colnames(perc85)=methodsList
  perc85=data.frame(perc85)
  
  ncores=min(ncores,length(inputDataList))
  
  
  x=list(seq(1,length(inputDataList)))
  res=lapply(x,.myClassificationBenchmarkDetailed,inputDataList,methodsList,expClassName)
  
  resROCdf=list()
  if(length(res)>0){
    
    if(length(res[[1]])>0){
    for(i in 1:length(res)){
      for(j in 1:length(res[[i]])){
        resROCdf=c(resROCdf,list(list(res[[i]][[j]]$ROCdf,method=paste0(res[[i]][[j]]$method,"_",res[[i]][[j]]$classifier))))
        PRresults[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$aucPR
        ROCresults[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$aucROC
        results85[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$res85
        results9[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$res9
        results95[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$res95
        perc85[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$perc85
        perc9[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$perc9
        perc95[res[[i]][[j]]$method,res[[i]][[j]]$classifier]=res[[i]][[j]]$perc95
      }
    }
    }
  }
  
  
  return(list(PRresults=PRresults,ROCresults=ROCresults,results85=results85,results95=results95,results9=results9,perc85=perc85,perc95=perc95,perc9=perc9,ROCdf=resROCdf))
}

.myClassificationEvaluationFn=function(testLabelsHat,testLabels,trainingLabels,trainingLabelsHat,testColNames){
  require(PRROC,quietly=T)
  require(pROC,quietly=T)
  
  ROCdf=data.frame(sampleName=testColNames,real=testLabels,estimated=testLabelsHat,stringsAsFactors = F)
  
  finalResPR=pr.curve(scores.class0=testLabelsHat,weights.class0=testLabels)$auc.integral
  finalResROC=as.numeric(roc(response=testLabels,predictor=as.vector(testLabelsHat))$auc)
  
  o=order(trainingLabelsHat,decreasing=T)
  tmpHat=trainingLabelsHat[o]
  tmp=trainingLabels[o]
  tmp=sapply(1:length(tmp),function(x) {sum(tmp[1:x]==1)/x})
  
  thr0.9=1
  if(length(which(tmp>=0.9))>0){
    thr0.9=max(which(tmp>=0.9))
    thr0.9=tmpHat[thr0.9]
  }
  
  thr0.85=1
  if(length(which(tmp>=0.85))>0){
    thr0.85=max(which(tmp>=0.85))
    thr0.85=tmpHat[thr0.85]
  }
  
  thr0.95=1
  if(length(which(tmp>=0.95))>0){
    thr0.95=max(which(tmp>=0.95))
    thr0.95=tmpHat[thr0.95]
  }
  
  thr0.85=max(thr0.85,0.5)
  thr0.9=max(thr0.9,0.5)
  thr0.95=max(thr0.95,0.5)
  
  length0.95=0
  length0.85=0
  length0.9=0
  
  pr0.85=which(testLabelsHat>=thr0.85)
  pr0.85=testLabels[pr0.85]
  length0.85=sum(pr0.85==1)/sum(testLabels==1)
  pr0.85=sum(pr0.85==1)/length(pr0.85)
  
  pr0.9=which(testLabelsHat>=thr0.9)
  pr0.9=testLabels[pr0.9]
  length0.9=sum(pr0.9==1)/sum(testLabels==1)
  pr0.9=sum(pr0.9==1)/length(pr0.9)
  
  pr0.95=which(testLabelsHat>=thr0.95)
  pr0.95=testLabels[pr0.95]
  length0.95=sum(pr0.95==1)/sum(testLabels==1)
  pr0.95=sum(pr0.95==1)/length(pr0.95)
  
  finalRes=list(aucPR=finalResPR,aucROC=finalResROC,pr0.85=pr0.85,pr0.9=pr0.9,pr0.95=pr0.95,perc0.85=length0.85,perc0.9=length0.9,perc0.95=length0.95,ROCdf=ROCdf)
  #print("area under the Precision-recall and ROC curves are being reported")
  return(finalRes)
}

.myClassificationBenchmarkDetailed=function(x,inputDataList,methodsList,expClassName){
  result=list()
  for(i in x){
    
    trainFeatures=inputDataList[[i]][["features"]]
    
    if(!is.null(dim(trainFeatures))&&dim(trainFeatures)[2]>0&&(ncol(trainFeatures)/nrow(trainFeatures))<10){
      trainLabels=inputDataList[[i]][["labels"]]
      tstFeatures=inputDataList[[i]][["testInputData"]]
      tstLabels=inputDataList[[i]][["testLabels"]]
      
      ##Correcting the column names
      colnames(trainFeatures)=gsub(" ", "", colnames(trainFeatures))
      colnames(tstFeatures)=gsub(" ", "", colnames(tstFeatures))
      colnames(trainFeatures)=gsub(":", "", colnames(trainFeatures))
      colnames(tstFeatures)=gsub(":", "", colnames(tstFeatures))
      
      indx=(regexpr("\\d", colnames(trainFeatures))==1)
      colnames(trainFeatures)[indx]=paste0("X",colnames(trainFeatures)[indx])
      indx=(regexpr("\\d", colnames(tstFeatures))==1)
      colnames(tstFeatures)[indx]=paste0("X",colnames(tstFeatures)[indx])
      
      
      for(im in 1:length(methodsList)){
        tmpRes=list()
        try({
        tmpRes=.myClassificationBenchmartSwitch(trainFeatures,trainLabels,tstFeatures,tstLabels,expClassName,methodsList[[im]])
        },silent = T)
        if(length(tmpRes)>0){
          
          #print(paste(i,":",inputDataList[[i]]$method))
          
          if(length(result)==0){
            result=list(list(method=inputDataList[[i]]$method,classifier=methodsList[[im]],aucPR=tmpRes$aucPR,aucROC=tmpRes$aucROC,res9=tmpRes$pr0.9,res95=tmpRes$pr0.95,res85=tmpRes$pr0.85,perc9=tmpRes$perc0.9,perc95=tmpRes$perc0.95,perc85=tmpRes$perc0.85,ROCdf=tmpRes$ROCdf))
          }
          else{
            result=c(result,list(list(method=inputDataList[[i]]$method,classifier=methodsList[[im]],aucPR=tmpRes$aucPR,aucROC=tmpRes$aucROC,res9=tmpRes$pr0.9,res95=tmpRes$pr0.95,res85=tmpRes$pr0.85,perc9=tmpRes$perc0.9,perc95=tmpRes$perc0.95,perc85=tmpRes$perc0.85,ROCdf=tmpRes$ROCdf)))
          }
        }
      }
    } else {
      if(!(!is.null(dim(trainFeatures))&&dim(trainFeatures)[2]>0)){
        print(paste0("error in ",inputDataList[[i]]$method))
      }
      
    }
  }
  
  return(result)
}

.myClassificationBenchmartSwitch=function(trainingInputData,trainingLabels,testInputData,testLabels,eClassName,method){
  switch(method,
         reg=.NeuralTreeBasedWrapper("reg",trainingInputData,trainingLabels,testInputData,testLabels,eClassName=eClassName),
         logReg=.NeuralTreeBasedWrapper("logReg",trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         lda=.myLDA(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         qda=.myQDA(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         sqda=.mySQDA(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         ridgeReg=.myRidgeReg(trainingInputData,trainingLabels,testInputData,testLabels,classification=T,experimentClassName=eClassName),
         lassoReg=.myLassoReg(trainingInputData,trainingLabels,testInputData,testLabels,classification=T,experimentClassName=eClassName),
         ridgeLogReg=.myLogRegRidge(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         lassoLogReg=.myLogRegLasso(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         elasticNetLogReg=.myLogRegElasticNet(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         randomForest=.NeuralTreeBasedWrapper('randomForest',trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         boosting=.NeuralTreeBasedWrapper('boosting',trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         bagging=.NeuralTreeBasedWrapper('bagging',trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
         neuralNet=.NeuralTreeBasedWrapper('neuralNet',trainingInputData,trainingLabels,testInputData,testLabels,eClassName))
}

.NeuralTreeBasedWrapper=function(classificationMethod,trainingInputData,trainingLabels,testInputData,testLabels,eClassName){
  if((ncol(trainingInputData)/nrow(trainingInputData))<2){
    switch(classificationMethod,
           randomForest=.myRandomForest(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
           boosting=.myBoosing(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
           bagging=.myBagging(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
           neuralNet=.myNeuralNet(trainingInputData,trainingLabels,testInputData,testLabels,eClassName),
           reg=.myReg(trainingInputData,trainingLabels,testInputData,testLabels,classification=T,experimentClassName=eClassName),
           logReg=.myLogReg(trainingInputData,trainingLabels,testInputData,testLabels,eClassName))
  } else {
    return(list())
  }
}
.myRegOriginal=function(trainingInputData,trainingLabels,testInputData,testLabels,classification=F,experimentClassName=NULL){
  require(MASS,quietly=T)
  require(car,quietly=T)
  
  cotrol="A"
  if(classification)
    cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  if(classification){
    if(is.null(experimentClassName))
      stop("Experiment class name should be provided")
    if(length(levels(as.factor(trainingLabels)))>2){
      print("error: more than two levels exit in the training labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
      trainingLabels=trainingLabels-1
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    }
  } else{
    if(is.numeric(trainingLabels)){
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    } else {
      print("input training labels are supposed to be numeric")
      stop()
    }
  }
  
  if(classification){
    if(length(levels(as.factor(testLabels)))>2){
      print("error: more than two levels exit in the test labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
      testLabels=testLabels-1
    }
  } else{
    if(!is.numeric(testLabels)){
      print("input test labels are supposed to be numeric")
      stop()
    }
  } 
  
  colnames(dataTrainingPrepared)[1]="response"
  formulaStr=paste0("response~",paste(colnames(dataTrainingPrepared)[-1],collapse = "+"))
  
  resLm=lm(formula = formula(formulaStr),as.data.frame(dataTrainingPrepared))
  #print(summary(resLm))
  vifRes=car::vif(resLm)
  if(max(vifRes)>10){
    print("collinearity was detected among covariates!")
    print("variance inflation factors:")
    print(car::vif(resLm))
  } else{
    print("no strong collinearity was detected")
  }
  
  rstudentRes=rstudent(resLm)
  if(max(abs(rstudentRes))<3){
    print("no outlier was detected")
  } else {
    print(paste("outliers were detected:",length(which(abs(rstudentRes)>3))))
    print(paste(row.names(trainingInputData)[which(abs(rstudentRes)>3)],collapse=", "))
  }
  
  leverageRes=hatvalues(resLm)
  if(max(abs(leverageRes))<3){
    print("no high leverage point was detected")
  } else {
    print(paste("high leverage points were detected:",length(which(abs(leverageRes)>3))))
    print(paste(row.names(data)[which(abs(leverageRes)>3)],collapse=", "))
  }
  
  testLabelsHat=predict(resLm,as.data.frame(testInputData))
  if(classification){
    testLabelsHat[testLabelsHat<0]=0
    testLabelsHat[testLabelsHat>1]=1
    finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,testColNames=row.names(testInputData))
  } else {
    finalRes <- mean((testLabels - testLabelsHat)^2)
    print("RMSE values are being reported")
  }
  
  return(finalRes)
}

.myReg=function(trainingInputData,trainingLabels,testInputData,testLabels,classification=F,experimentClassName=NULL){
  require(MASS,quietly=T)
  require(car,quietly=T)
  
  cotrol="A"
  if(classification)
    cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  if(classification){
    if(is.null(experimentClassName))
      stop("Experiment class name should be provided")
    if(length(levels(as.factor(trainingLabels)))>2){
      print("error: more than two levels exit in the training labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
      trainingLabels=trainingLabels-1
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    }
  } else{
    if(is.numeric(trainingLabels)){
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    } else {
      print("input training labels are supposed to be numeric")
      stop()
    }
  }
  
  if(classification){
    if(length(levels(as.factor(testLabels)))>2){
      print("error: more than two levels exit in the test labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
      testLabels=testLabels-1
    }
  } else{
    if(!is.numeric(testLabels)){
      print("input test labels are supposed to be numeric")
      stop()
    }
  } 
  
  colnames(dataTrainingPrepared)[1]="response"
  formulaStr=paste0("response~",paste(colnames(dataTrainingPrepared)[-1],collapse = "+"))
  
  resLm=lm(formula = formula(formulaStr),as.data.frame(dataTrainingPrepared))
  #print(summary(resLm))
  #vifRes=car::vif(resLm)
  #if(max(vifRes)>10){
  #  print("collinearity was detected among covariates!")
  #  print("variance inflation factors:")
  #  print(car::vif(resLm))
  #} else{
  #  print("no strong collinearity was detected")
  #}
  
  
  
  testLabelsHat=predict(resLm,as.data.frame(testInputData))
  trainingLabelsHat=predict(resLm,as.data.frame(dataTrainingPrepared))
  if(classification){
    testLabelsHat[testLabelsHat<0]=0
    testLabelsHat[testLabelsHat>1]=1
    testLabelsHat[testLabelsHat<0]=0
    testLabelsHat[testLabelsHat>1]=1
    finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  } else {
    finalRes <- mean((testLabels - testLabelsHat)^2)
    print("RMSE values are being reported")
  }
  
  return(finalRes)
}

.myLogReg=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
   
  
  colnames(dataTrainingPrepared)[1]="response"
  formulaStr=paste0("response~",paste(colnames(dataTrainingPrepared)[-1],collapse = "+"))
  
  
  
  resGlm=glm(formula = formula(formulaStr),data=as.data.frame(dataTrainingPrepared),family = binomial)
  #print(summary(resGlm))
  
  testLabelsHat=predict(resGlm,as.data.frame(testInputData),type="response")
  
  trainingLabelsHat=predict(resGlm,as.data.frame(dataTrainingPrepared),type="response")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
 
  return(finalRes)
}

.myLDA=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  require(MASS,quietly=T)
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  
  
  colnames(dataTrainingPrepared)[1]="response"
  formulaStr=paste0("response~",paste(colnames(dataTrainingPrepared)[-1],collapse = "+"))
  
  finalRes=list()
  try({
    
    
  resLda=lda(formula = formula(formulaStr),data=as.data.frame(dataTrainingPrepared))
  #print(resLda)
  
  testLabelsHat=predict(resLda,as.data.frame(testInputData))
  testLabelsHat=testLabelsHat$posterior[,2]
  
  trainingLabelsHat=predict(resLda,as.data.frame(dataTrainingPrepared))
  trainingLabelsHat=trainingLabelsHat$posterior[,2]
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  },silent=T)
  return(finalRes)
}

.myQDA=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  require(MASS,quietly=T)
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  
  
  colnames(dataTrainingPrepared)[1]="response"
  formulaStr=paste0("response~",paste(colnames(dataTrainingPrepared)[-1],collapse = "+"))
  
  finalRes=list()
  try({
    
    
  resQda=qda(formula = formula(formulaStr),data=as.data.frame(dataTrainingPrepared))
  #print(resQda)
  
  testLabelsHat=predict(resQda,as.data.frame(testInputData))
  testLabelsHat=testLabelsHat$posterior[,2]
  
  
  trainingLabelsHat=predict(resQda,as.data.frame(dataTrainingPrepared))
  trainingLabelsHat=trainingLabelsHat$posterior[,2]
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  },silent = T)
  return(finalRes)
}

.mySQDA=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  require(SQDA,quietly=T)
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=t(.matrixExtraction(trainingInputData))
  testInputData=t(.matrixExtraction(testInputData))
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
  }
  
  colnames(trainingInputData)=trainingLabels
  rownames(trainingInputData)=paste0("X",rownames(trainingInputData))
  colnames(testInputData)=testLabels
  rownames(testInputData)=paste0("X",rownames(testInputData))
  
  res=.mysQDA(trainingInputData,testInputData,lams=0.2,presel=FALSE)
  testLabelsHat=res$pred
  testLabelsHat=as.numeric(factor(testLabelsHat,levels=c(cotrol,experimentClassName)))-1
  testLabels=as.numeric(factor(testLabels,levels=c(cotrol,experimentClassName)))-1
  
  ##needs correction
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.myRidgeReg=function(trainingInputData,trainingLabels,testInputData,testLabels,classification=F,experimentClassName=NULL){
  require(MASS,quietly=T)
  require(glmnet,quietly=T)
  
  
  cotrol="A"
  if(classification)
    cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  if(classification){
    if(is.null(experimentClassName))
      stop("Experiment class name should be provided")
    if(length(levels(as.factor(trainingLabels)))>2){
      print("error: more than two levels exit in the training labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
      trainingLabels=trainingLabels-1
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    }
  } else{
    if(is.numeric(labels)){
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    } else {
      print("input training labels are supposed to be numeric")
      stop()
    }
  }
  
  if(classification){
    if(length(levels(as.factor(testLabels)))>2){
      print("error: more than two levels exit in the test labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
      testLabels=testLabels-1
    }
  } else{
    if(!is.numeric(labels)){
      print("input test labels are supposed to be numeric")
      stop()
    }
  } 
  
  grid=10^seq(10,-6,length=400)
  cv.out=cv.glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=0,lambda=grid,type.measure="deviance")
  
  
  
  ridge.mod=glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=0,lambda=cv.out$lambda.min,thresh=1e-12)
  testLabelsHat=predict(ridge.mod,s=cv.out$lambda.min,newx=testInputData,type = "response")
  trainingLabelsHat=predict(ridge.mod,s=cv.out$lambda.min,newx=dataTrainingPrepared[,-1],type = "response")
  if(classification){
    testLabelsHat[testLabelsHat<0]=0
    testLabelsHat[testLabelsHat>1]=1
    trainingLabelsHat[trainingLabelsHat<0]=0
    trainingLabelsHat[trainingLabelsHat>1]=1
    finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
    } else {
    finalRes <- mean((testLabels - testLabelsHat)^2)
    print("RMSE values are being reported")
  }
  
  return(finalRes)
}

.myLassoReg=function(trainingInputData,trainingLabels,testInputData,testLabels,classification=F,experimentClassName=NULL){
  require(MASS,quietly=T)
  require(glmnet,quietly=T)
  
  cotrol="A"
  if(classification)
    cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  if(classification){
    if(is.null(experimentClassName))
      stop("Experiment class name should be provided")
    if(length(levels(as.factor(trainingLabels)))>2){
      print("error: more than two levels exit in the training labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
      trainingLabels=trainingLabels-1
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    }
  } else{
    if(is.numeric(labels)){
      dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
    } else {
      print("input training labels are supposed to be numeric")
      stop()
    }
  }
  
  if(classification){
    if(length(levels(as.factor(testLabels)))>2){
      print("error: more than two levels exit in the test labels")
      stop()
    } else {
      if(length(cotrol)>1)
        stop("Provided experiment class name couldn't be found")
      testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
      testLabels=testLabels-1
    }
  } else{
    if(!is.numeric(labels)){
      print("input test labels are supposed to be numeric")
      stop()
    }
  } 
  
  grid=10^seq(10,-6,length=400)
  cv.out=cv.glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=1,lambda=grid)
  
  
  lasso.mod=glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=1,lambda=cv.out$lambda.min,thresh=1e-12)
  testLabelsHat=predict(lasso.mod,s=cv.out$lambda.min,newx=testInputData,type = "response")
  trainingLabelsHat=predict(lasso.mod,s=cv.out$lambda.min,newx=dataTrainingPrepared[,-1],type = "response")
  if(classification){
    testLabelsHat[testLabelsHat<0]=0
    testLabelsHat[testLabelsHat>1]=1
    
    trainingLabelsHat[trainingLabelsHat<0]=0
    trainingLabelsHat[trainingLabelsHat>1]=1
    
    finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
    
    } else {
    finalRes <- mean((testLabels - testLabelsHat)^2)
    print("RMSE values are being reported")
  }
  
  return(finalRes)
}

.myLogRegLasso=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(glmnet,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  
  
  grid=10^seq(10,-6,length=400)
  
  cv.out=cv.glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=1,lambda=grid,family="binomial")
  
  lasso.mod=glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=1,lambda=cv.out$lambda.min,thresh=1e-12,family="binomial")
  testLabelsHat=predict(lasso.mod,s=cv.out$lambda.min,newx=testInputData,type = "response")
  trainingLabelsHat=predict(lasso.mod,s=cv.out$lambda.min,newx=dataTrainingPrepared[,-1],type = "response")
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  finalRes$interanlROC=(-1)
  return(finalRes)
}

.myLogRegRidge=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(glmnet,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  
  
  grid=10^seq(10,-6,length=400)
  
  cv.out=cv.glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=0,lambda=grid,family="binomial")
  
  
  
  ridge.mod=glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=0,lambda=cv.out$lambda.min,thresh=1e-12,family="binomial")
  testLabelsHat=predict(ridge.mod,s=cv.out$lambda.min,newx=testInputData,type = "response")
  trainingLabelsHat=predict(ridge.mod,s=cv.out$lambda.min,newx=dataTrainingPrepared[,-1],type = "response")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.myLogRegElasticNet=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  require(glmnet,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  
  
  grid=10^seq(10,-6,length=400)
  alphasOfInterest<-seq(0,1,by=0.05) 
  
  cvs<-lapply(alphasOfInterest, function(curAlpha){
    cv.glmnet(trainingInputData, trainingLabels, alpha=curAlpha,lambda=grid,family="binomial")
  })
  
  indx=which.min(unlist(lapply(cvs,function(x){mean(x$cvm)})))
  alphaSel=alphasOfInterest[indx]
  lambdaSel=cvs[[indx]]$lambda.min
  
  elasticNet.mod=glmnet(dataTrainingPrepared[,-1],dataTrainingPrepared[,1],alpha=alphaSel,lambda=lambdaSel,thresh=1e-12,family="binomial")
  
  testLabelsHat=predict(elasticNet.mod,s=lambdaSel,newx=testInputData,type = "response")
  trainingLabelsHat=predict(elasticNet.mod,s=lambdaSel,newx=dataTrainingPrepared[,-1],type = "response")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
  
}

.myRandomForest=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(randomForest,quietly=T)
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  colnames(dataTrainingPrepared)[1]="response"
  colnames(dataTrainingPrepared)=paste0("X",colnames(dataTrainingPrepared))
  colnames(testInputData)=paste0("X",colnames(testInputData))
  
  forest.fit=randomForest(as.factor(Xresponse)~.,data=dataTrainingPrepared,mtry=round(sqrt(dim(trainingInputData)[2])),ntree=1000,importance=T,type=classification, norm.votes=T)
  testLabelsHat=predict(forest.fit,newdata=testInputData, type = "prob")
  trainingLabelsHat=predict(forest.fit,newdata=dataTrainingPrepared, type = "prob")
  finalRes=.myClassificationEvaluationFn(testLabelsHat[,2],testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.myBoosing=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(gbm,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  colnames(dataTrainingPrepared)[1]="response"
  
  gbm.fit=gbm(as.character(response)~.,data=as.data.frame(dataTrainingPrepared),distribution = "bernoulli",interaction.depth =max(round(sqrt(dim(trainingInputData)[2]-1))-2,3),n.trees=5000)
  testLabelsHat=predict(gbm.fit,newdata=as.data.frame(testInputData), n.trees=5000,type="response")
  trainingLabelsHat=predict(gbm.fit,newdata=as.data.frame(dataTrainingPrepared), n.trees=5000,type="response")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.myBagging=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  
  require(randomForest,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  colnames(dataTrainingPrepared)[1]="response"
  
  forest.fit=randomForest(as.factor(response)~.,data=dataTrainingPrepared,mtry=dim(trainingInputData)[2],ntree=1000,importance=T,type=classification, norm.votes=T)
  testLabelsHat=predict(forest.fit,newdata=testInputData, type = "prob")
  trainingLabelsHat=predict(forest.fit,newdata=dataTrainingPrepared, type = "prob")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat[,2],testLabels,trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.mySVMLinear=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(e1071,quietly=T)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    trainingLabels=as.numeric(factor(as.character(trainingLabels),levels=c(cotrol,experimentClassName)))
    trainingLabels=trainingLabels-1
    dataTrainingPrepared=cbind(trainingLabels,trainingInputData)
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    if(length(cotrol)>1)
      stop("Provided experiment class name couldn't be found")
    testLabels=as.numeric(factor(as.character(testLabels),levels=c(cotrol,experimentClassName)))
    testLabels=testLabels-1
  }
  colnames(dataTrainingPrepared)[1]="response"
  
  forest.fit=randomForest(as.factor(response)~.,data=dataTrainingPrepared,mtry=dim(trainingInputData)[2],ntree=1000,importance=T,type=classification, norm.votes=T)
  testLabelsHat=predict(forest.fit,newdata=testInputData, type = "prob")
  
  finalRes=.myClassificationEvaluationFn(testLabelsHat,testLabels,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.myNeuralNet=function(trainingInputData,trainingLabels,testInputData,testLabels,experimentClassName){
  require(nnet,quietly=T)
  
  cotrol=setdiff(unique(as.character(trainingLabels)),experimentClassName)
  
  if(dim(trainingInputData)[1]!=length(trainingLabels)||dim(testInputData)[1]!=length(testLabels))
    stop("wrong dimensions")
  
  if(dim(trainingInputData)[2]!=dim(testInputData)[2])
    stop("dimensions of training and test input data do not match!")
  
  for(i in 1:dim(trainingInputData)[2]){
    if(is.character(trainingInputData[,i]))
      trainingInputData[,i]=as.factor(trainingInputData[,i])
    
    if(is.character(testInputData[,i]))
      testInputData[,i]=as.factor(testInputData[,i])
  }
  
  for(i in 1:dim(testInputData)[2]){
    if(!all(levels(trainingInputData[,i])==levels(testInputData[,i])))
      stop(paste("levels of training and test data do not match for",colnames(trainingInputData)[i]))
  }
  
  
  trainingInputData=.matrixExtraction(trainingInputData)
  testInputData=.matrixExtraction(testInputData)
  
  
  if(length(levels(as.factor(trainingLabels)))>2){
    print("error: more than two levels exit in the training labels")
    stop()
  } else {
    trainingLabels=factor(trainingLabels,levels=c(cotrol,experimentClassName))
    
    dataTrainingPrepared=cbind(trainingLabels,data.frame(trainingInputData))
  }
  
  if(length(levels(as.factor(testLabels)))>2){
    print("error: more than two levels exit in the test labels")
    stop()
  } else {
    testLabels=factor(as.character(testLabels),levels=c(cotrol,experimentClassName))
  }
  colnames(dataTrainingPrepared)[1]="response"
  
  sizeGrid=seq(2,51,3)
  decayGrid=10^seq(10,-6,length=50)
  
  dy=data.frame(var=dataTrainingPrepared[,1])
  model=model.matrix(~0+var,data=dy)
  res=list()
  for(i in 1:length(sizeGrid)){
    cvs<-list(lapply(decayGrid, function(curDecay){
      folds=sample(1:10,dim(trainingInputData)[1],replace=T)
      meanRes=0
      for(j in 1:10){
        
        capture.output(ir1 <- nnet(dataTrainingPrepared[folds!=j,-1], model[folds!=j,], size = sizeGrid[i], rang = 0.1,
                    decay = curDecay, maxit = 1000,softmax=TRUE))
        yhat=predict(ir1,dataTrainingPrepared[folds==j,-1])
        meanRes=meanRes+pr.curve(scores.class0=yhat[,2],weights.class0=(as.numeric(dataTrainingPrepared[folds==j,1])-1))$auc.integral
      }
      return(meanRes/10)
    }))
    res=c(res,cvs)
  }
  
  indx=unlist(lapply(res,function(x){which.max(unlist(x))}))
  selInd=1
  maxScore=0
  for(i in 1:length(indx)){
    if(res[[i]][[indx[i]]]>maxScore){
      maxScore=res[[i]][[indx[i]]]
      selInd=i
    }
  }
  
  capture.output(irF <- nnet(dataTrainingPrepared[,-1], model[,], size = sizeGrid[selInd], rang = 0.1,
                             decay = decayGrid[indx[selInd]], maxit = 1000,softmax=TRUE))
  
  yhat=predict(irF,testInputData)
  trainingLabelsHat=predict(irF,dataTrainingPrepared)
  finalRes=.myClassificationEvaluationFn(yhat[,2],(as.numeric(testLabels)-1),trainingLabels=dataTrainingPrepared[,1],trainingLabelsHat=trainingLabelsHat,testColNames=row.names(testInputData))
  
  return(finalRes)
}

.mysQDA=function (train.data = NULL, test.data = NULL, len = 100, lams = seq(0.02, 1, length = 10), presel = T, prelam = 0.2, margin = 0.05) 
{
  require(SQDA,quietly = T)
  group <- unique(colnames(train.data))
  table <- SQDA:::sortgene(train.data)
  seq <- rownames(table)
  dat.sel = train.data[seq, ]
  dat.sel.ts = test.data[seq, ]
  times = as.integer(nrow(train.data)/len)
  k <- (1:times)
  if (presel) {
    res <- SQDA:::simpleAGG3(data = dat.sel, data.new = dat.sel.ts, 
                             len = len, times = times, lam = prelam)
    k <- (1:times)[res$cv.error <= min(res$cv.error) + margin]
  }
  lik <- NULL
  cverror = NULL
  for (j in k) {
    seq = (1:len) + len * (j - 1)
    dat.sel2 <- dat.sel[seq, ]
    dat.sel.ts2 = dat.sel.ts[seq, ]
    res <- SQDA:::simpleAGG3(data = dat.sel2, data.new = dat.sel.ts2, 
                             len = len, times = 1, lam = lams)
    lik <- rbind(lik, res$lik)
    cverror <- c(cverror, res$cv.error)
  }
  j <- (1:length(k))
  s <- sort(c(2 * j, 2 * j - 1))
  lik <- lik[s, ]
  p = 1
  if (nrow(lik) > 2) {
    p = nrow(lik)/2
    seq1 = seq(1, nrow(lik), by = 2)
    seq2 = seq(2, nrow(lik), by = 2)
    res <- rbind(apply(lik[seq1, ], 2, sum), apply(lik[seq2, 
                                                       ], 2, sum))
  }
  else {
    res <- lik
  }
  rownames(res) <- group
  colnames(res) <- colnames(test.data)
  pred <- group[max.col(t(res))]
  return(list(pred = pred, p = p))
}

.matrixExtraction=function (object) 
{
  y <- list()
  if (is(object, "list")) {
    if (is(object, "EList")) {
      y$exprs <- as.matrix(object$E)
      y$Amean <- rowMeans(y$exprs, na.rm = TRUE)
    }
    else {
      if (is(object, "EListRaw")) 
        stop("EListRaw object: please normalize first")
      if (is(object, "RGList")) 
        stop("RGList object: please normalize first")
      y$printer <- object$printer
      if (is.null(object$M)) 
        stop("data object isn't of a recognized data class")
      y$exprs <- as.matrix(object$M)
      if (!is.null(object$A)) 
        y$Amean <- rowMeans(as.matrix(object$A), na.rm = TRUE)
    }
    y$weights <- object$weights
    y$probes <- object$genes
    y$design <- object$design
  }
  else {
    if (is(object, "ExpressionSet")) {
      if (!requireNamespace("Biobase", quietly = TRUE)) 
        stop("Biobase package required but is not available")
      y$exprs <- Biobase::exprs(object)
    }
    else {
      if (is(object, "PLMset")) {
        y$exprs <- object@chip.coefs
        if (length(y$exprs) == 0) 
          stop("chip.coefs has length zero")
      }
      else {
        if (is(object, "marrayNorm")) {
          y$exprs <- object@maM
        }
        else {
          if (is(object, "eSet")) {
            if (!requireNamespace("Biobase", quietly = TRUE)) 
              stop("Biobase package required but is not available")
            y$exprs <- Biobase::assayDataElement(object, 
                                                 "exprs")
            if (length(object@featureData@data)) 
              y$probes <- object@featureData@data
            y$Amean <- rowMeans(y$exprs, na.rm = TRUE)
            if ("weights" %in% Biobase::assayDataElementNames(object)) 
              y$weights <- Biobase::assayDataElement(object, 
                                                     "weights")
          }
          else {
            if (is.vector(object)) 
              y$exprs <- matrix(object, nrow = 1)
            else y$exprs <- data.matrix(object)
          }
        }
      }
    }
  }
  return(y$exprs)
}