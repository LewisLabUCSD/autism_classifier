#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#setwd("~/classificationModuleJabba/classificationCode/3- mainCode")
source("WGCNAandGeneFilterationMethods.R")
source("pipelines.R")
source("ClassificationModule.R")

#rm(list=ls())
#args="~/Desktop/"
#args="classificationSet_subject_scores/out1/"
#args="~/OneDrive - UC San Diego/Vahid/OneDrive - UC San Diego/classificationPaper/newData/server_results/classificationDiscoverySet/out1/"
.myfoldRunnerFn=function(){
  require(Biobase,quietly = T)
  #arg1: directory address for each fold
  
  path=gsub("\"", "", args[1])
  itrCounter=0
  while((!file.exists(paste0(path,"results.rda")))|(itrCounter<1)){
    itrCounter=itrCounter+1
    load(paste0(path,'data.rda'))
    #middleFns=list("pipeline6")
    #inputData=list(data)
    #methodList=middleFns
    #ncores=1
    resultsMiddle=.myMiddleFn(list(data),methodList=middleFns,ncores=1,expClassName=expClassName)
    #posteriorFns=list('wgcna')
    resultsPosterior=.myPosteriorFn(resultsMiddle,methodList=posteriorFns,ncores=1,expClassName=expClassName)
    
    #classificationMethodsList2=classificationMethodsList
    #classificationMethodsList=list('bagging')
    result=.myClassificationBenchmarking(resultsPosterior,ncores=1,expClassName=expClassName,methodsList=classificationMethodsList)
    save(result,file=paste0(path,'results.rda'))
  }
}

.myFoldMakerFn=function(inputLabels,sampleNames,nfold){
  inputLabels=as.factor(inputLabels)
  levLabels=levels(inputLabels)
  lev1=which(labels==levLabels[1])
  lev2=which(labels==levLabels[2])
  nlev1=rep(1,ceiling(length(lev1)/nfold))
  for(i in 2:nfold)
    nlev1=c(nlev1,rep(i,ceiling(length(lev1)/nfold)))
  folds1=sample(nlev1,length(lev1),replace = F)
  
  nlev2=rep(1,ceiling(length(lev2)/nfold))
  for(i in 2:nfold)
    nlev2=c(nlev2,rep(i,ceiling(length(lev2)/nfold)))
  folds2=sample(nlev2,length(lev2),replace = F)
  
  df=data.frame(sampleName=sampleNames[lev1],label=labels[lev1],fold=folds1,stringsAsFactors = F)
  df=rbind(df,data.frame(sampleName=sampleNames[lev2],label=labels[lev2],fold=folds2,stringsAsFactors = F))
  return(df)
}

.myBenchmarkingWrapperDetailed=function(x,inputExpSet,tstLabelsList,expClassName,ncores,initializerFns,middleFns,posteriorFns,classificationMethodsList){
  
  result=list()
  
  for(i in x){
    trainingInput=inputExpSet[,!tstLabelsList[[i]]$testFlag]
    trainingLabels=tstLabelsList[[i]]$inputLabels[!tstLabelsList[[i]]$testFlag]
    testInput=inputExpSet[,tstLabelsList[[i]]$testFlag]
    testLabels=tstLabelsList[[i]]$inputLabels[tstLabelsList[[i]]$testFlag]
    
    inputDataList=.myFeatureSelectionBenchmarking(data=trainingInput,labels=trainingLabels,testInputData=testInput,testsLabels=testLabels,ncores=ncores,initializerFns=initializerFns,middleFns=middleFns,posteriorFns=posteriorFns)
    tmpRes=.myClassificationBenchmarking(inputDataList,ncores=ncores,expClassName=expClassName,methodsList=classificationMethodsList)
    if(length(result)==0){
      result=list(list(PRresult=tmpRes$PRresults,ROCresults=tmpRes$ROCresults,cv=x$status))
    } else {
      result=c(result,list(list(PRresult=tmpRes$PRresults,ROCresults=tmpRes$ROCresults,cv=x$status)))
    }
  }
  return(list(PRresult=result$PRresults,ROCresults=result$ROCresults,cv=x$status))
}

.myfoldRunnerFn()
