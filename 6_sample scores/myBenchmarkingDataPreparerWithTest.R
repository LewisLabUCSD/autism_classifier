#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args=c("inputExpData_ageBalanced_TD39.rda","classificationSet_subject_scores/")
source('pipelines.R')
.myBenchmarkingDataPreparer=function(){
  require(Biobase,quietly = T)
  
  load(args[1])
  output.directory=args[2]
  if(!file.exists(output.directory)){
    dir.create(output.directory)
  } else {
    unlink(output.directory, recursive = T, force = T)
    dir.create(output.directory)
  }
  
  
  outputfilesList=""
  
  ncores=1
  counter=1
  outputfilesList=""
    
  trainingInput=inputExpData
  trainingLabels=labels
  testInput=testInputExpData
  
  results=.myInitializer(trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,method=initializerFns,prevMethod="")
  for(j in 1:length(results)){
    dir.create(paste0(output.directory,"out",counter))
    data=results[[j]]
    save(data,middleFns,posteriorFns,classificationMethodsList,expClassName,file=paste0(output.directory,"out",counter,"/data.rda"))
    outputfilesList=c(outputfilesList,paste0(output.directory,"out",counter,"/"))
      
    counter=counter+1
  }
  
  
  outputfilesList=outputfilesList[-1]
  write.table(outputfilesList,file=paste0(output.directory,"dataList.txt"),row.names = F,col.names = F)
}

.myFoldMakerFn=function(inputLabels,sampleNames,nfold){
  inputLabels=as.factor(inputLabels)
  levLabels=levels(inputLabels)
  lev1=which(inputLabels==levLabels[1])
  lev2=which(inputLabels==levLabels[2])
  nlev1=rep(1,ceiling(length(lev1)/nfold))
  for(i in 2:nfold)
    nlev1=c(nlev1,rep(i,ceiling(length(lev1)/nfold)))
  folds1=sample(nlev1,length(lev1),replace = F)
  
  nlev2=rep(1,ceiling(length(lev2)/nfold))
  for(i in 2:nfold)
    nlev2=c(nlev2,rep(i,ceiling(length(lev2)/nfold)))
  folds2=sample(nlev2,length(lev2),replace = F)
  
  df=data.frame(sampleName=sampleNames[lev1],label=inputLabels[lev1],fold=folds1,stringsAsFactors = F)
  df=rbind(df,data.frame(sampleName=sampleNames[lev2],label=inputLabels[lev2],fold=folds2,stringsAsFactors = F))
  return(df)
}

.myBenchmarkingDataPreparer()
