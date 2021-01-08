source("WGCNAandGeneFilterationMethods.R")
source("pipelines.R")
source("ClassificationModule.R")
.myBenchmarkingWrapper=function(inputExpData,
                                labels,
                                expClassName="proband",
                                ncores=1,
                                nfold=5,
                                npermTest=1,
                                initializerFns=list("no","cov","var","cov_var","varImportance"),
                                middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline4","pipeline5"),
                                posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls'),
                                classificationMethodsList=list("reg","logReg","lda","qda","ridgeReg","lassoReg","ridgeLogReg","lassoLogReg","elasticNetLogReg","randomForest","boosting","bagging","neuralNet")){
  require(parallel)

  
  folds=.myFoldMakerFn(labels,colnames(inputExpData),nfold)
  testlabelsList=list(list(status="real",testFlag=folds$fold==1,inputLabels=labels))
  for (i in 2:nfold){
    testlabelsList=c(testlabelsList,list(list(status=paste0("real"),testFlag=folds$fold==i,inputLabels=labels)))
  }
  
  if(npermTest>0){
    for(i in 1:npermTest){
      tmpLabels=sample(labels,size=length(labels),replace = F)
      folds=.myFoldMakerFn(tmpLabels,colnames(inputExpData),nfold)
      for (j in 1:nfold){
        testlabelsList=c(testlabelsList,list(list(status=paste0("perm",i),testFlag=folds$fold==i,inputLabels=tmpLabels)))
      }
    }
  }
  
  
  for(i in 1:length(testlabelsList)){
   
    trainingInput=inputExpData[,!testlabelsList[[i]]$testFlag]
    trainingLabels=testlabelsList[[i]]$inputLabels[!testlabelsList[[i]]$testFlag]
    testInput=inputExpData[,testlabelsList[[i]]$testFlag]
    testLabels=testlabelsList[[i]]$inputLabels[testlabelsList[[i]]$testFlag]
    
    inputDataList=.myFeatureSelectionBenchmarking(data=trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,ncores=ncores,initializerFns=initializerFns,middleFns=middleFns,posteriorFns=posteriorFns)
    
    if(length(result)==0){
      result=list(inputDataList)
    } else {
      result=c(result,list(inputDataList))
    }
  }
  
  result=.myClassificationBenchmarking(inputDataList,ncores=ncores,expClassName=expClassName,methodsList=classificationMethodsList)
  return(result)
}

.myfoldRunnerFn=function(){
  require(parallel)
  
  trainingInput=inputExpData[,!testLabels$testFlag]
  trainingLabels=testLabels$inputLabels[!testLabels$testFlag]
  testInput=inputExpData[,testLabels$testFlag]
  testLabels=testLabels$inputLabels[testLabels$testFlag]
    
  inputDataList=.myFeatureSelectionBenchmarking(data=trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,ncores=ncores,initializerFns=initializerFns,middleFns=middleFns,posteriorFns=posteriorFns)
    
  if(length(result)==0){
    result=list(inputDataList)
  } else {
    result=c(result,list(inputDataList))
  }
  
  result=.myClassificationBenchmarking(inputDataList,ncores=ncores,expClassName=expClassName,methodsList=classificationMethodsList)
  return(result)
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


print(".myBenchmarkingWrapper(inputExpData,
      labels,
      expClassName=proband,
      ncores=1,
      nfold=10,
      npermTest=1,  ##for each permutation, nfold cv will be performed. should be non-negative
      posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls'),
      middleFns=list(no,pipeline1,pipeline2,pipeline3,pipeline4,pipeline5),
      initializerFns=list(no,cov,var,cov_var,varImportance),
      classificationMethodsList=list(reg,logReg,lda,qda,ridgeReg,lassoReg,ridgeLogReg,lassoLogReg,elasticNetLogReg,randomForest,boosting,bagging,neuralNet)")



