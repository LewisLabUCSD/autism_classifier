#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#rm(list=ls())
#setwd("~/classificationModuleJabba/classificationCode/3- mainCode/")
#args=c("inputExpData_ageBalanced_TD39.rda","classificationSet_wSbjScores/")
#runReal=T
#args=c("inputExpData_ageBalanced_TD39.rda","classificationSet_subject_scores/")

source('pipelines.R')

.myFoldMakerFn=function(inputLabels,sampleNames,nfold){
  df=data.frame(lbl=inputLabels,sampleName=sampleNames,stringsAsFactors = F)
  df=df[sample(nrow(df)),]
  
  lev1=which(df$lbl==unique(df$lbl)[1])
  lev2=which(df$lbl==unique(df$lbl)[2])
  foldLev1=split(1:length(lev1),cut(1:length(lev1),nfold))
  foldLev2=split(1:length(lev2),cut(1:length(lev2),nfold))
  
  df1=df[df$lbl==unique(df$lbl)[1],]
  df2=df[df$lbl==unique(df$lbl)[2],]
  
  
  df1$fold="no"
  df2$fold="no"
  for(i in 1:nfold){
    df1$fold[foldLev1[[i]]]=i
    df2$fold[foldLev2[[i]]]=i
  }
  
  df=rbind(df1,df2)
  df=df[,c("sampleName","lbl","fold")]
  colnames(df)=c("sampleName","label","fold")
  df[,1]=as.character(df[,1])
  df[,2]=as.character(df[,2])
  df[,3]=as.numeric(df[,3])
  
  
  return(df)
}

.myFoldCheckerFn=function(inputPath){
  require(data.table)
  x=fread(inputPath,header=F)
  
  
  for(i in 1:nrow(x)){
    load(paste0(x[i],"data.rda"))
    tst=summary(as.factor(data$testLabels))
    tst=data.table(name=names(tst),value=tst,fold=i)
    if(i==1){
      res=tst
    } else {
      res=rbind(res,tst)
    }
    rm(data,classificationMethodsList,middleFns,posteriorFns,tst,expClassName)
  }
  
  res=dcast(data = res,formula = fold~name,fun.aggregate = sum,value.var = "value")
  slclmns=setdiff(colnames(res),"fold")
  res$fraction=as.numeric(res[,get(slclmns[1])])/as.numeric(res[,get(slclmns[2])])
  
  print(paste0("assessing the label distributions in ",nrow(x), " files"))
  if(sum(res$fraction>2|res$fraction<0.5)>0){
    print("error in ratio of class labels for test sets")
  } else {
    print("no error was found in test sets")
  }
  
  
  for(i in 1:nrow(x)){
    load(paste0(x[i],"data.rda"))
    tst=summary(as.factor(data$labels))
    tst=data.table(name=names(tst),value=tst,fold=i)
    if(i==1){
      res=tst
    } else {
      res=rbind(res,tst)
    }
    rm(data,classificationMethodsList,middleFns,posteriorFns,tst,expClassName)
  }
  
  res=dcast(data = res,formula = fold~name,fun.aggregate = sum,value.var = "value")
  slclmns=setdiff(colnames(res),"fold")
  res$fraction=as.numeric(res[,get(slclmns[1])])/as.numeric(res[,get(slclmns[2])])
  
  if(sum(res$fraction>2|res$fraction<0.5)>0){
    print("error in ratio of class labels for training sets")
  } else {
    print("no error was found in training sets")
  }
}


.myBenchmarkingDataPreparer=function(){
  #arg1: inputData
  #arg2: output directory
  require(Biobase,quietly = T)
  
  #.myPackageInstaller()
  
  .myPackageChecker()
  
  set.seed(12345)
  load(args[1])
  output.directory=args[2]
  if(!file.exists(output.directory)){
    dir.create(output.directory)
  } else {
    unlink(output.directory, recursive = T, force = T)
    dir.create(output.directory)
  }
  
  
  folds=.myFoldMakerFn(inputLabels=labels,sampleNames=colnames(inputExpData),nfold=nfold)
  folds=folds[match(colnames(inputExpData),folds$sampleName),]
  if(!all(colnames(inputExpData)==folds$sampleName)){
  print("error in matching col names!")
  }
  testlabelsList=list(list(status="real_1",testFlag=folds$fold==1,inputLabels=labels))
  if(runReal==TRUE){
  for (i in 1:nfold){
    testlabelsList=c(testlabelsList,list(list(status=paste0("real_",i),testFlag=folds$fold==i,inputLabels=labels)))
  }
  }
  
  
  if(npermTest>0){
    for(i in 1:npermTest){
      tmpLabels=sample(labels,size=length(labels),replace = F)
      folds=.myFoldMakerFn(tmpLabels,colnames(inputExpData),nfold)
      for (j in 1:nfold){
        testlabelsList=c(testlabelsList,list(list(status=paste0("perm",i,"_",j),testFlag=folds$fold==i,inputLabels=tmpLabels)))
      }
    }
  }
  
  testlabelsList=testlabelsList[-1]
  
  outputfilesList=""
  
  ncores=1
  counter=1
  outputfilesList=""
  
  
  for(i in 1:length(testlabelsList)){
    tmplbls=testlabelsList[[i]]
    trainingInput=inputExpData[,!tmplbls$testFlag]
    trainingLabels=tmplbls$inputLabels[!tmplbls$testFlag]
    testInput=inputExpData[,tmplbls$testFlag]
    testLabels=tmplbls$inputLabels[tmplbls$testFlag]
    results=.myInitializer(trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,expClassName = expClassName,method=initializerFns,prevMethod=tmplbls$status)
    for(j in 1:length(results)){
      dir.create(paste0(output.directory,"out",counter))
      data=results[[j]]
      save(data,middleFns,posteriorFns,classificationMethodsList,expClassName,file=paste0(output.directory,"out",counter,"/data.rda"))
      
      outputfilesList=c(outputfilesList,paste0(output.directory,"out",counter,"/"))
      
      counter=counter+1
    }
  }
  
  outputfilesList=outputfilesList[-1]
  write.table(outputfilesList,file=paste0(output.directory,"dataList.txt"),row.names = F,col.names = F,quote=F)
  .myFoldCheckerFn(paste0(output.directory,"dataList.txt"))
}


.myBenchmarkingDataPreparer()


