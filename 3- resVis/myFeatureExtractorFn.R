#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#args=c("inputDataJabba.rda","selectedPathList.rda",5)

source("WGCNAandGeneFilterationMethods.R")
source("pipelines.R")

.myFeatureExtractorRunner=function(){
  library(parallel)
  
  load(args[1])
  load(args[2])
  ncores=args[3]
  
  for(i in 1:length(pathList)){
    if(i==1){
      inputArranged=list(list(path=pathList[i],inputData=inputExpData,inputLabels=labels))
    } else {
      inputArranged=c(inputArranged,list(list(path=pathList[i],inputData=inputExpData,inputLabels=labels)))
    }
  }
  
  results=mclapply(inputArranged,.myFeatureExtractorFn,mc.cores = ncores)
  
  save(results,file="selectedFeatureSets.rda")
}

.myFeatureExtractorFn=function(xInput){
  require(Biobase,quietly = T)
  
  path=xInput$path
  print(path)
  inputData=xInput$inputData
  inputLabels=xInput$inputLabels
  
  #inputData=inputArranged[[1]]$inputData
  #inputLabels=inputArranged[[1]]$inputLabels
  
  path=tolower(path)
  results=list(inputExpData=inputData,inputLabels=inputLabels)
  
  indx=regexpr("_",path)
  preNames=substr(path,0,indx-1)
  affixNames=substr(path,indx+1,nchar(path))
  
  indx=regexpr("_",affixNames)
  while(indx>(-1)){
    method=substr(affixNames,0,indx-1)
    affixNames=substr(affixNames,indx+1,nchar(affixNames))
    indx=regexpr("_",affixNames)
    results=.myFeatureSelSwitch(results,method = method)
    
  }
  
  results=.myFeatureSelSwitch(results,method = affixNames)
  
  return(list(list(name=path,genes=row.names(results$inputExpData))))
}

.myFeatureSelSwitch=function(inputExpSet,method){
  #available methods: cov, var, cor, WGCNA, grn, grn1, grn3, zscore, SIS, lm, selectV, PLSR, CPPLS, PCR, cvplogistic, cvplogisticFast, svm, varImportance, Safs
  switch(method,
         no=.myNoFeatSl(inputExpSet),
         cov=.myCovGeneFilterationFeatSl(inputExpSet = inputExpSet),
         var=.myVarGeneFilterationFeatSl(inputExpSet = inputExpSet),
         cor=.myCorFilterationFeatSl(inputExpSet = inputExpSet),
         grn1=.myGrn1FeatSl(inputExpSet = inputExpSet),
         grn2=.myGrn2FeatSl(inputExpSet = inputExpSet),
         grn3=.myGrn3FeatSl(inputExpSet = inputExpSet),
         zscore=.myZscoreFeatSl(inputExpSet = inputExpSet),
         sis=.mySisFeatSl(inputExpSet = inputExpSet),
         selectvcomexphc=.mySelectVcomExpHCFeatSl(inputExpSet = inputExpSet),
         selectvcomhc=.mySelectVcomHCFeatSl(inputExpSet = inputExpSet),
         selectvcomfair=.mySelectVcomFairFeatSl(inputExpSet = inputExpSet),
         selectvnotcomexphc=.mySelectVnotcomExpHCFeatSl(inputExpSet = inputExpSet),
         selectvnotcomhc=.mySelectVnotcomHCFeatSl(inputExpSet = inputExpSet),
         selectvnotcomfair=.mySelectVnotcomFairFeatSl(inputExpSet = inputExpSet),
         plsr=.myNoFeatSl(inputExpSet = inputExpSet),
         cppls=.myNoFeatSl(inputExpSet = inputExpSet),
         pcr=.myNoFeatSl(inputExpSet = inputExpSet),
         cvplogistic=.myCvplogisticFilterationFn(inputExpSet = inputExpSet,labels = labels),
         cvplogisticfast=.myCvplogisticFastFilterationFn(inputExpSet = inputExpSet,labels = labels),
         svm=.mySVMfilterationFn(inputExpSet = inputExpSet,labels = labels),
         vif=.myVIFfilterationFn(inputExpSet = inputExpSet,labels = labels),
         logisticfwd=.myLogisticFwdFeatSl(inputExpSet = inputExpSet),
         foba=.myFobaFilterationFn(inputExpSet = inputExpSet,labels = labels),
         lm=.myLmFilterationFn(inputExpSet=inputExpSet,labels = labels,outputsize = outputsize),
         varimportance=.myVarImportanceFeatSl(inputExpSet=inputExpSet),
         fafs=.mySafsFilterationFn(inputExpSet=inputExpSet,labels = labels),
         wgcna=.myNoFeatSl(inputExpSet=inputExpSet),
         gseaall=.myGseaAllFeatSl(inputExpSet=inputExpSet),
         gseamed=.myGseaMedFeatSl(inputExpSet=inputExpSet),
         gseamean=.myGseaMeanFeatSl(inputExpSet=inputExpSet))
}

.myCovGeneFilterationFeatSl=function(inputExpSet){
  tmpRes=.myCovGeneFilterationFn(inputExpSet = inputExpSet$inputExpData,testInputData = inputExpSet$inputExpData)
  return((list(inputExpData=tmpRes[[1]],inputLabelsabels=inputExpSet$inputLabels)))
}

.myVarGeneFilterationFeatSl=function(inputExpSet){
  tmpRes=.myVarGeneFilterationFn(inputExpSet = inputExpSet$inputExpData,testInputData = inputExpSet$inputExpData)
  return((list(inputExpData=tmpRes[[1]],inputLabels=inputExpSet$inputLabels)))
}

.myCorFilterationFeatSl=function(inputExpSet){
  tmpRes=.myCorFilterationFn(inputExpSet = inputExpSet$inputExpData)
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes[[1]],],inputLabels=inputExpSet$inputLabels)))
}

.myGrn2FeatSl=function(inputExpSet){
  tmpRes=.myGeneFilteration(inputExpSet$inputExpData,inputExpSet$inputLabels,testInputData=inputExpSet$inputExpData,testLabels=inputExpSet$inputLabels,method="grn2",outputsize=NULL)
  tmpRes=tmpRes[["grn2"]]
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes,],inputLabels=inputExpSet$inputLabels)))
}

.myGrn1FeatSl=function(inputExpSet){
  tmpRes=.myGRNFilterationFn(inputExpSet = inputExpSet$inputExpData,method="grn1")
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes[[1]]$Gene,],inputLabels=inputExpSet$inputLabels)))
}

.myGrn3FeatSl=function(inputExpSet){
  tmpRes=.myGRNFilterationFn(inputExpSet = inputExpSet$inputExpData,method="grn3")
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes[[1]]$Gene,],inputLabels=inputExpSet$inputLabels)))
}

.myZscoreFeatSl=function(inputExpSet){
  tmpRes=.myZscoreFilterationFn(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels)
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes,],inputLabels=inputExpSet$inputLabels)))
}

.myNoFeatSl=function(inputExpSet){
  return(inputExpSet)
}

.mySisFeatSl=function(inputExpSet){
  tmpRes=.mySISFilterationFn(inputExpSet = inputExpSet$inputExpData,labels = inputExpSet$inputLabels)
  return((list(inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% tmpRes,],inputLabels=inputExpSet$inputLabels)))
}

.myGseaAllFeatSl=function(inputExpSet){
  dataP=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData=inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="GSEA")
  if(length(dataP)>1|!is.na(dataP)){
    return((list(inputExpData=inputExpSet$inputExpData[dataP$geneNames,],inputLabels=inputExpSet$inputLabels)))
  
  } else {
    return((list()))
  }
}

.myGseaMedFeatSl=function(inputExpSet){
  dataP=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData=inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="GSEA")
  if(length(dataP)>1|!is.na(dataP)){
    GSEAresMed=dataP[dataP$count>=(summary(dataP$count)[3]),]
    return((list(inputExpData=inputExpSet$inputExpData[GSEAresMed$geneNames,],inputLabels=inputExpSet$inputLabels)))
    
  } else {
    return((list()))
  }
}

.myGseaMeanFeatSl=function(inputExpSet){
  dataP=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData=inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="GSEA")
  if(length(dataP)>1|!is.na(dataP)){
    GSEAresMean=dataP[dataP$count>=(summary(dataP$count)[4]),]
    return((list(inputExpData=inputExpSet$inputExpData[GSEAresMean$geneNames,],inputLabels=inputExpSet$inputLabels)))
    
  } else {
    return((list()))
  }
}

.myWGCNAfeatSl=function(inputExpSet){
  resWgcna=list()
  try({
    msgs=capture.output({
      resWgcna=.myGeneFilteration(inputExpSet=inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData = inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method = "WGCNA")
    })
  }, silent = T)
  
  if(length(resWgcna)>0)
    resWgcna=(list(inputExpData=resWgcna[["eigenGenes"]],inputLabels=inputExpSet$inputLabels))
  
  return(resWgcna)
}

.mySISfeatSl=function(inputExpSet){
  resSIS=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData = inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="SIS")
  resSIS=(list(inputLabels=inputExpSet$labels,inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% resSIS,]))
  return(resSIS)
}

.myLogisticFwdFeatSl=function(inputExpSet){
  resLogFwd=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData = inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="logisticFwd")
  if(!is.na(resLogFwd)){
    indx=regexpr("X",resLogFwd)
    resLogFwd=substr(resLogFwd,indx+1,nchar(resLogFwd))
    resLogFwd=list(inputLabels=inputExpSet$inputLabels,inputExpData=inputExpSet$inputExpData[row.names(inputExpSet$inputExpData) %in% resLogFwd,])
    
    return((resLogFwd))
  } else {
    return(list())
  }
}

.myPcrFeatSl=function(inputExpSet){
  resPcr=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,inputExpSet$inputLabels,testInputData = inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method = "PCR")
  resPcr=list(inputLabels=inputExpSet$inputLabels,inputExpData=resPcr[["results"]])
  return(resPcr)
}

.myPlsrFeatSl=function(inputExpSet){
  resPlsr=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,inputExpSet$inputLabels,inputExpSet$inputExpData,inputExpSet$inputLabels,method = "PLSR")
  resPlsr=list(inputLabels=labels,inputExpData=resPlsr[["results"]])
  return(resPlsr)
}

.myCpplsFeatSl=function(inputExpSet){
  resCppls=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,inputExpSet$inputLabels,inputExpSet$inputExpData,inputExpSet$inputLabels,method = "CPPLS")
  if("results" %in% names(resCppls)){
    resCppls=(list(inputLabels=inputExpSet$inputLabels,inputExpData=resCppls[["results"]]))
  } else {
    resCppls=(list())
  }
  return(resCppls)
}

.mySelectVcomExpHCFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = 'ExpHC')
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.mySelectVcomHCFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = "HC")
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.mySelectVcomFairFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = "Fair")
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.mySelectVnotcomExpHCFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = 'ExpHC',comvar=F)
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.mySelectVnotcomHCFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = "HC",comvar=F)
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.mySelectVnotcomFairFeatSl=function(inputExpSet){
  require(HiDimDA,quietly=T)
  
  expHCcv=SelectV(t(.matrixExtraction(inputExpSet$inputExpData)),as.factor(inputExpSet$inputLabels),Selmethod = "Fair",comvar=F)
  
  return((list(inputExpData=inputExpSet$inputExpData[expHCcv$vkptInd,],inputLabels=inputExpSet$inputLabels)))
}

.myVarImportanceFeatSl=function(inputExpSet){
  res=.myGeneFilteration(inputExpSet = inputExpSet$inputExpData,labels=inputExpSet$inputLabels,testInputData = inputExpSet$inputExpData,testLabels = inputExpSet$inputLabels,method="varImportance")
  return((list(inputLabels=inputExpSet$inputLabels,inputExpData=res$features)))
}

.myFeatureExtractorRunner()
