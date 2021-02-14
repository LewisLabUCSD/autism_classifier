#!/usr/bin/env Rscript

#bokan works
args = commandArgs(trailingOnly=TRUE)

#rm(list=ls())
#setwd("~")
#args='~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/classificationDiscoverySet/'
#args='~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/server_results/independent_dataSets/WG6/independentSetCorrect_complete/'
#args='~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/server_results/independent_dataSets/RNAseq/RNAseqIndependent_complete/'
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test/results/nondeconvoluted_probes/testDataset_HT12/testDataset_HT12/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/WG6_scaled/WG6Scale/newWG6testScale/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/WG6_kernel/WG6Kernel/newWG6test/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_scaled/data/deconvoluted_genes/mainDataset_HT12_WG6/classificationSet/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/mainDataset_HT12_WG6/results/mainDataset_HT12_WG6/classificationSet/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v3/data/deconvoluted_genes/mainDataset_HT12_WG6/results/mainDataset_HT12_WG6/classificationSet/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/LDtestDataset_HT12/testDataset_HT12/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/testDataset_HT12/testDataset_HT12/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2perm1/data/deconvoluted_genes/testDataset_HT12/testDataset_HT12/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/sampleScores/Longitudinaltest/testDataset_HT12/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/sampleScores_correct/mainDataset_HT12_WG6/classificationSet/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/age_balanced/nondeconvoluted_genes/sampleBased/mainDataset_HT12_WG6_v2/classificationSetInternalROC/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/nondeconvoluted_genes/combined_ageBalanced/combined_ageBalanced/classificationSetInternalROC/"
#args="/Users/gazestani/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/combined_ageBalanced/wgcnaFixed/classificationSetInternalROC/"
#args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/combined_ageBalanced/gseaFixed/combined_ageBalanced/classificationSetInternalROC/"
#args="~/classification/ageBalanced_TD39/validation/testDataset_HT12/"
#args="~/classification/ageBalanced_TD39/validation/Longitudinaltest/"
#args="~/classification/ageBalanced_TD39/validation/LDtestDataset_HT12/"
#args="~/classification/ageBalanced_TD39/mainDataset_HT12_WG6_deconv_genes/classificationSet_subject_scores/"
#setwd("~/classificationModuleJabba/classificationCode/3- mainCode/")
#args="~/Desktop/ageBalanced_decov_genes_wF/mainDataset_HT12/mainDataset_HT12/"
#args="testDataset_HT12/"
#args="/data/vahid/classification/ageBalanced_TD39/perm5/classificationSet_subject_scores/"
.myDataCollector=function(args){
  # args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSet/"
  # args="/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/"
  
  require(plyr)
  path=gsub("\"", "", args[1])
  dataList=read.table(paste0(path,"dataList.txt"))
  path=unlist(strsplit(args,"/"))
  path=path[-length(path)]
  path=paste(path,collapse = "/")
  if(length(path)>0){
  path=paste0(path,"/")
  } else {
      path=""
  }
  PRres=""
  # line_count=0
  for(i in 1:nrow(dataList)){
    if(file.exists(paste0(path,dataList$V1[i],'results.rda'))){
      #load(paste0(path,dataList$V1[i],'results.rda'))
      load(paste0(path,dataList$V1[i],'results.rda'))
      if(all(length(PRres)==1&PRres=="")){
        PRres=data.frame(name=row.names(result[["PRresults"]]),result[["PRresults"]])
        ROCres=data.frame(name=row.names(result[["ROCresults"]]),result[["ROCresults"]])
        res85=data.frame(name=row.names(result[["results85"]]),result[["results85"]])
        res9=data.frame(name=row.names(result[["results9"]]),result[["results9"]])
        res95=data.frame(name=row.names(result[["results95"]]),result[["results95"]])
        perc85=data.frame(name=row.names(result[["perc85"]]),result[["perc85"]])
        perc9=data.frame(name=row.names(result[["perc9"]]),result[["perc9"]])
        perc95=data.frame(name=row.names(result[["perc95"]]),result[["perc95"]])
      } else {
        PRres=rbind.fill(PRres,data.frame(name=row.names(result[["PRresults"]]),result[["PRresults"]]))
        ROCres=rbind.fill(ROCres,data.frame(name=row.names(result[["ROCresults"]]),result[["ROCresults"]]))
        res85=rbind.fill(res85,data.frame(name=row.names(result[["results85"]]),result[["results85"]]))
        res9=rbind.fill(res9,data.frame(name=row.names(result[["results9"]]),result[["results9"]]))
        res95=rbind.fill(res95,data.frame(name=row.names(result[["results95"]]),result[["results95"]]))
        perc85=rbind.fill(perc85,data.frame(name=row.names(result[["perc85"]]),result[["perc85"]]))
        perc9=rbind.fill(perc9,data.frame(name=row.names(result[["perc9"]]),result[["perc9"]]))
        perc95=rbind.fill(perc95,data.frame(name=row.names(result[["perc95"]]),result[["perc95"]]))
      }
      # line_count=line_count+1
    }
  }
  
  
  
  if(sum(duplicated(PRres$name))>0){
    tst=as.character(PRres$name[1])
    tst=unlist(strsplit(tst,"_"))
    options(warn = 0)
    if(is.na(as.numeric(tst[2]))){
      PRres$name=.myIdCorrectorFn(PRres$name)
      ROCres$name=.myIdCorrectorFn(ROCres$name)
      res85$name=.myIdCorrectorFn(res85$name)
      res9$name=.myIdCorrectorFn(res9$name)
      res95$name=.myIdCorrectorFn(res95$name)
      perc95$name=.myIdCorrectorFn(perc95$name)
      perc9$name=.myIdCorrectorFn(perc9$name)
      perc85$name=.myIdCorrectorFn(perc85$name)
    }
    options(warn = 1)
  }
  
  row.names(PRres)=PRres$name
  row.names(ROCres)=ROCres$name
  row.names(res85)=res85$name
  row.names(res9)=res9$name
  row.names(res95)=res95$name
  row.names(perc85)=perc85$name
  row.names(perc9)=perc9$name
  row.names(perc95)=perc95$name
  PRres=PRres[,-which(colnames(PRres)=="name")]
  ROCres=ROCres[,-which(colnames(ROCres)=="name")]
  res85=res85[,-which(colnames(res85)=="name")]
  res9=res9[,-which(colnames(res9)=="name")]
  res95=res95[,-which(colnames(res95)=="name")]
  perc85=perc85[,-which(colnames(perc85)=="name")]
  perc9=perc9[,-which(colnames(perc9)=="name")]
  perc95=perc95[,-which(colnames(perc95)=="name")]
  
  PRres=.myMatrixArranger(PRres)
  ROCres=.myMatrixArranger(ROCres)
  res85=.myMatrixArranger(res85)
  res9=.myMatrixArranger(res9)
  res95=.myMatrixArranger(res95)
  perc85=.myMatrixArranger(perc85)
  perc9=.myMatrixArranger(perc9)
  perc95=.myMatrixArranger(perc95)
  
  save(PRres,ROCres,res85,res9,res95,perc85,perc9,perc95,file=paste0(args[1],"ResultsArranged.rda"))
}

.myMatrixArranger=function(inputData){
  inputData[is.na(inputData)]=0
  tmp=strsplit(row.names(inputData),"_")
  options(warn=-1)
  if(!is.na(as.numeric(tmp[[1]][2]))){
    tmp=unlist(lapply(tmp,function(x) paste0(x[-2],collapse = "_")))
  } else {
    tmp=row.names(inputData)
  }
  options(warn = 0)
  
  tmp=data.frame(name=row.names(inputData),category=tmp,stringsAsFactors = F)
  
  resMean=matrix(0,nrow=length(unique(tmp$category)),ncol=ncol(inputData))
  resCount=matrix(0,nrow=length(unique(tmp$category)),ncol=ncol(inputData))
  
  row.names(resMean)=unique(tmp$category)
  row.names(resCount)=unique(tmp$category)
  colnames(resMean)=colnames(inputData)
  colnames(resCount)=colnames(inputData)
  
  if(!all(row.names(inputData)==tmp$name)){
    print("Error")
  }
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(inputData)){
      tmpCount=0
      if(inputData[i,j]!=0){
        tmpCount=1
      }
      resMean[row.names(resMean)==tmp$category[i],colnames(resMean)==colnames(inputData)[j]]=resMean[row.names(resMean)==tmp$category[i],colnames(resMean)==colnames(inputData)[j]]+inputData[i,j]
      resCount[row.names(resMean)==tmp$category[i],colnames(resMean)==colnames(inputData)[j]]=resCount[row.names(resMean)==tmp$category[i],colnames(resMean)==colnames(inputData)[j]]+tmpCount
    }
  }
  
  tmpCount=resCount
  tmpCount[tmpCount==0]=1
  resMean=resMean/tmpCount
  
  if(sum(resMean[resCount==0])>0){
    print("Error")
  }
  
  return(list(resMeans=resMean,resConts=resCount))
}

.myIdCorrectorFn=function(inputIdList){
  inputIdList=as.character(inputIdList)
  tst=as.character(inputIdList[1])
  tst=unlist(strsplit(tst,"_"))
  tst=c(tst[1],"1",tst[2:length(tst)])
  tst=paste(tst,collapse = "_")
  res=tst
  prevId=0
  for(i in 2:length(inputIdList)){
    counter=sum(inputIdList[1:(i-1)] %in% inputIdList[i])+1
    prevId=max(prevId,counter)
    tst=as.character(inputIdList[i])
    tst=unlist(strsplit(tst,"_"))
    tst=c(tst[1],prevId,tst[2:length(tst)])
    tst=paste(tst,collapse = "_")
    res=c(res,tst)
    
  }
  return(res)
}

# suppressWarnings({
  # .myDataCollector(args)
# })
