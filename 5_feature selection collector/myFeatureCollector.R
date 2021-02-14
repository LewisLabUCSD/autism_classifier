#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

rm(list=ls())
setwd("~")
args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/featureSelection/testDataset_HT12/"

.myDataCollector=function(){
  
  path=gsub("\"", "", args[1])
  dataList=read.table(paste0(path,"dataList.txt"))
  path=unlist(strsplit(args,"/"))
  path=path[-length(path)]
  path=paste(path,collapse = "/")
  path=paste0(path,"/")
  res=list()
  for(i in 1:nrow(dataList)){
    if(file.exists(paste0(path,dataList$V1[i],'results.rda'))){
      load(paste0(path,dataList$V1[i],'results.rda'))
      for(j in 1:length(resultsPosterior)){
        if(nrow(resultsPosterior[[j]]$features)!=length(resultsPosterior[[j]]$labels)){
          print("error")
        }
        res=c(res,list(list(features=colnames(resultsPosterior[[j]]$features),method=resultsPosterior[[j]]$method)))
      }
        
    }
  }
  
  save(res,file=paste0(args[1],"feature.rda"))
}


suppressWarnings({
  .myDataCollector()
})
