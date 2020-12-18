#analysis of sample scores
rm(list=ls())
args="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/HT_as_test_v2/data/deconvoluted_genes/sampleScores_correct/mainDataset_HT12_WG6/classificationSet/"

.myResultOrganizerFn=function(inputData,countThr=max(inputData[[2]]),replaceZero=F){
  #results with number of folds less than the countThr are set to zero
  #zeros can be replaced with 0.5 (ROC and PR)
  
  inputData[[1]][inputData[[2]]<countThr]=0
  if(replaceZero){
    inputData[[1]][inputData[[1]]==0]=0.5
  }
  res=inputData[[1]]
  
  thr=0
  if(replaceZero){
    thr=0.5
  }
  x=apply(res,1,function(x) sum(x==thr))
  res=res[which(x<(ncol(inputData[[1]])-1)),]
  return(res)
}

.myWeightFn=function(inputData=1,inputROCres=ROCres,calculateThr=T){
  res=2^(1.5*(inputData-mean(as.vector(inputROCres)))/sd(as.vector(inputROCres)))/2^(1.5*((0.5)-mean(as.vector(inputROCres)))/sd(as.vector(inputROCres)))
  if(calculateThr){
    if(sum(res<quantile(.myWeightFn(as.vector(inputROCres),inputROCres,calculateThr=F),0.5))>0){
      res[res<quantile(.myWeightFn(as.vector(inputROCres),inputROCres,calculateThr=F),0.5)]=0
    }
  }
  
  return(res)
}

.myWeightFn=function(inputData=1,inputROCres=ROCres,calculateThr=T){
  
  
  return(inputData)
}

.myAddSubjectId=function(inputResult,inputPath){
  load(paste0(inputPath,'data.rda'))
  pdata=data.frame(colname=colnames(data$testInputData),pData(data$testInputData))
  
  return(pdata)
}

.myDataCollector=function(){
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
  
  ROCres=""
  resSampleScores=list()
  for(i in 1:nrow(dataList)){
    if(file.exists(paste0(path,dataList$V1[i],'results.rda'))){
      load(paste0(path,dataList$V1[i],'results.rda'))
      tmpAnno=.myAddSubjectId(result,paste0(path,dataList$V1[i]))
      resSampleScores=c(resSampleScores,list(result$ROCdf))
      if(all(length(ROCres)==1&ROCres=="")){
        ROCres=data.frame(name=row.names(result[["ROCresults"]]),result[["ROCresults"]])
        anno=tmpAnno
      } else {
        ROCres=rbind.fill(ROCres,data.frame(name=row.names(result[["ROCresults"]]),result[["ROCresults"]]))
        anno=rbind(anno,tmpAnno)
      }
      rm(result)
    }
  }
  rm(dataList,i)
  anno=anno[!duplicated(anno$colname),]
  
  if(sum(duplicated(ROCres$name))>0){
    tst=as.character(ROCres$name[1])
    tst=unlist(strsplit(tst,"_"))
    options(warn = 0)
    if(is.na(as.numeric(tst[2]))){
      ROCres$name=.myIdCorrectorFn(ROCres$name)
    }
    options(warn = 1)
  }
  
  row.names(ROCres)=ROCres$name
  
  ROCres=ROCres[,-which(colnames(ROCres)=="name")]
  
  
  ROCres=.myMatrixArranger(ROCres)
  
  ROCres=.myResultOrganizerFn(ROCres,replaceZero = T)
  ROCres[ROCres<0.5]=0.5
  ROCres=ROCres-0.5
  
  #################################
  #is main? ie folds are separated by a number?
  isMain=T
  ###########################
  
  ROCroutes=row.names(ROCres)
  
  ROCroutes=which(ROCres>=0.3,arr.ind = T)
  ROCroutes=paste0(row.names(ROCres)[ROCroutes[,1]],"_",colnames(ROCres)[ROCroutes[,2]])
  
  resSamplesCombined=list()
  for(i in 1:length(resSampleScores)){
    print(i)
    for(j in 1:length(resSampleScores[[i]])){
      tmpMethod=resSampleScores[[i]][[j]]$method
      tmpx=unlist(strsplit(tmpMethod,"_"))
      if(isMain){
        tmpx=tmpx[-2]
        tmpMethod=paste(tmpx,collapse = "_")
      }
      tmpx=paste(tmpx[-length(tmpx)],collapse = "_")
      
      
      if(sum(ROCroutes %in% tmpMethod)>0){
        #tmpMethod=paste0(tmpMethod,"_",sample(colnames(ROCres)[-5],1))
        tmpDf=resSampleScores[[i]][[j]][[1]]
        if(length(resSamplesCombined)==0){
          resSamplesCombined=c(resSamplesCombined,list(list(method=tmpMethod,df=resSampleScores[[i]][[i]][[1]])))
        } else{
          if(sum(unlist(lapply(resSamplesCombined,function(x) x$method==tmpMethod)))>0){
            tmpIndex=which(unlist(lapply(resSamplesCombined,function(x) x$method==tmpMethod)))
            tmpDfM=resSamplesCombined[[tmpIndex]]$df
            tmpDfM=merge(tmpDfM,tmpDf[,-2],by="sampleName",all=T)
            resSamplesCombined[[tmpIndex]]$df=tmpDfM
            rm(tmpIndex,tmpDfM)
          } else {
            resSamplesCombined=c(resSamplesCombined,list(list(method=tmpMethod,df=resSampleScores[[i]][[i]][[1]])))
          }
          
        }
        rm(tmpDf)
      }
      rm(tmpMethod)
    }
  }
  #rm(i,j,resSampleScores)
  summary(unlist(lapply(resSamplesCombined,function(x) ncol(x$df))))
  
  
  plot(as.vector(ROCres),.myWeightFn(as.vector(ROCres),ROCres))
  
  for(i in 1:length(resSamplesCombined)){
    tmpDf=resSamplesCombined[[i]]$df
    tmpDf2=data.frame(tmpDf[,-which(colnames(tmpDf) %in% c("sampleName","real","combinedScore"))])
    if(sum(tmpDf2>1,na.rm = T)>0){
      tmpDf2[which(tmpDf2>1,arr.ind = T)]=1
    }
    if(sum(tmpDf2<0,na.rm = T)>0){
      tmpDf2[which(tmpDf2<0,arr.ind = T)]=0
    }
    
    
    tmpDf$combinedScore=rowMeans(tmpDf2,na.rm = T)
    resSamplesCombined[[i]]$df=tmpDf
    rm(tmpDf,tmpDf2)
    
    rowName=unlist(strsplit(resSamplesCombined[[i]]$method,split = "_"))
    colName=rowName[length(rowName)]
    rowName=paste0(rowName[-length(rowName)],collapse = "_")
    if(sum(row.names(ROCres) %in% rowName)>0 & sum(colnames(ROCres) %in% colName)>0){
      resSamplesCombined[[i]]$weight=.myWeightFn(ROCres[rowName,colName],ROCres)
    } else {
      resSamplesCombined[[i]]$weight=.myWeightFn(0,ROCres)
    }
    
  }
  
  weightSum=sum(unlist(lapply(resSamplesCombined,function(x) x$weight)))
  
  resFinalScores=data.frame(sampleName="a",real=-1,weightedSum=0,sumOfWeights=0,stringsAsFactors = F)
  for(i in 1:length(resSamplesCombined)){
    if(resSamplesCombined[[i]]$weight>0){
      #print(i)
      tmpDf=resSamplesCombined[[i]]$df
      tmpDf$scoreWeighted=tmpDf$combinedScore*resSamplesCombined[[i]]$weight
      for(j in 1:nrow(tmpDf)){
        if(sum(resFinalScores$sampleName %in% tmpDf$sampleName[j])>0){
          resFinalScores$weightedSum[resFinalScores$sampleName==tmpDf$sampleName[j]]=(tmpDf$scoreWeighted[j]+resFinalScores$weightedSum[resFinalScores$sampleName==tmpDf$sampleName[j]])
          resFinalScores$sumOfWeights[resFinalScores$sampleName==tmpDf$sampleName[j]]=(resSamplesCombined[[i]]$weight+resFinalScores$sumOfWeights[resFinalScores$sampleName==tmpDf$sampleName[j]])
        } else{
          resFinalScores=rbind(resFinalScores,data.frame(sampleName=tmpDf$sampleName[j],real=tmpDf$real[j],weightedSum=tmpDf$scoreWeighted[j],sumOfWeights=resSamplesCombined[[i]]$weight,stringsAsFactors = F))
        }
      }
    }
    
  }
  resFinalScores=resFinalScores[-1,]
  resFinalScores$score=resFinalScores$weightedSum/resFinalScores$sumOfWeights
  
  resFinalScores=merge(resFinalScores,anno,by.x="sampleName",by.y="colname",all.x=T)
  resFinalScores=resFinalScores[order(resFinalScores$score,decreasing = T),]
  write.table(resFinalScores,file=paste0(path,"resultsArranged_m2.txt"),sep="\t",row.names=F)
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

suppressWarnings({
  .myDataCollector()
})