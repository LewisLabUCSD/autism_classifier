
###filteration pipelines

source("WGCNAandGeneFilterationMethods.R")
#cov, var, cor, WGCNA, grn, zscore, SIS, lm, selectV, PLSR, CPPLS, PCR, cvplogistic, cvplogisticFast, svm, vif, logisticFwd, foba, varImportance, safs

#prior filteration methods: no,cov, var, cov-var,varImportance
#midle filteratin methods: pipelines
#postrior filteration methods: no,wgcna, logisticFwd, SIS, PCR, PLSR

#pipeline 1-GRN
.pipeline1=function(inputData,labels,testInputData,testLabels,expClassName){
  result=list()
  dataGRN=.myGeneFilteration(inputData,labels,testInputData,testLabels,expClassName=expClassName,method="grn")
  
  dataGRN1=dataGRN[["grn1"]]
  dataGRN2=dataGRN[["grn2"]]
  dataGRN3=dataGRN[["grn3"]]
  
  dataGRN1Cntl=dataGRN[["grn1Cntl"]]
  dataGRN2Cntl=dataGRN[["grn2Cntl"]]
  dataGRN3Cntl=dataGRN[["grn3Cntl"]]
  
  result=list()
  {
    #GRN4 (gene with deg>min(deg))
    tmp=as.data.frame(table(c(as.character(dataGRN1$net[,"Fnode"]),as.character(dataGRN1$net[,"Snode"]))))
    tmp=tmp[which(tmp$Freq>(min(tmp$Freq))),]
    dataGRN4=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(tmp$Var1)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(tmp$Var1)),],testLabels=testLabels,method="grn4")
    result=c(result,list(dataGRN4))
    
    tmp=as.data.frame(table(c(as.character(dataGRN1Cntl$net[,"Fnode"]),as.character(dataGRN1Cntl$net[,"Snode"]))))
    tmp=tmp[which(tmp$Freq>(min(tmp$Freq))),]
    dataGRN4Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(tmp$Var1)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(tmp$Var1)),],testLabels=testLabels,method="grn4Cntl")
    result=c(result,list(dataGRN4Cntl))
    
  }
  
  {
    #GRN5 (gene with deg in top 20ile)
    tmp=as.data.frame(table(c(as.character(dataGRN1$net[,"Fnode"]),as.character(dataGRN1$net[,"Snode"]))))
    tmp=tmp[which(tmp$Freq>=(tmp$Freq[nrow(tmp)/5])),]
    dataGRN5=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(tmp$Var1)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(tmp$Var1)),],testLabels=testLabels,method="grn5")
    result=c(result,list(dataGRN5))
    
    tmp=as.data.frame(table(c(as.character(dataGRN1Cntl$net[,"Fnode"]),as.character(dataGRN1Cntl$net[,"Snode"]))))
    tmp=tmp[which(tmp$Freq>=(tmp$Freq[nrow(tmp)/5])),]
    dataGRN5Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(tmp$Var1)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(tmp$Var1)),],testLabels=testLabels,method="grn5Cntl")
    result=c(result,list(dataGRN5Cntl))
    
  }
  
  dataGRN1=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN1$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN1$Gene)),],testLabels=testLabels,method="grn1")
  result=c(result,list(dataGRN1))
  dataGRN1Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN1Cntl$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN1Cntl$Gene)),],testLabels=testLabels,method="grn1Cntl")
  result=c(result,list(dataGRN1Cntl))
  
  
  dataGRN2=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN2$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN2$Gene)),],testLabels=testLabels,method="grn2")
  result=c(result,list(dataGRN2))
  dataGRN2Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN2Cntl$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN2Cntl$Gene)),],testLabels=testLabels,method="grn2Cntl")
  result=c(result,list(dataGRN2Cntl))
  
  
  if(length(dataGRN3$Gene)>0){
    if(length(dataGRN3$Gene)>5){
      dataGRN3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataGRN3$Gene)),],labels=labels,testInputData,testLabels,method="cor")
      
      dataGRN3.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN3.cor)),],testLabels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN3.cor)),],method="grn3_cor")
      result=c(result,list(dataGRN3.cor))
      
    } 
    dataGRN3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN3$Gene)),],testLabels=testLabels,method="grn3")
    result=c(result,list(dataGRN3))
  }
  
  if(length(dataGRN3Cntl$Gene)>0){
    if(length(dataGRN3Cntl$Gene)>5){
      dataGRN3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataGRN3Cntl$Gene)),],labels=labels,testInputData,testLabels,method="cor")
      
      dataGRN3.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN3.cor)),],testLabels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN3.cor)),],method="grn3Cntl_cor")
      result=c(result,list(dataGRN3.cor))
      
    } 
    dataGRN3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataGRN3Cntl$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataGRN3Cntl$Gene)),],testLabels=testLabels,method="grn3Cntl")
    result=c(result,list(dataGRN3))
  }
  
  return(result)
}



#inputExpSet=inputData[dataP,];inputLabels=labels;testInputData=testInputData[dataP,];outputsize=NULL

#pipeline 2-zscore
.pipeline2=function(inputData,labels,testInputData,testLabels,expClassName){
  
  dataP=.myGeneFilteration(inputData,labels,testInputData,testLabels,expClassName=expClassName,method="zscore")
  zscoreRes=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP)),],testLabels=testLabels,method="zscore")
  finalRes=list(zscoreRes)
  
  dataZscore.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP)),],labels=labels,testInputData,testLabels,expClassName=expClassName,method="cor")
  dataZscore.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataZscore.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataZscore.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore_cor")
  finalRes=c(finalRes,list(dataZscore.cor))
  dataZscore.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
  dataZscore.grn3=dataZscore.grn3org[["grn3"]]
  if(length(dataZscore.grn3)>0){
    dataZscore.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataZscore.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataZscore.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore_grn3")
    finalRes=c(finalRes,list(dataZscore.grn3))
  }
  
  if(sum(names(dataZscore.grn3org)=="grn3Cntl")>0){
    dataZscore.grn3=dataZscore.grn3org[["grn3Cntl"]]
    if(length(dataZscore.grn3)>0){
      dataZscore.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataZscore.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataZscore.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore_grn3Cntl")
      finalRes=c(finalRes,list(dataZscore.grn3))
    }
  }
  
  return(finalRes)
}

#pipeline 3-selectV
.pipeline3=function(inputData,labels,testInputData,testLabels,expClassName){
  dataP=.myGeneFilteration(inputData,labels,testInputData,testLabels,expClassName=expClassName,method="selectV")
  
  resComVarExpHc=dataP$resCommonVar$ExpHC$genes
  resComVarHC=dataP$resCommonVar$HC$genes
  resComVarFair=dataP$resCommonVar$Fair$genes
  #resComVarFDR=dataP$resCommonVar$Fdr$genes
  
  resNotComVarExpHC=dataP$ResNotCommonVar$ExpHC$genes
  resNotComVarHC=dataP$ResNotCommonVar$HC$genes
  resNotComVarFair=dataP$ResNotCommonVar$Fair$genes
  #resNotComVarFDR=dataP$ResNotCommonVar$Fdr$genes
  
  resFinal=list()
  
  if(length(resComVarExpHc)>0){
    resComVarExpHcList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc)),],testLabels=testLabels,method="selectVcomExpHC")
    resFinal=list(resComVarExpHcList)
    
    
    resComVarExpHc.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resComVarExpHc.grn3=resComVarExpHc.grn3org[["grn3"]]
    
    if(length(resComVarExpHc.grn3)>0){
      if(nrow(resComVarExpHc.grn3)>5){
        resComVarExpHc.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarExpHc.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3.zscore)),],testLabels=testLabels,method="selectVcomExpHC_grn3_zscore")
        resFinal=c(resFinal,list(resComVarExpHc.grn3.zscore))
        
        resComVarExpHc.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarExpHc.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3.cor)),],testLabels=testLabels,method="selectVcomExpHC_grn3_cor")
        resFinal=c(resFinal,list(resComVarExpHc.grn3.corList))
      }
      
      resComVarExpHc.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,method="selectVcomExpHC_grn3")
      resFinal=c(resFinal,list(resComVarExpHc.grn3))
    }
    
    resComVarExpHc.grn3=resComVarExpHc.grn3org[["grn3Cntl"]]
    if(length(resComVarExpHc.grn3)>0){
      if(nrow(resComVarExpHc.grn3)>5){
        resComVarExpHc.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarExpHc.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3.zscore)),],testLabels=testLabels,method="selectVcomExpHC_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resComVarExpHc.grn3.zscore))
        
        resComVarExpHc.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarExpHc.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3.cor)),],testLabels=testLabels,method="selectVcomExpHC_grn3Cntl_cor")
        resFinal=c(resFinal,list(resComVarExpHc.grn3.corList))
      }
      
      resComVarExpHc.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.grn3$Gene)),],testLabels=testLabels,method="selectVcomExpHC_grn3Cntl")
      resFinal=c(resFinal,list(resComVarExpHc.grn3))
    }
    
    
    resComVarExpHc.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    
    resComVarExpHc.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.cor)),],testLabels=testLabels,method="selectVcomExpHC_cor")
    resFinal=c(resFinal,list(resComVarExpHc.corList))
    
    
    resComVarExpHc.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.cor)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
    resComVarExpHc.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarExpHc.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarExpHc.cor.zscore)),],testLabels=testLabels,method="selectVcomExpHC_cor_zscore")
    resFinal=c(resFinal,list(resComVarExpHc.cor.zscore))
  }
  ######################
  if(length(resComVarHC)>0){
    resComVarHCList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC)),],testLabels=testLabels,method="selectVcomHC")
    resFinal=c(resFinal,list(resComVarHCList))
    
    resComVarHC.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resComVarHC.grn3=resComVarHC.grn3org[["grn3"]]
    if(length(resComVarHC.grn3)>0){
      if(nrow(resComVarHC.grn3)>5){
        resComVarHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3.cor)),],testLabels=testLabels,method="selectVcomHC_grn3_cor")
        resFinal=c(resFinal,list(resComVarHC.grn3.corList))
        
        resComVarHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3.zscore)),],testLabels=testLabels,method="selectVcomHC_grn3_zscore")
        resFinal=c(resFinal,list(resComVarHC.grn3.zscore))
      }
      resComVarHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomHC_grn3")
      resFinal=c(resFinal,list(resComVarHC.grn3))
    }
    
    resComVarHC.grn3=resComVarHC.grn3org[["grn3Cntl"]]
    if(length(resComVarHC.grn3)>0){
      if(nrow(resComVarHC.grn3)>5){
        resComVarHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3.cor)),],testLabels=testLabels,method="selectVcomHC_grn3Cntl_cor")
        resFinal=c(resFinal,list(resComVarHC.grn3.corList))
        
        resComVarHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3.zscore)),],testLabels=testLabels,method="selectVcomHC_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resComVarHC.grn3.zscore))
      }
      resComVarHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomHC_grn3Cntl")
      resFinal=c(resFinal,list(resComVarHC.grn3))
    }
    
    resComVarHC.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    
    resComVarHC.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomHC_cor")
    resFinal=c(resFinal,list(resComVarHC.corList))
    
    
    resComVarHC.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarHC.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.cor)),],testLabels=testLabels,labels,expClassName=expClassName,method="zscore")
    resComVarHC.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarHC.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarHC.cor.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomHC_cor_zscore")
    resFinal=c(resFinal,list(resComVarHC.cor.zscore))
  }
  ########################
  if(length(resComVarFair)>0){
    resComVarFairList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair")
    resFinal=c(resFinal,list(resComVarFairList))
    
    resComVarFair.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resComVarFair.grn3=resComVarFair.grn3org[["grn3"]]
    if(length(resComVarFair.grn3)>0){
      if(nrow(resComVarFair.grn3)>5){
        
        resComVarFair.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarFair.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3_cor")
        resFinal=c(resFinal,list(resComVarFair.grn3.corList))
        
        resComVarFair.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarFair.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3_zscore")
        resFinal=c(resFinal,list(resComVarFair.grn3.zscore))
      }
      resComVarFair.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3")
      resFinal=c(resFinal,list(resComVarFair.grn3))
    }
    
    resComVarFair.grn3=resComVarFair.grn3org[["grn3Cntl"]]
    if(length(resComVarFair.grn3)>0){
      if(nrow(resComVarFair.grn3)>5){
        
        resComVarFair.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resComVarFair.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3Cntl_cor")
        resFinal=c(resFinal,list(resComVarFair.grn3.corList))
        
        resComVarFair.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resComVarFair.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resComVarFair.grn3.zscore))
      }
      resComVarFair.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_grn3Cntl")
      resFinal=c(resFinal,list(resComVarFair.grn3))
    }
    
    resComVarFair.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resComVarFair.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_cor")
    resFinal=c(resFinal,list(resComVarFair.corList))
    
    resComVarFair.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resComVarFair.cor)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
    resComVarFair.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resComVarFair.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resComVarFair.cor.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVcomFair_cor_zscore")
    resFinal=c(resFinal,list(resComVarFair.cor.zscore))
  }
  ##########################
  ##########################
  if(length(resNotComVarExpHC)>0){
    resNotComVarExpHCList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC")
    resFinal=c(resFinal,list(resNotComVarExpHCList))
    
    resNotComVarExpHC.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resNotComVarExpHC.grn3=resNotComVarExpHC.grn3org[["grn3"]]
    if(length(resNotComVarExpHC.grn3)>0){
      
      if(nrow(resNotComVarExpHC.grn3)>5){
        resNotComVarExpHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resNotComVarExpHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3_cor")
        resFinal=c(resFinal,list(resNotComVarExpHC.grn3.corList))
        
        resNotComVarExpHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resNotComVarExpHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3_zscore")
        resFinal=c(resFinal,list(resNotComVarExpHC.grn3.zscore))
      }
      resNotComVarExpHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3")
      resFinal=c(resFinal,list(resNotComVarExpHC.grn3))
    }
    
    resNotComVarExpHC.grn3=resNotComVarExpHC.grn3org[["grn3Cntl"]]
    if(length(resNotComVarExpHC.grn3)>0){
      
      if(length(resNotComVarExpHC.grn3)>5){
        resNotComVarExpHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resNotComVarExpHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3Cntl_cor")
        resFinal=c(resFinal,list(resNotComVarExpHC.grn3.corList))
        
        resNotComVarExpHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resNotComVarExpHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resNotComVarExpHC.grn3.zscore))
      }
      resNotComVarExpHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_grn3Cntl")
      resFinal=c(resFinal,list(resNotComVarExpHC.grn3))
    }
    
    resNotComVarExpHC.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resNotComVarExpHC.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_cor")
    resFinal=c(resFinal,list(resNotComVarExpHC.corList))
    
    resNotComVarExpHC.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.cor)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
    resNotComVarExpHC.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarExpHC.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarExpHC.cor.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComExpHC_cor_zscore")
    resFinal=c(resFinal,list(resNotComVarExpHC.cor.zscore))
  }
  ######################
  if(length(resNotComVarHC)>0){
    resNotComVarHCList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC")
    resFinal=c(resFinal,list(resNotComVarHCList))
    
    resNotComVarHC.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resNotComVarHC.grn3=resNotComVarHC.grn3org[["grn3"]]
    if(length(resNotComVarHC.grn3)>0){
      if(nrow(resNotComVarHC.grn3)>5){
        resNotComVarHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,labels=labels,expClassName=expClassName,method="cor")
        resNotComVarHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3_cor")
        resFinal=c(resFinal,list(resNotComVarHC.grn3.corList))
        
        resNotComVarHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,labels=labels,expClassName=expClassName,method="zscore")
        resNotComVarHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3_zscore")
        resFinal=c(resFinal,list(resNotComVarHC.grn3.zscore))
      }
      resNotComVarHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3")
      resFinal=c(resFinal,list(resNotComVarHC.grn3))
    }
    
    resNotComVarHC.grn3=resNotComVarHC.grn3org[["grn3Cntl"]]
    if(length(resNotComVarHC.grn3)>0){
      if(length(resNotComVarHC.grn3)>5){
        resNotComVarHC.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,labels=labels,expClassName=expClassName,method="cor")
        resNotComVarHC.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3Cntl_cor")
        resFinal=c(resFinal,list(resNotComVarHC.grn3.corList))
        
        resNotComVarHC.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,labels=labels,expClassName=expClassName,method="zscore")
        resNotComVarHC.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resNotComVarHC.grn3.zscore))
      }
      resNotComVarHC.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_grn3Cntl")
      resFinal=c(resFinal,list(resNotComVarHC.grn3))
    }
    
    resNotComVarHC.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resNotComVarHC.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_cor")
    resFinal=c(resFinal,list(resNotComVarHC.corList))
    
    resNotComVarHC.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.cor)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
    resNotComVarHC.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarHC.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarHC.cor.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComHC_cor_zscore")
    resFinal=c(resFinal,list(resNotComVarHC.cor.zscore))
  }
  ########################
  if(length(resNotComVarFair)>0){
    resNotComVarFairList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair")
    resFinal=c(resFinal,list(resNotComVarFairList))
    
    resNotComVarFair.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resNotComVarFair.grn3=resNotComVarFair.grn3org[["grn3"]]
    
    if(length(resNotComVarFair.grn3)>0){
      if(nrow(resNotComVarFair.grn3)>5){
        resNotComVarFair.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resNotComVarFair.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3_cor")
        resFinal=c(resFinal,list(resNotComVarFair.grn3.corList))
        
        resNotComVarFair.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resNotComVarFair.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3_zscore")
        resFinal=c(resFinal,list(resNotComVarFair.grn3.zscore))
      }
      resNotComVarFair.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3")
      resFinal=c(resFinal,list(resNotComVarFair.grn3))
    }
    
    resNotComVarFair.grn3=resNotComVarFair.grn3org[["grn3Cntl"]]
    
    if(length(resNotComVarFair.grn3)>0){
      if(nrow(resNotComVarFair.grn3)>5){
        resNotComVarFair.grn3.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="cor")
        resNotComVarFair.grn3.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3Cntl_cor")
        resFinal=c(resFinal,list(resNotComVarFair.grn3.corList))
        
        resNotComVarFair.grn3.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
        resNotComVarFair.grn3.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3Cntl_zscore")
        resFinal=c(resFinal,list(resNotComVarFair.grn3.zscore))
      }
      resNotComVarFair.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_grn3Cntl")
      resFinal=c(resFinal,list(resNotComVarFair.grn3))
    }
    
    resNotComVarFair.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resNotComVarFair.corList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.cor)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_cor")
    resFinal=c(resFinal,list(resNotComVarFair.corList))
    
    resNotComVarFair.cor.zscore=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.cor)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.cor)),],testLabels=testLabels,expClassName=expClassName,method="zscore")
    resNotComVarFair.cor.zscore=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resNotComVarFair.cor.zscore)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resNotComVarFair.cor.zscore)),],testLabels=testLabels,expClassName=expClassName,method="selectVnotComFair_cor_zscore")
    resFinal=c(resFinal,list(resNotComVarFair.cor.zscore))
  }
  return(resFinal)
}

#pipeline 4-svm
.pipeline4=function(inputData,labels,testInputData,testLabels,expClassName){
  dataP=.myGeneFilteration(inputData,labels,testInputData=testInputData,testLabels=testLabels,method="svm")
  resultScad=list(labels=labels,features=inputData[dataP$resultScad,],testInputData=testInputData[dataP$resultScad,],testLabels=testLabels,method="svmScad")
  resultL1=list(labels=labels,features=inputData[dataP$resultL1,],testInputData=testInputData[dataP$resultL1,],testLabels=testLabels,method="svmL1")
  resultElasticNet=list(labels=labels,features=inputData[dataP$resultElasticNet,],testInputData=testInputData[dataP$resultElasticNet,],testLabels=testLabels,method="svmElasticNet")
  resultScadL2=list(labels=labels,features=inputData[dataP$resultScadL2,],testInputData=testInputData[dataP$resultScadL2,],testLabels=testLabels,method="svmScadL2")
  return(list(resultScad,resultL1,resultScadL2,resultElasticNet))
}

#pipeline 5-GSEA
.pipeline5=function(inputData,labels,testInputData,testLabels,expClassName){
  
  dataPall=.myGeneFilteration(inputData,labels=labels,testInputData=testInputData,testLabels = testLabels,expClassName=expClassName,method="GSEA")
  resList=list()
  dataP=dataPall[[1]]
  if(length(dataP)>1|sum(!is.na(dataP))>0){
    GSEAres=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall")
    
    dataP=dataP[order(dataP$count,decreasing = T),]
    
    GSEADecile=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/10)],],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/10)],],testLabels=testLabels,expClassName=expClassName,method="GSEADecile")
    
    GSEA20ile=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/5)],],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/5)],],testLabels=testLabels,expClassName=expClassName,method="GSEA20ile")
    
    GSEAresMed=dataP[dataP$count>=(summary(dataP$count)[3]),]
    GSEAresMedList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(GSEAresMed$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMed$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmed")
    
    GSEAresMean=dataP[dataP$count>=(summary(dataP$count)[4]),]
    GSEAresMeanList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(GSEAresMean$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMean$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmean")
    
    resGSEAall.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAall.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall_cor")
    
    resGSEAMean.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(GSEAresMean$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMean$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAMean.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAMean.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAMean.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmean_cor")
    
    resGSEAMed.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(GSEAresMed$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMed$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAMed.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAMed.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAMed.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmed_cor")
    
    
    resGSEAall.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resGSEAall.grn3=resGSEAall.grn3org[["grn3"]]
    resGSEAall.grn3Cntl=resGSEAall.grn3org[["grn3Cntl"]]
    if(length(resGSEAall.grn3)>0 & length(resGSEAall.grn3Cntl)>0){
      resGSEAall.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall_grn3")
      resGSEAall.grn3Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3Cntl$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3Cntl$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall_grn3Cntl")
      
      resList=c(resList,list(resGSEAall.cor,resGSEAall.grn3,resGSEAall.grn3Cntl,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor,GSEADecile,GSEA20ile))
    } else {
      if(length(resGSEAall.grn3)>0){
        resGSEAall.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall_grn3")
        
        resList=c(resList,list(resGSEAall.cor,resGSEAall.grn3,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor))
      } else {
        resList= c(resList,list(resGSEAall.cor,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor,GSEADecile,GSEA20ile)) 
      }
      
    }
  }
  
  dataP=dataPall[[2]]
  if(length(dataP)>1|sum(!is.na(dataP))>0){
    GSEAres=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall0.1")
    
    dataP=dataP[order(dataP$count,decreasing = T),]
    
    GSEADecile=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/10)],],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/10)],],testLabels=testLabels,expClassName=expClassName,method="GSEADecile0.1")
    
    GSEA20ile=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/5)],],testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames))[1:(nrow(dataP)/5)],],testLabels=testLabels,expClassName=expClassName,method="GSEA20ile0.1")
    
    GSEAresMed=dataP[dataP$count>=(summary(dataP$count)[3]),]
    GSEAresMedList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(GSEAresMed$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMed$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmed0.1")
    
    GSEAresMean=dataP[dataP$count>=(summary(dataP$count)[4]),]
    GSEAresMeanList=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(GSEAresMean$geneNames)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMean$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmean0.1")
    
    resGSEAall.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAall.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall0.1_cor")
    
    resGSEAMean.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(GSEAresMean$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMean$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAMean.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAMean.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAMean.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmean0.1_cor")
    
    resGSEAMed.cor=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(GSEAresMed$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(GSEAresMed$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="cor")
    resGSEAMed.cor=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAMed.cor)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAMed.cor)),],testLabels=testLabels,expClassName=expClassName,method="GSEAmed0.1_cor")
    
    
    resGSEAall.grn3org=.myGeneFilteration(inputData[which(row.names(inputData) %in% as.character(dataP$geneNames)),],labels=labels,testInputData=testInputData[which(row.names(testInputData) %in% as.character(dataP$geneNames)),],testLabels=testLabels,expClassName=expClassName,method="grn3")
    resGSEAall.grn3=resGSEAall.grn3org[["grn3"]]
    resGSEAall.grn3Cntl=resGSEAall.grn3org[["grn3Cntl"]]
    if(length(resGSEAall.grn3)>0 & length(resGSEAall.grn3Cntl)>0){
      resGSEAall.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall0.1_grn3")
      resGSEAall.grn3Cntl=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3Cntl$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3Cntl$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall0.1_grn3Cntl")
      
      resList=c(resList,list(resGSEAall.cor,resGSEAall.grn3,resGSEAall.grn3Cntl,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor,GSEADecile,GSEA20ile))
    } else {
      if(length(resGSEAall.grn3)>0){
        resGSEAall.grn3=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(resGSEAall.grn3$Gene)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(resGSEAall.grn3$Gene)),],testLabels=testLabels,expClassName=expClassName,method="GSEAall0.1_grn3")
        
        resList=c(resList,list(resGSEAall.cor,resGSEAall.grn3,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor))
      } else {
        resList= c(resList,list(resGSEAall.cor,GSEAres,GSEAresMeanList,GSEAresMedList,resGSEAMean.cor,resGSEAMed.cor,GSEADecile,GSEA20ile)) 
      }
      
    }
  }
  
  return(resList)
}


#pipeline 6-DE analysis
.pipeline6=function(inputData,labels,testInputData,testLabels,expClassName){
  
  dataDE=.myGeneFilteration(inputExpSet=inputData,labels,testInputData,testLabels,expClassName=expClassName,method="DE",outputsize=NULL)
  dataDE=dataDE[[1]]
  resFinal=list()
  if(sum(dataDE$adj.P.Val<0.05,na.rm = T)>0){
    FDRall=row.names(dataDE)[dataDE$adj.P.Val<0.05]
    resFDRall=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(FDRall)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(FDRall)),],testLabels=testLabels,method="FDRall")
    resFinal=c(resFinal,list(resFDRall))
  }
  
  if(sum(dataDE$adj.P.Val<0.05 & dataDE$logFC>0,na.rm = T)>0){
    FDRup=row.names(dataDE)[dataDE$adj.P.Val<0.05 & dataDE$logFC>0]
    resFDRup=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(FDRup)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(FDRup)),],testLabels=testLabels,method="FDRup")
    resFinal=c(resFinal,list(resFDRup))
  }
  
  if(sum(dataDE$adj.P.Val<0.05 & dataDE$logFC<0,na.rm = T)>0){
    FDRdown=row.names(dataDE)[dataDE$adj.P.Val<0.05 & dataDE$logFC<0]
    resFDRdown=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(FDRdown)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(FDRdown)),],testLabels=testLabels,method="FDRdown")
    resFinal=c(resFinal,list(resFDRdown))
  }
  
  if(sum(dataDE$adj.P.Val<0.01,na.rm = T)>0){
    Pvalall=row.names(dataDE)[dataDE$adj.P.Val<0.01]
    resPvalall=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(Pvalall)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(Pvalall)),],testLabels=testLabels,method="Pvalall")
    resFinal=c(resFinal,list(resPvalall))
  }
  
  if(sum(dataDE$adj.P.Val<0.01 & dataDE$logFC>0,na.rm = T)>0){
    Pvalup=row.names(dataDE)[dataDE$adj.P.Val<0.01 & dataDE$logFC>0]
    resPvalup=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(Pvalup)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(Pvalup)),],testLabels=testLabels,method="Pvalup")
    resFinal=c(resFinal,list(resPvalup))
  }
  
  if(sum(dataDE$adj.P.Val<0.01 & dataDE$logFC<0,na.rm=T)>0){
    Pvaldown=row.names(dataDE)[dataDE$adj.P.Val<0.01 & dataDE$logFC<0]
    resPvaldown=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(Pvaldown)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(Pvaldown)),],testLabels=testLabels,method="Pvaldown")
    resFinal=c(resFinal,list(resPvaldown))
  }
  
  dataDE=dataDE[order(dataDE$P.Value,decreasing = F),]
  if(length(row.names(dataDE))>100){
    top100all=row.names(dataDE)[1:100]
    res100all=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top100all)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top100all)),],testLabels=testLabels,method="top100all")
    resFinal=c(resFinal,list(res100all))
  }
  
  if(length(row.names(dataDE))>500){
    top500all=row.names(dataDE)[1:500]
    resTop500all=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top500all)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top500all)),],testLabels=testLabels,method="top500all")
    resFinal=c(resFinal,list(resTop500all))
  }
  
  if(length(row.names(dataDE))>1000){
    top1000all=row.names(dataDE)[1:1000]
    resTop1000all=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top1000all)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top1000all)),],testLabels=testLabels,method="top1000all")
    resFinal=c(resFinal,list(resTop1000all))
  }
  
  if(length(row.names(dataDE[dataDE$logFC>0,]))>100){
    top100up=row.names(dataDE[dataDE$logFC>0,])[1:100]
    resTop100up=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top100up)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top100up)),],testLabels=testLabels,method="top100up")
    resFinal=c(resFinal,list(resTop100up))
  }
  
  
  if(length(row.names(dataDE[dataDE$logFC>0,]))>500){
    top500up=row.names(dataDE[dataDE$logFC>0,])[1:500]
    resTop500up=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top500up)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top500up)),],testLabels=testLabels,method="top500up")
    resFinal=c(resFinal,list(resTop500up))
  }
  
  
  if(length(row.names(dataDE[dataDE$logFC>0,]))>1000){
    top1000up=row.names(dataDE[dataDE$logFC>0,])[1:1000]
    resTop1000up=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top1000up)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top1000up)),],testLabels=testLabels,method="top1000up")
    resFinal=c(resFinal,list(resTop1000up))
  }
  
  
  if(length(row.names(dataDE[dataDE$logFC<0,]))>100){
    top100down=row.names(dataDE[dataDE$logFC<0,])[1:100]
    resTop100down=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top100down)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top100down)),],testLabels=testLabels,method="top100down")
    resFinal=c(resFinal,list(resTop100down))
  }
  
  
  if(length(row.names(dataDE[dataDE$logFC<0,]))>500){
    top500down=row.names(dataDE[dataDE$logFC<0,])[1:500]
    resTop500down=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top500down)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top500down)),],testLabels=testLabels,method="top500down")
    resFinal=c(resFinal,list(resTop500down))
  }
  
  
  if(length(row.names(dataDE[dataDE$logFC<0,]))>1000){
    top1000down=row.names(dataDE[dataDE$logFC<0,])[1:1000]
    resTop1000down=list(labels=labels,features=inputData[which(row.names(inputData) %in% as.character(top1000down)),],testInputData=testInputData[which(row.names(testInputData) %in% as.character(top1000down)),],testLabels=testLabels,method="top1000down")
    resFinal=c(resFinal,list(resTop1000down))
  }
  
  return(resFinal)
}

.myFeatureSelectionBenchmarking=function(data,labels,testInputData,testLabels,ncores=1,prevMethod,initializerFns=list("no","cov","var","cov_var","varImportance"),middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline4","pipeline5","pipeline6"),posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls'),expClassName){
  .myPackageInstaller()
  .myPackageChecker()
  results=list();
  results=.myInitializer(data,labels=labels,testInputData=testInputData,testLabels=testLabels,method=initializerFns,prevMethod)
  resultsMiddle=.myMiddleFn(results,methodList=middleFns,ncores=ncores,expClassName=expClassName)
  resultsPosterior=.myPosteriorFn(resultsMiddle,methodList=posteriorFns,ncores=ncores)
  return(resultsPosterior)
}

.myPackageInstaller=function(){
  packagesRequired=c("plyr","ddalpha","parallel","SIS","HiDimDA","pls","parallel","cvplogistic","penalizedSVM","VIF","foba","leaps","caret","flashClust","WGCNA","car","PRROC","glmnet","randomForest","gbm","e1071","nnet","pROC","SQDA","car")
  bioconductorPackages=c("BiocInstaller","minet","Biobase","edgeR","AnnotationDbi", "impute", "GO.db", "preprocessCore")
  
  source("http://bioconductor.org/biocLite.R")
  for(i in 1:length(bioconductorPackages)){
    if((bioconductorPackages[i] %in% rownames(installed.packages())) == FALSE) {biocLite(bioconductorPackages[i])}
  }
  
  for(i in 1:length(packagesRequired)){
    if((packagesRequired[i] %in% rownames(installed.packages())) == FALSE) {install.packages(packagesRequired[i])}
  }
}

.myPackageChecker=function(){
  packagesRequired=c("plyr","parallel","SIS","HiDimDA","pls","parallel","cvplogistic","penalizedSVM","VIF","foba","leaps","caret","flashClust","WGCNA","SQDA","car")
  bioconductorPackages=c("minet","Biobase","edgeR","AnnotationDbi", "impute", "GO.db", "preprocessCore")
  
  for(i in 1:length(bioconductorPackages)){
    if((bioconductorPackages[i] %in% rownames(installed.packages())) == FALSE) {print(paste0("needs to be installed: ",bioconductorPackages[i]))}
  }
  
  for(i in 1:length(packagesRequired)){
    if((packagesRequired[i] %in% rownames(installed.packages())) == FALSE) {print(paste0("needs to be installed: ",packagesRequired[i]))}
  }
}

.myInitializer=function(data,labels,testInputData,testLabels,expClassName,method=list("no","cov","var","cov_var","varImportance"),prevMethod=""){
  print("features data would be in expressionSet format")
  res=list();
  for(i in 1:length(method)){
    tmpRes=.myInitializerSwitch(data,labels=labels,testInputData = testInputData,testLabels=testLabels,expClassName = expClassName,method = method[[i]])
    tmpRes[[1]]$method=paste0(prevMethod,"_",tmpRes[[1]]$method)
    if(i==1)
      res=tmpRes
    else
      res=c(res,tmpRes)
  }
  return(res) 
}

#inputData=list(data)
#methodsList=methodList=middleFns
#ncores=1
#expClassName=expClassName
#j=1
#i=3

.myMiddleFn=function(inputData,methodList=list("no","pipeline1","pipeline2","pipeline3","pipeline5","pipeline6"),ncores=1,expClassName){
  print("features data would be in expressionSet format")
  
  ncores=min(ncores,length(inputData))
  
  #require(parallel,quietly=T)
  if(ncores>1){
    x=split(1:length(inputData), cut(1:length(inputData), ncores))
  } else {
    x=1:length(inputData)
  }
  #res=mclapply(x,.myMiddleFnDetailed,inputData,methodList,mc.preschedule=F,mc.cores=ncores)
  res=lapply(x,.myMiddleFnDetailed,inputData,methodList,expClassName=expClassName)
  result=res[[1]]
  if(length(res)>1){
    for(ic in 2:length(res))
      result=c(result,res[[ic]])
  }
  tst=lapply(result,function(x) dim(x$features)[1])
  tst=unlist(tst)
  tst=tst>0
  result=result[tst]
  
  tst=lapply(result,function(x) dim(x$testInputData)[1])
  tst=unlist(tst)
  tst=tst>0
  result=result[tst]
  return(result)
}

.myMiddleFnDetailed=function(x,inputData,methodsList,expClassName){
  result=list()
  for(j in x){
    labels=inputData[[j]]$labels
    dataP=inputData[[j]]$features
    testInputData=inputData[[j]]$testInputData
    testLabels=inputData[[j]]$testLabels
    method=inputData[[j]]$method
    #pipeline3=.pipeline3(dataP,labels=labels,testInputData=testInputData,testLabels=testLabels)
    for(i in 1:length(methodsList)){
      tmpRes=list()
      try({
      tmpRes=.myMiddleFnSwitch(inputData=dataP,labels=labels,testInputData = testInputData,testLabels=testLabels,method=methodsList[[i]],expClassName=expClassName)
      },silent = T)
      print(paste0(length(tmpRes)," results were collected from ",method,"_",methodsList[[i]]," route"))
      if(length(tmpRes)>0){
        for(counter in 1:length(tmpRes)){
          tmpMethod=paste0(method,"_",tmpRes[[counter]][["method"]])
          tmpList=list(list(labels=labels,features=tmpRes[[counter]]$features,testInputData=tmpRes[[counter]]$testInputData,testLabels=tmpRes[[counter]]$testLabels,method=tmpMethod))
          if(length(result)==0){
            result=tmpList
          } else {
            result=c(result,tmpList)
          }
        }
      }
    }
  }
  
  return(result)
}

#inputData=resultsMiddle
#methodList=posteriorFns;ncores=1;expClassName=expClassName

.myPosteriorFn=function(inputData,methodList=list("no","pcr","plsr","cppls","logisticFwd","wgcna","sis"),ncores=1,expClassName){
  print("feature data would be transposed[nrow(data)=length(labels)]")
  require(parallel,quietly=T)
  result=list()
  ncores=min(ncores,length(inputData))
  indSelParMethods=methodList %in% c("logisticFwd","sis")
  
  if(sum(indSelParMethods)>0){
    if(ncores>1&length(inputData)>1){
      x=1:length(inputData)
    } else {
      if(length(inputData)>0){
        x=1:length(inputData)
      } else {
        x=NA
      }
    }
    if(!is.na(x)){
      res=lapply(x,.myPosteriorFnDetailed,inputData,methodList[indSelParMethods],expClassName=expClassName)
      #res=lapply(x,.myPosteriorFnDetailed,inputData,methodList[indSelParMethods])
      result=res[[1]]
      if(length(res)>1){
        for(i in 2:length(res)){
          result=c(result,res[[i]])
        }
      }
    }
  }
  
  if(sum(!indSelParMethods)>0){
    res=.myPosteriorFnDetailed(1:length(inputData),inputData,methodList=methodList[!indSelParMethods],expClassName=expClassName)
    if(length(result)==0){
      result=res
    } else {
      result=c(result,res)
    }
  }
  return(result)
}

.myPosteriorFnDetailed=function(x,inputData,methodList,expClassName){
  results=list()
  dfErrors=data.frame(inputId=0,method="a",stringsAsFactors = F)
  #dfNew=data.frame(fcounter=0,scounter=0,internalcounter=0,method="",stringsAsFactors = F)
  for(i in x){
    tmpData=inputData[[i]]$features
    tmpLabels=inputData[[i]]$labels
    tmpMethod=inputData[[i]]$method
    tmpTstInputData=inputData[[i]]$testInputData
    tmpTstLabels=inputData[[i]]$testLabels
    for(j in 1:length(methodList)){
      tmpRes=list()
      try({
      tmpRes=.myPosteriorFnSwitch(tmpData,tmpLabels,testInputData = tmpTstInputData,testLabels = tmpTstLabels,expClassName=expClassName,methodList[[j]])
      },silent = T)
      print(paste0(length(tmpRes)," results were collected from ",tmpMethod,"_",methodList[[j]]," route"))
      if(length(tmpRes)>0){
          tmpResMethod=paste0(tmpMethod,"_",methodList[[j]])
          tmpResList=list(list(labels=tmpLabels,features=tmpRes[[1]]$features,testInputData=tmpRes[[1]]$testInputData,testLabels=tmpRes[[1]]$testLabels,method=tmpResMethod))
          #dfNew=rbind(dfNew,c(i,j,counter,tmpResMethod))
          if(length(results)==0){
            results=tmpResList
          } else {
            results=c(results,tmpResList)
          }
      } else {
        dfErrors=rbind(dfErrors,c(i,methodList[j]))
      }
    }
  }
  
  return(results)
}

.myInitializerSwitch=function(inputData,labels,testInputData,testLabels,expClassName,method){
  switch(method,
         no=.myNoWrapperFn(inputData,labels=labels,testInputData = testInputData,testLabels=testLabels,expClassName=expClassName),
         cov=.myCovWrapperFn(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName),
         var=.myVarWrapperFn(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName),
         cov_var=.myCovVarWrapperFn(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName),
         varImportance=.myVarImportanceWrapperFn(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName))
}

.myMiddleFnSwitch=function(inputData,labels,testInputData,testLabels,method,expClassName){
  switch(method,
         no=.myNoWrapperFn(inputData,labels = labels,testInputData=testInputData,testLabels = testLabels,expClassName=expClassName),
         pipeline1=.pipeline1(inputData,labels=labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName),
         pipeline2=.pipeline2(inputData,labels=labels,testInputData = testInputData,testLabels=testLabels,expClassName=expClassName),
         pipeline3=.pipeline3(inputData,labels=labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName),
         pipeline4=.pipeline4(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName),
         pipeline5=.pipeline5(inputData,labels=labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName),
         pipeline6=.pipeline6(inputData,labels=labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))
}

.myPosteriorFnSwitch=function(inputData,labels,testInputData,testLabels,expClassName,method){
  switch(method,
         no={return(.myPosteriorNoFn(inputData,labels = labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         wgcna={return(.myWGCNAwrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         pcr={return(.myPCRwrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         plsr={return(.myPLSRwrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         cppls={return(.myCPPLSwrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         sis={return(.mySISwrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))},
         logisticFwd={return(.myLogisticFwdWrapperFn(inputData,labels,testInputData=testInputData,testLabels=testLabels,expClassName=expClassName))})
}



.myWGCNAwrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  
  res=list(list(labels=labels,features=t(.matrixExtraction(inputData)),testInputData=t(.matrixExtraction(testInputData)),testLabels=testLabels,method="wgcna"))
  return(res)
}


.mySISwrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  resSIS=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="SIS")
  resSIS=list(list(labels=labels,features=t(.matrixExtraction(inputData)[resSIS,]),testInputData=t(.matrixExtraction(testInputData)[resSIS,]),testLabels=testLabels,method="SIS"))
  return(resSIS)
}

.myLogisticFwdWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  resLogFwd=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="logisticFwd")
  if(!is.na(resLogFwd)){
    indx=regexpr("X",resLogFwd)
    resLogFwd=substr(resLogFwd,indx+1,nchar(resLogFwd))
    resLogFwd=list(labels=labels,features=t(.matrixExtraction(inputData)[resLogFwd,]),testInputData=t(.matrixExtraction(testInputData)[resLogFwd,]),testLabels=testLabels,expClassName=expClassName,method="logFwd")
    
    return(list(resLogFwd))
  } else {
    return(list())
  }
}

.myPCRwrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=list(list(labels=labels,features=t(.matrixExtraction(inputData)),testInputData=t(.matrixExtraction(testInputData)),testLabels=testLabels,method="pcr"))
  return(res)
}

.myPLSRwrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=list(list(labels=labels,features=t(.matrixExtraction(inputData)),testInputData=t(.matrixExtraction(testInputData)),testLabels=testLabels,method="plsr"))
  return(res)
}

.myCPPLSwrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=list(list(labels=labels,features=t(.matrixExtraction(inputData)),testInputData=t(.matrixExtraction(testInputData)),testLabels=testLabels,method="cppl"))
  return(res)
}

.myPosteriorNoFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=list(list(labels=labels,features=t(.matrixExtraction(inputData)),testInputData=t(.matrixExtraction(testInputData)),testLabels=testLabels,method="no"))
  return(res)
}

.myNoWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  return(list(list(labels=labels,features=inputData,testInputData=testInputData,testLabels=testLabels,method="no")))
}

.myCovWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="cov")
  return(list(list(labels=labels,features=res$inputExpSet,testInputData=res$testInputData,testLabels=testLabels,method="cov")))
}

.myVarImportanceWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="varImportance")
  return(list(list(labels=labels,features=res$features,testInputData=res$testInputData,testLabels=testLabels,method="varImportance")))
}

.myVarWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName=expClassName){
  res=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="var")
  return(list(list(labels=labels,features=res$inputExpSet,testInputData=res$testInputData,testLabels=testLabels,method="var")))
}

.myCovVarWrapperFn=function(inputData,labels,testInputData,testLabels,expClassName){
  res=.myGeneFilteration(inputData,labels=labels,testInputData = testInputData,testLabels = testLabels,expClassName=expClassName,method="cov")
  res=.myGeneFilteration(res$inputExpSet,labels=labels,testInputData = res$testInputData,testLabels = testLabels,expClassName=expClassName,method="var")
  return(list(list(labels=labels,features=res$inputExpSet,testInputData=res$testInputData,testLabels=testLabels,method="cov_var")))
}


