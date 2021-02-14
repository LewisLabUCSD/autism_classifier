#New definition of training and test sets:
##batch 1 and WG6 are used for training
##batch 2 is used as the test 

#batch correction was done only on ASD, TD, and LD samples

#Note that the data could be generated with or without collapsing the probe level data to gene level
#Note that the data can be deconvoluted or not to be deconvoluted

library(Biobase)
.myFoldCheckerFn=function(inputPath){
  require(data.table)
  x=read.table(inputPath,stringsAsFactors=F,sep="\t")
  
  
  for(i in 1:nrow(x)){
    load(paste0(x$V1[i],"data.rda"))
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
    load(paste0(x$V1[i],"data.rda"))
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

.myFinalAdderFn=function(inputData){
  final_ADOS_CoSoTot=rep(NA,nrow(inputData))
  final_ADOS_CoSoTot[is.na(inputData$ados_CoSoTot_5)]=inputData$ados_CoSoTot_4[is.na(inputData$ados_CoSoTot_5)]
  final_ADOS_CoSoTot[is.na(inputData$ados_CoSoTot_4)]=inputData$ados_CoSoTot_3[is.na(inputData$ados_CoSoTot_4)]
  final_ADOS_CoSoTot[is.na(inputData$ados_CoSoTot_3)]=inputData$ados_CoSoTot_2[is.na(inputData$ados_CoSoTot_3)]
  final_ADOS_CoSoTot[!is.na(inputData$ados_CoSoTot_5)]=inputData$ados_CoSoTot_5[!is.na(inputData$ados_CoSoTot_5)]
  final_ADOS_CoSoTot[is.na(inputData$ados_CoSoTot_2)]=inputData$ados_CoSoTot_1[is.na(inputData$ados_CoSoTot_2)]
  inputData$final_ADOS_CoSoTot=final_ADOS_CoSoTot
  
  final_ADOS_RRTot=rep(NA,nrow(inputData))
  final_ADOS_RRTot[is.na(inputData$ados_RRTot_5)]=inputData$ados_RRTot_4[is.na(inputData$ados_RRTot_5)]
  final_ADOS_RRTot[is.na(inputData$ados_RRTot_4)]=inputData$ados_RRTot_3[is.na(inputData$ados_RRTot_4)]
  final_ADOS_RRTot[is.na(inputData$ados_RRTot_3)]=inputData$ados_RRTot_2[is.na(inputData$ados_RRTot_3)]
  final_ADOS_RRTot[!is.na(inputData$ados_RRTot_5)]=inputData$ados_RRTot_5[!is.na(inputData$ados_RRTot_5)]
  final_ADOS_RRTot[is.na(inputData$ados_RRTot_2)]=inputData$ados_RRTot_1[is.na(inputData$ados_RRTot_2)]
  if(sum(is.na(final_ADOS_RRTot))>0){
    print("Some ADOS RR scores are missing")
  }
  inputData$final_ADOS_RRTot=final_ADOS_RRTot
  
  final_ADOS_CoSoRRTot=rep(NA,nrow(inputData))
  final_ADOS_CoSoRRTot[is.na(inputData$ados_CoSoTotRRTot_5)]=inputData$ados_CoSoTotRRTot_4[is.na(inputData$ados_CoSoTotRRTot_5)]
  final_ADOS_CoSoRRTot[is.na(inputData$ados_CoSoTotRRTot_4)]=inputData$ados_CoSoTotRRTot_3[is.na(inputData$ados_CoSoTotRRTot_4)]
  final_ADOS_CoSoRRTot[is.na(inputData$ados_CoSoTotRRTot_3)]=inputData$ados_CoSoTotRRTot_2[is.na(inputData$ados_CoSoTotRRTot_3)]
  final_ADOS_CoSoRRTot[!is.na(inputData$ados_CoSoTotRRTot_5)]=inputData$ados_CoSoTotRRTot_5[!is.na(inputData$ados_CoSoTotRRTot_5)]
  final_ADOS_CoSoRRTot[is.na(inputData$ados_CoSoTotRRTot_2)]=inputData$ados_CoSoTotRRTot_1[is.na(inputData$ados_CoSoTotRRTot_2)]
  if(sum(is.na(final_ADOS_CoSoRRTot))>0){
    print("Some ADOS CoSoRR scores are missing")
  }
  inputData$final_ADOS_CoSoRRTot=final_ADOS_CoSoRRTot
  
  final_mullen_ELC=rep(NA,nrow(inputData))
  final_mullen_ELC[is.na(inputData$mullen_ELC_Std_5)]=inputData$mullen_ELC_Std_4[is.na(inputData$mullen_ELC_Std_5)]
  final_mullen_ELC[is.na(inputData$mullen_ELC_Std_4)]=inputData$mullen_ELC_Std_3[is.na(inputData$mullen_ELC_Std_4)]
  final_mullen_ELC[is.na(inputData$mullen_ELC_Std_3)]=inputData$mullen_ELC_Std_2[is.na(inputData$mullen_ELC_Std_3)]
  final_mullen_ELC[!is.na(inputData$mullen_ELC_Std_5)]=inputData$mullen_ELC_Std_5[!is.na(inputData$mullen_ELC_Std_5)]
  final_mullen_ELC[is.na(inputData$mullen_ELC_Std_2)]=inputData$mullen_ELC_Std_1[is.na(inputData$mullen_ELC_Std_2)]
  final_mullen_ELC[final_mullen_ELC==0]=NA
  if(sum(is.na(final_mullen_ELC))>0){
    print("Some mullen_ELC scores are missing")
  }
  inputData$final_mullen_ELC=final_mullen_ELC
  
  #adding final DSM version
  final_DSM=rep(NA,nrow(inputData))
  final_DSM[is.na(inputData$DxJ_DSM_Version_5)]=inputData$DxJ_DSM_Version_4[is.na(inputData$DxJ_DSM_Version_5)]
  final_DSM[is.na(inputData$DxJ_DSM_Version_4)]=inputData$DxJ_DSM_Version_3[is.na(inputData$DxJ_DSM_Version_4)]
  final_DSM[is.na(inputData$DxJ_DSM_Version_3)]=inputData$DxJ_DSM_Version_2[is.na(inputData$DxJ_DSM_Version_3)]
  final_DSM[!is.na(inputData$DxJ_DSM_Version_5)]=inputData$DxJ_DSM_Version_5[!is.na(inputData$DxJ_DSM_Version_5)]
  final_DSM[is.na(inputData$DxJ_DSM_Version_2)]=inputData$DxJ_DSM_Version_1[is.na(inputData$DxJ_DSM_Version_2)]
  inputData$final_DSM=final_DSM
  
  final_vine_AdapBehav_DomStd=rep(NA,nrow(inputData))
  final_vine_AdapBehav_DomStd[is.na(inputData$vine_AdapBehav_DomStd_5)]=inputData$vine_AdapBehav_DomStd_4[is.na(inputData$vine_AdapBehav_DomStd_5)]
  final_vine_AdapBehav_DomStd[is.na(inputData$vine_AdapBehav_DomStd_4)]=inputData$vine_AdapBehav_DomStd_3[is.na(inputData$vine_AdapBehav_DomStd_4)]
  final_vine_AdapBehav_DomStd[is.na(inputData$vine_AdapBehav_DomStd_3)]=inputData$vine_AdapBehav_DomStd_2[is.na(inputData$vine_AdapBehav_DomStd_3)]
  final_vine_AdapBehav_DomStd[!is.na(inputData$vine_AdapBehav_DomStd_5)]=inputData$vine_AdapBehav_DomStd_5[!is.na(inputData$vine_AdapBehav_DomStd_5)]
  final_vine_AdapBehav_DomStd[is.na(inputData$vine_AdapBehav_DomStd_2)]=inputData$vine_AdapBehav_DomStd_1[is.na(inputData$vine_AdapBehav_DomStd_2)]
  if(sum(is.na(final_vine_AdapBehav_DomStd))>0){
    print("Some vine_AdapBehav_DomStd scores are missing")
  }
  inputData$final_vine_AdapBehav_DomStd=final_vine_AdapBehav_DomStd
  
  final_vine_MtrGross_AgeEq_mo=rep(NA,nrow(inputData))
  final_vine_MtrGross_AgeEq_mo[is.na(inputData$vine_MtrGross_AgeEq_mo_5)]=(inputData$vine_MtrGross_AgeEq_mo_4-inputData$vine_agemo_4)[is.na(inputData$vine_MtrGross_AgeEq_mo_5)]
  final_vine_MtrGross_AgeEq_mo[is.na(inputData$vine_MtrGross_AgeEq_mo_4)]=(inputData$vine_MtrGross_AgeEq_mo_3-inputData$vine_agemo_3)[is.na(inputData$vine_MtrGross_AgeEq_mo_4)]
  final_vine_MtrGross_AgeEq_mo[is.na(inputData$vine_MtrGross_AgeEq_mo_3)]=(inputData$vine_MtrGross_AgeEq_mo_2-inputData$vine_agemo_2)[is.na(inputData$vine_MtrGross_AgeEq_mo_3)]
  final_vine_MtrGross_AgeEq_mo[!is.na(inputData$vine_MtrGross_AgeEq_mo_5)]=(inputData$vine_MtrGross_AgeEq_mo_5-inputData$vine_agemo_5)[!is.na(inputData$vine_MtrGross_AgeEq_mo_5)]
  final_vine_MtrGross_AgeEq_mo[is.na(inputData$vine_MtrGross_AgeEq_mo_2)]=(inputData$vine_MtrGross_AgeEq_mo_1-inputData$vine_agemo_1)[is.na(inputData$vine_MtrGross_AgeEq_mo_2)]
  if(sum(is.na(final_vine_MtrGross_AgeEq_mo))>0){
    print("Some vine_MtrGross_AgeEq_mo scores are missing")
  }
  inputData$final_vine_AdapBehav_DomStd=final_vine_AdapBehav_DomStd
  
  final_mullen_RLT=rep(NA,nrow(inputData))
  final_mullen_RLT[is.na(inputData$mullen_RLT_5)]=(inputData$mullen_RLT_4)[is.na(inputData$mullen_RLT_5)]
  final_mullen_RLT[is.na(inputData$mullen_RLT_4)]=(inputData$mullen_RLT_3)[is.na(inputData$mullen_RLT_4)]
  final_mullen_RLT[is.na(inputData$mullen_RLT_3)]=(inputData$mullen_RLT_2)[is.na(inputData$mullen_RLT_3)]
  final_mullen_RLT[is.na(inputData$mullen_RLT_2)]=(inputData$mullen_RLT_1)[is.na(inputData$mullen_RLT_2)]
  final_mullen_RLT[!is.na(inputData$mullen_RLT_5)]=(inputData$mullen_RLT_5)[!(is.na(inputData$mullen_RLT_5))]
  if(sum(is.na(final_mullen_RLT))>0){
    print("Some mullen_RLT scores are missing")
  }
  inputData$final_mullen_RLT=final_mullen_RLT
  
  final_mullen_ELT=rep(NA,nrow(inputData))
  final_mullen_ELT[is.na(inputData$mullen_ELT_5)]=(inputData$mullen_ELT_4)[is.na(inputData$mullen_ELT_5)]
  final_mullen_ELT[is.na(inputData$mullen_ELT_4)]=(inputData$mullen_ELT_3)[is.na(inputData$mullen_ELT_4)]
  final_mullen_ELT[is.na(inputData$mullen_ELT_3)]=(inputData$mullen_ELT_2)[is.na(inputData$mullen_ELT_3)]
  final_mullen_ELT[is.na(inputData$mullen_ELT_2)]=(inputData$mullen_ELT_1)[is.na(inputData$mullen_ELT_2)]
  final_mullen_ELT[!is.na(inputData$mullen_ELT_5)]=(inputData$mullen_ELT_5)[!(is.na(inputData$mullen_ELT_5))]
  if(sum(is.na(final_mullen_ELT))>0){
    print("Some mullen_ELT scores are missing")
  }
  inputData$final_mullen_ELT=final_mullen_ELT
  
  diagnosis_Seq=inputData$DxJ_DxGroup_1
  diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_2)]=paste0(diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_2)],";",inputData$DxJ_DxGroup_2[!is.na(inputData$DxJ_DxGroup_2)])
  diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_3)]=paste0(diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_3)],";",inputData$DxJ_DxGroup_3[!is.na(inputData$DxJ_DxGroup_3)])
  diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_4)]=paste0(diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_4)],";",inputData$DxJ_DxGroup_4[!is.na(inputData$DxJ_DxGroup_4)])
  diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_5)]=paste0(diagnosis_Seq[!is.na(inputData$DxJ_DxGroup_5)],";",inputData$DxJ_DxGroup_5[!is.na(inputData$DxJ_DxGroup_5)])
  inputData$final_diagnosis_Seq=diagnosis_Seq
  
  inputData=inputData[,c("subjectid",colnames(inputData)[grepl("final",colnames(inputData))],colnames(inputData)[grepl("recent",colnames(inputData))])]
  return(inputData)
  
}
.myDxNormalizingFn=function(inputData,inputColName){
  # inputData=LWdata
  # inputColName ='DxJ_DxGroup_3'
  if(sum(regexpr("Tests Typica",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Tests Typica",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("Test Typica",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Test Typica",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("Prev Tested",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Prev Tested",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("PrevDDTyp",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("PrevDDTyp",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("Prev tested",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Prev tested",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("Tests typical",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Tests typical",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("tests typical",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("tests typical",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("PrevLDDTyp",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("PrevLDDTyp",inputData[,inputColName])>(-1)),inputColName]="TD"
  }
  
  if(sum(regexpr("Preemie No Delay",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Preemie No Delay",inputData[,inputColName])>(-1)),inputColName]="PND"
  }
  
  if(sum(regexpr("MD",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("MD",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  if(sum(regexpr("Premie Delay",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Premie Delay",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  if(sum(regexpr("FMD",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("FMD",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  if(sum(regexpr("PDDNOS",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("PDDNOS",inputData[,inputColName])>(-1)),inputColName]="ASD"
  }
  
  if(sum(regexpr("Typ Sib ASD",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Typ Sib ASD",inputData[,inputColName])>(-1)),inputColName]="TypSibASD"
  }
  
  if(sum(regexpr("DD",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("DD",inputData[,inputColName])>(-1)),inputColName]="GDD"
  }
  
  if(sum(regexpr("ASD Features",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("ASD Features",inputData[,inputColName])>(-1)),inputColName]="ASD Features"
  }
  
  if(sum(regexpr("Fragile X",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Fragile X",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  if(sum(regexpr("Aspergers",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("Aspergers",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  if(sum(regexpr("ADHD",inputData[,inputColName])>(-1),na.rm = T)>0){
    inputData[which(regexpr("ADHD",inputData[,inputColName])>(-1)),inputColName]="Other"
  }
  
  return(inputData)
}
.myDxReaderFn=function(){
  # LWdata=read.table("/Volumes/Work/Vahid_work/classification_newcode_data/final_result_plot/LWReport_04182020_updated.csv",row.names=1,sep=",",header=T,stringsAsFactors = F,comment.char = "",dec=NULL,quote ="",na.strings = "NULL")
  LWdata <- read.csv("/Volumes/Work/Vahid_work/classification_newcode_data/final_result_plot/LWReport_04182020_updated.csv")
  LWdata=.myDxNormalizingFn(LWdata,"DxJ_DxGroup_1")
  LWdata=.myDxNormalizingFn(LWdata,"DxJ_DxGroup_2")
  LWdata=.myDxNormalizingFn(LWdata,"DxJ_DxGroup_3")
  LWdata=.myDxNormalizingFn(LWdata,"DxJ_DxGroup_4")
  LWdata=.myDxNormalizingFn(LWdata,"DxJ_DxGroup_5")
  LWdata=.myDxNormalizingFn(LWdata,"recentDxJ_dxCode")
  LWdata=LWdata[!(LWdata$recentDxJ_dxCode %in% c("Dropped","DROPPED","DROPPRD","No Dx; Abondaned","P7P6D","Test typical; sig fam psychopathology bipolar; depression; OCD")),]
  LWdata=.myFinalAdderFn(LWdata)
  return(LWdata)
}

.myBatchCorrectionFn=function(inputData,doCollapsing=F,doDeconvolution=T,usedDxBinary=F){
  require(sva)
  require(CellCODE)
  require(biomaRt)
  require(limma)
  
  inputData=inputData[,inputData$diagnosis_binary %in% c("ASD","LD","TD")]
  inputData$diagnosis_binary=as.factor(as.character(inputData$diagnosis_binary))
  
  #batch = pData(inputData)$batch
  #modcombat = model.matrix(~1, data=pData(inputData))
  #modexperiment = model.matrix(~diagnosis_binary, data=pData(inputData))
  #combat_edata = ComBat(dat=exprs(inputData), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)
  #tmp=inputData
  #exprs(tmp)=combat_edata
  
  model = model.matrix(~0+diagnosis_binary + as.factor(batch)+as.factor(sex),data=pData(inputData))
  fit = lmFit(inputData,design=model)
  beta1 <- fit$coefficients[, 4:6, drop = FALSE]
  beta1[is.na(beta1)] <- 0
  expAdj=exprs(inputData) - beta1 %*% t(model[,4:6])
  tmp=inputData
  exprs(tmp)=expAdj
  
  
  
  
  #selecting the top 1000 most variable probes based on variance metric
  selProbes=apply(exprs(tmp),1,var)
  selProbes=data.frame(name=names(selProbes),val=selProbes)
  o=order(selProbes,decreasing = T)
  selProbes=selProbes[o,]
  selProbes=selProbes$name[1:1000]
  selProbes=as.character(selProbes)
  
  #clustering of samples based on top 1000 most variable probes
  hc=hclust(as.dist((1-cor(exprs(tmp)[row.names(tmp) %in% selProbes,]))/2),method="average")
  x=colnames(tmp)[hc$order]
  
  #removal of potential outlier
  tmp=tmp[,-which(colnames(tmp) %in% c(x[1:2],x[(length(x)-2):length(x)]))]
  
  fData(tmp)$Entrez_Gene_ID=fData(tmp)$ENTREZ_GENE_ID
  
  if(doCollapsing){
    x=rowMeans(exprs(tmp))
    o=order(x,decreasing=T)
    tmp=tmp[o,]
    tmp=tmp[!duplicated(fData(tmp)$Entrez_Gene_ID),]
    tmp=tmp[!is.na(fData(tmp)$Entrez_Gene_ID),]
    row.names(tmp)=fData(tmp)$Entrez_Gene_ID
  }
  
  if(doDeconvolution){
    mart=useMart("ENSEMBL_MART_ENSEMBL")
    hensemble=useDataset("hsapiens_gene_ensembl",mart)
    gsIds=getBM(attributes = c("entrezgene","external_gene_name")
                ,filters="entrezgene",values=unique(fData(tmp)$Entrez_Gene_ID),mart=hensemble)
    
    gsIds=gsIds[!duplicated(gsIds$entrezgene),]
    tmp2=data.frame(entrezgene=fData(tmp)$Entrez_Gene_ID,id=row.names(tmp),stringsAsFactors = F)
    tmp2=merge(tmp2,gsIds,by="entrezgene",all.x=T)
    
    tmp=tmp[order(row.names(tmp)),]
    tmp2=tmp2[order(tmp2$id),]
    if(all(row.names(tmp)==tmp2$id)){
      fData(tmp)$geneSymbol=tmp2$external_gene_name
    } else {
      print("error in gene symbol mapping")
    }
    
    tmp2=tmp
    o=rowMeans(exprs(tmp2))
    tmp2=tmp2[order(o,decreasing = T),]
    tmp2=tmp2[!duplicated(fData(tmp2)$geneSymbol)]
    tmp2=tmp2[!is.na(fData(tmp2)$geneSymbol),]
    row.names(tmp2)=fData(tmp2)$geneSymbol
    
    data("IRIS")
    irisTag=tagData(IRIS[,c("Neutrophil-Resting","CD4Tcell-N0", 
                            "Monocyte-Day0", "Bcell-naïve",
                            "NKcell-control", "PlasmaCell-FromPBMC", "DendriticCell-LPSstimulated")],2, max=50, 
                    ref=tmp2, ref.mean=F);
    colnames(irisTag)=c("Neutrophil","Tcell", "Monocyte", 
                        "Bcell", "NKcell", "PlasmaCell", "DendriticCell" )
    if(usedDxBinary){
      SPVs=getAllSPVs(exprs(tmp2), pData(tmp2)$diagnosis_binary, irisTag, "mixed", T) 
    } else {
      SPVs=getAllSPVs(exprs(tmp2), pData(tmp2)$diagnosis_multi1, irisTag, "mixed", T) 
    }
    
    
    if(all(colnames(tmp)==colnames(tmp2))){
      pData(tmp)=cbind(pData(tmp),as.data.frame(SPVs))
    } else {
      print("Error in adding SPVs")
    }
    
    
    model= model.matrix(~0+diagnosis_binary + Neutrophil+Tcell+Monocyte+Bcell+NKcell,data=pData(tmp))
    colnames(model)=unlist(lapply(colnames(model),function(x) gsub(" ","",x)))
    #colnames(model)=c("ASD","LD","TD", "Neutrophil", "Tcell", "Monocyte", "Bcell","NKcell")
    
    fit = lmFit(tmp,design=model)
    beta1 <- fit$coefficients[, (ncol(model)-4):ncol(model), drop = FALSE]
    beta1[is.na(beta1)] <- 0
    expAdj=exprs(tmp) - beta1 %*% t(model[,(ncol(model)-4):ncol(model)])
    exprs(tmp)=expAdj
    
  }
  
  tmp=tmp[,tmp$diagnosis_binary %in% c("ASD","LD","TD")]
  tmp$diagnosis_binary=droplevels(as.factor(tmp$diagnosis_binary))
  
  return(tmp)
}

.myMainTrainingFn=function(inputData,outPath=""){
  inputData=inputData[,inputData$diagnosis_binary %in% c("TD","ASD","proband")]
  inputExpData=inputData
  labels=as.character(pData(inputData)$diagnosis_binary)
  labels[labels=="ASD"]="proband"
  
  expClassName="proband";
  ncores=20;
  nfold=5;
  npermTest=0;
  runReal=TRUE #run data without label permutation
  posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls')
  #middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5")
  middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5",'pipeline6')
  initializerFns=list("no","cov","var","cov_var","varImportance")
  classificationMethodsList=list("reg","logReg","lda","qda","sqda","ridgeReg","lassoReg","ridgeLogReg","lassoLogReg","elasticNetLogReg","randomForest","boosting","bagging")
  save(inputExpData,labels,expClassName,ncores,nfold,npermTest,posteriorFns,middleFns,initializerFns,classificationMethodsList,runReal,file=paste0(outPath,"inputDataJabba.rda"))
}

.myTestWriterFn=function(inputTraining,inputTest,outPath=""){
  inputTraining=inputTraining[,inputTraining$diagnosis_binary %in% c("TD","ASD","proband")]
  inputTest=inputTest[,inputTest$diagnosis_binary %in% c("TD","ASD","proband")]
  
  labelsTraining=as.character(pData(inputTraining)$diagnosis_binary)
  labelsTraining[labelsTraining=="ASD"]="proband"
  
  labelsTest=as.character(pData(inputTest)$diagnosis_binary)
  labelsTest[labelsTest=="ASD"]="proband"
  
  expClassName="proband";
  ncores=20;
  nfold=5;
  npermTest=1;
  runReal=TRUE #run data without label permutation
  posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls')
  #middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5")
  middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5",'pipeline6')
  initializerFns=list("no","cov","var","cov_var","varImportance")
  classificationMethodsList=list("reg","logReg","lda","qda","sqda","ridgeReg","lassoReg","ridgeLogReg","lassoLogReg","elasticNetLogReg","randomForest","boosting","bagging")
  
  results=.myInitializer(inputTraining,labels=labelsTraining,testInputData=inputTest,testLabels=labelsTest,expClassName = "proband",prevMethod="independent")
  
  ncores=1
  counter=1
  outputfilesList=""
  for(j in 1:length(results)){
    dir.create(paste0(outPath,"out",counter))
    data=results[[j]]
    save(data,middleFns,posteriorFns,classificationMethodsList,expClassName,file=paste0(outPath,"out",counter,"/data.rda"))
    
    outputfilesList=c(outputfilesList,paste0(outPath,"out",counter,"/"))
    
    counter=counter+1
  }
  
  outputfilesList=outputfilesList[-1]
  write.table(outputfilesList,file=paste0(outPath,"dataList.txt"),row.names = F,col.names = F,quote=F)
  .myFoldCheckerFn(paste0(outPath,"dataList.txt"))
}

.myLDTestWriterFn=function(inputTraining,inputTest,outPath=""){
  inputTraining=inputTraining[,inputTraining$diagnosis_binary %in% c("TD","ASD","proband")]
  inputTest=inputTest[,inputTest$diagnosis_binary %in% c("TD","LD")]
  
  labelsTraining=as.character(pData(inputTraining)$diagnosis_binary)
  labelsTraining[labelsTraining=="ASD"]="proband"
  
  labelsTest=as.character(pData(inputTest)$diagnosis_binary)
  labelsTest[labelsTest=="LD"]="proband"
  
  expClassName="proband";
  ncores=20;
  nfold=5;
  npermTest=1;
  runReal=TRUE #run data without label permutation
  posteriorFns=list('no','wgcna', 'logisticFwd', 'sis', 'pcr', 'plsr','cppls')
  #middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5")
  middleFns=list("no","pipeline1","pipeline2","pipeline3","pipeline5")
  initializerFns=list("no","cov","var","cov_var","varImportance")
  classificationMethodsList=list("reg","logReg","lda","qda","sqda","ridgeReg","lassoReg","ridgeLogReg","lassoLogReg","elasticNetLogReg","randomForest","boosting","bagging")
  
  results=.myInitializer(inputTraining,labels=labelsTraining,testInputData=inputTest,testLabels=labelsTest,expClassName = "proband",prevMethod="independent")
  
  ncores=1
  counter=1
  outputfilesList=""
  for(j in 1:length(results)){
    dir.create(paste0(outPath,"out",counter))
    data=results[[j]]
    save(data,middleFns,posteriorFns,classificationMethodsList,expClassName,file=paste0(outPath,"out",counter,"/data.rda"))
    
    outputfilesList=c(outputfilesList,paste0(outPath,"out",counter,"/"))
    
    counter=counter+1
  }
  
  outputfilesList=outputfilesList[-1]
  write.table(outputfilesList,file=paste0(outPath,"dataList.txt"),row.names = F,col.names = F,quote=F)
  .myFoldCheckerFn(paste0(outPath,"dataList.txt"))
}

setwd("~/classificationModuleJabba/classificationCode/3- mainCode/")
source("WGCNAandGeneFilterationMethods.R")
source("pipelines.R")

setwd("~")

rm(list=ls())

load("~/OneDrive - UC San Diego/vahid_mac/Documents/data/dropbox_shared_tiziano/WG6/WG6_HT12_Complete.rda")

lwdata=.myDxReaderFn()
lwdata=merge(pData(dataQuantile),lwdata,by.x="subjectId",by.y="subjectid",all.x=T)

dataQuantile=dataQuantile[,order(dataQuantile$subjectId)]
lwdata=lwdata[order(lwdata$subjectId),]

all(lwdata$subjectId==dataQuantile$subjectId)
pData(dataQuantile)$diagnosis_binary=lwdata$recentDxJ_dxCode

#data could be generated with or without collapsing the probes to genes
dataQuantileCollapsedDeconvoluted=.myBatchCorrectionFn(inputData=dataQuantile,doCollapsing = T,doDeconvolution = T,usedDxBinary=T)


table(dataQuantileCollapsedDeconvoluted$batch,dataQuantileCollapsedDeconvoluted$diagnosis_binary)
rm(dataQuantile)

data1WGcollapsedDeconvoluted=dataQuantileCollapsedDeconvoluted[,dataQuantileCollapsedDeconvoluted$batch %in% c("1","WG6")]
data2collapsedDeconvoluted=dataQuantileCollapsedDeconvoluted[,dataQuantileCollapsedDeconvoluted$batch==2]

data1WGcollapsedDeconvoluted=data1WGcollapsedDeconvoluted[,-which(data1WGcollapsedDeconvoluted$subjectId %in% data2collapsedDeconvoluted$subjectId & data1WGcollapsedDeconvoluted$diagnosis_binary=="ASD")]

tmp=data1WGcollapsedDeconvoluted[,order(data1WGcollapsedDeconvoluted$age,decreasing = T)]
tmp=tmp[,tmp$diagnosis_binary %in% c("ASD","TD")]
tmp=tmp[,!(tmp$subjectId %in% data2collapsedDeconvoluted$subjectId & tmp$diagnosis_binary=="ASD")]
tmp=tmp[,-which(duplicated(tmp$subjectId)& tmp$diagnosis_binary=="ASD")]
tmp=tmp[,order(tmp$age,decreasing = F)]
tmp=tmp[,-which(duplicated(tmp$subjectId) & tmp$diagnosis_binary=="TD")]

dupSbj=data1WGcollapsedDeconvoluted$subjectId[duplicated(data1WGcollapsedDeconvoluted$subjectId)]
dupSbj=unique(as.character(dupSbj))
data1WGcollapsedDeconvoluted_c=data1WGcollapsedDeconvoluted[,data1WGcollapsedDeconvoluted$subjectId %in% dupSbj]
data1WGcollapsedDeconvoluted_c=data1WGcollapsedDeconvoluted_c[,(!colnames(data1WGcollapsedDeconvoluted_c) %in% colnames(tmp))]
data1WGcollapsedDeconvoluted=tmp
rm(dupSbj,tmp)

tmp=data2collapsedDeconvoluted[,order(data2collapsedDeconvoluted$age,decreasing = T)]
tmp=tmp[,!tmp$subjectId %in% data1WGcollapsedDeconvoluted$subjectId]
tmp=tmp[,-which(duplicated(data2collapsedDeconvoluted$subjectId)& data2collapsedDeconvoluted$diagnosis_binary=="ASD")]
tmp=tmp[,order(tmp$age,decreasing = F)]
tmp=tmp[,-which(duplicated(tmp$subjectId) & tmp$diagnosis_binary %in% c("TD","LD"))]
dupSbj=data2collapsedDeconvoluted$subjectId[duplicated(data2collapsedDeconvoluted$subjectId)]
dupSbj=unique(as.character(dupSbj))
data2collapsedDeconvoluted_c=data2collapsedDeconvoluted[,data2collapsedDeconvoluted$subjectId %in% dupSbj & (!colnames(data2collapsedDeconvoluted) %in% colnames(tmp))]
data2collapsedDeconvoluted=tmp
rm(tmp)

data1WGcollapsedDeconvoluted_c=data1WGcollapsedDeconvoluted_c[,!duplicated(data1WGcollapsedDeconvoluted_c$subjectId)]

xtmp=data1WGcollapsedDeconvoluted[,(data1WGcollapsedDeconvoluted$age>15&data1WGcollapsedDeconvoluted$diagnosis_binary=="TD")|(data1WGcollapsedDeconvoluted$age<40&data1WGcollapsedDeconvoluted$diagnosis_binary=="ASD")]
table(xtmp$diagnosis_binary,xtmp$batch)

t.test(xtmp$age~xtmp$diagnosis_binary)
fisher.test(table(xtmp$sex,xtmp$diagnosis_binary))
data1WGcollapsedDeconvoluted=xtmp
rm(xtmp)

print(paste0("number of samples in data2: ",ncol(data2collapsedDeconvoluted)))
print(paste0("number of samples in data1WG: ",ncol(data1WGcollapsedDeconvoluted)))

print(paste0("number of samples in dup. data1WG: ",ncol(data1WGcollapsedDeconvoluted_c)))


#assigning the subjects to the training and the test sets


sum(data1WGcollapsedDeconvoluted$subjectId %in% data2collapsedDeconvoluted$subjectId)


#making the training and the test sets

#Deconvoluted + collapsed

.myMainTrainingFn(data1WGcollapsedDeconvoluted,outPath="~/OneDrive - UC San Diego/Vahid/OneDrive - UC San Diego/classificationPaper/HT_as_test_05142019/ageBalanced_decov_genes_wF_v2/mainDataset_HT12/")
.myTestWriterFn(data1WGcollapsedDeconvoluted,data2collapsedDeconvoluted,outPath="~/OneDrive - UC San Diego/Vahid/OneDrive - UC San Diego/classificationPaper/HT_as_test_05142019/ageBalanced_decov_genes_wF_v2/testDataset_HT12/")
.myLDTestWriterFn(data1WGcollapsedDeconvoluted,data2collapsedDeconvoluted,outPath="~/OneDrive - UC San Diego/Vahid/OneDrive - UC San Diego/classificationPaper/HT_as_test_05142019/ageBalanced_decov_genes_wF_v2/LDtestDataset_HT12/")
.myTestWriterFn(data1WGcollapsedDeconvoluted,data1WGcollapsedDeconvoluted_c,outPath="~/OneDrive - UC San Diego/Vahid/OneDrive - UC San Diego/classificationPaper/HT_as_test_05142019/ageBalanced_decov_genes_wF_v2/Longitudinaltest/")



