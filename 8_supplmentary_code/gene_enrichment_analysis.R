#get the gene information

# LD
ld_name = rownames(data[["testInputData"]]@phenoData@data)[data[["testInputData"]]@phenoData@data[["diagnosis_binary"]]=='LD']
weightedE_table=dataLD[["data"]][["concensus"]]
weightedE_table[weightedE_table$sampleName %in% ld_name,]
hist(weightedE_table[weightedE_table$sampleName %in% ld_name,]$weightedEstimate)
args = commandArgs(trailingOnly=TRUE)

rm(list=ls())
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/6_sample scores/")
args="/Volumes/Work/Vahid_work/classification_bokan/new_runner/microarray_test_data/"

# Making test + main
{
load('/Volumes/Work/Vahid_work/classification_newcode_data/inputDataJabba_main.rda')
  
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/out1/data.rda")
# test_set=list()
testInputExpData=data[["testInputData"]]
testLabels=data[["testLabels"]]
testLabels=replace(testLabels, testLabels=='proband','ASD')
middleFns=middleFns
posteriorFns=posteriorFns
classificationMethodsList=classificationMethodsList
labels=as.vector(data[["features"]]@phenoData@data[["diagnosis_binary"]])
inputExpData=data[["features"]]
expClassName='ASD'
initializerFns=initializerFns

annotation <- fData(testInputExpData)
metaData <- data.frame(labelDescription=colnames(annotation)) 
probeinfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=metaData)
phenoData_ <- rbind(pData(testInputExpData), pData(inputExpData))
phenoData_ <- new("AnnotatedDataFrame", data=data.frame(phenoData_), varMetadata=data.frame(labelDescription=colnames(phenoData_)) )
expMat <- cbind(exprs(testInputExpData), exprs(inputExpData))
colnames(expMat)=rownames(phenoData_)
testInputExpData=new("ExpressionSet", exprs=expMat, featureData = probeinfo, phenoData=phenoData_)
testLabels = c(testLabels, labels)

save(inputExpData,labels,testInputExpData,testLabels,middleFns,initializerFns,posteriorFns,classificationMethodsList,expClassName,
     file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_main_with_240_tests.rda")

trainingInput=inputExpData
trainingLabels=labels
testInput=testInputExpData
source('pipelines.R')

results=.myInitializer(trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,method=initializerFns,prevMethod="")

counter=1
output.directory="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_test_main_added/"
outputfilesList=""
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


# Making ASD+TD+LD 
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/out1/data.rda")

# test_set=list()
testInputExpData=data[["testInputData"]]
testLabels=data[["testLabels"]]
middleFns=middleFns
posteriorFns=posteriorFns
classificationMethodsList=classificationMethodsList
# load('/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/out1/data.rda')
load("/Volumes/Work/Vahid_work/classification_newcode_data/inputDataJabba_main.rda")

testInputExpData_asd=data[["testInputData"]][, which(data[["testLabels"]]=='proband')]

labels=replace(labels, labels=='proband','ASD')

testLabels=c(testLabels, data[["testLabels"]][which(data[["testLabels"]]=='proband')])
testLabels=replace(testLabels, testLabels=='proband','ASD')

expClassName='ASD'
initializerFns=initializerFns


annotation <- fData(testInputExpData)
metaData <- data.frame(labelDescription=colnames(annotation)) 
probeinfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=metaData)
phenoData_ <- rbind(pData(testInputExpData), pData(testInputExpData_asd))
phenoData_ <- new("AnnotatedDataFrame", data=data.frame(phenoData_), varMetadata=data.frame(labelDescription=colnames(phenoData_)) )
expMat <- cbind(exprs(testInputExpData), exprs(testInputExpData_asd))
colnames(expMat) = rownames(phenoData_)
testInputExpData=new("ExpressionSet", exprs=expMat, featureData = probeinfo, phenoData=phenoData_)


save(inputExpData,labels,testInputExpData,testLabels,middleFns,initializerFns,posteriorFns,classificationMethodsList,expClassName,
     file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_main_with_ASD_and_LD.rda")


# Making ASD+LD 
rm(list=ls())
load("/Volumes/Work/Vahid_work/classification_newcode_data/inputDataJabba_main.rda")
rm(classificationMethodsList, posteriorFns, middleFns, expClassName, runReal, npermTest, nfold, ncores)
labels=replace(labels, labels=='proband','ASD')

load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/out1/data.rda")

# test_set=list()
testInputExpData=data[["testInputData"]][, which(data[["testLabels"]]=='proband')]
testLabels=data[["testLabels"]][which(data[["testLabels"]]=='proband')]
testLabels=replace(testLabels, testLabels=='proband','ASD')

# load('/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/out1/data.rda')
load("/Volumes/Work/Vahid_work/classification_newcode_data/data_LD.rda")
testInputExpData_ld=data[["testInputData"]][, which(data[["testLabels"]]=='proband')]

testLabels=c(testLabels, data[["testLabels"]][which(data[["testLabels"]]=='proband')])
testLabels=replace(testLabels, testLabels=='proband','TD')

expClassName='ASD'

annotation <- fData(testInputExpData)
metaData <- data.frame(labelDescription=colnames(annotation)) 
probeinfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=metaData)
phenoData_ <- rbind(pData(testInputExpData), pData(testInputExpData_ld))
phenoData_ <- new("AnnotatedDataFrame", data=data.frame(phenoData_), varMetadata=data.frame(labelDescription=colnames(phenoData_)) )
expMat <- cbind(exprs(testInputExpData), exprs(testInputExpData_ld))
colnames(expMat) = rownames(phenoData_)
testInputExpData=new("ExpressionSet", exprs=expMat, featureData = probeinfo, phenoData=phenoData_)


save(inputExpData,labels,testInputExpData,testLabels,middleFns,initializerFns,posteriorFns,classificationMethodsList,expClassName,
     file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_main_with_ASD_and_LD.rda")

############
# than go to the jabba and do 
# sudo chmod -R a+rwx /data/bokan
# nohup Rscript myFoldRunnerFn.R /data/bokan/classification_mainCode/new_runner/out3/ > /data/bokan/classification_mainCode/new_runner/out3/log.txt &
# 


### making the LD data
{
  
  ### do the initialization for main+test
  setwd('/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/')
  source('pipelines.R')
  require(Biobase,quietly = T)
  
  load("/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_main_with_ASD_and_LD.rda")
  output.directory='/Volumes/Work/Vahid_work/classification_newcode_data/ld_runner/'
  if(!file.exists(output.directory)){
    dir.create(output.directory)
  } else {
    unlink(output.directory, recursive = T, force = T)
    dir.create(output.directory)
  }
  
  
  trainingInput=inputExpData
  trainingLabels=labels
  testInput=testInputExpData
  method=initializerFns
  testLabels=testLabels
  testInputData=testInput
  results=.myInitializer(trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,method=initializerFns,prevMethod="")
  
  counter=1
  outputfilesList=""
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

############################
# processing the ld_result #
############################
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/3_resVis/")
source('myDataCollector.R')
args="/Volumes/Work/Vahid_work/classification_newcode_data/ld_runner/"
# args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSetInternalROC/"
.myFeatureCollector(args)
.myDataCollector(args)



### making the longitude data
{
  load("/Volumes/Work/Vahid_work/classification_newcode_data/data_Longitudinal.rda")
  # test_set=list()
  testInputExpData=data[["testInputData"]]
  testLabels=data[["testLabels"]]
  testLabels=replace(testLabels, testLabels=='proband','ASD')
  middleFns=middleFns
  posteriorFns=posteriorFns
  classificationMethodsList=classificationMethodsList
  
  ###
  load('/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170.rda')
  labels=as.vector(inputExpData@phenoData@data[["diagnosis_binary"]])
  inputExpData=inputExpData
  expClassName='ASD'
  initializerFns=initializerFns
  
  save(inputExpData,labels,testInputExpData,testLabels,middleFns,initializerFns,posteriorFns,classificationMethodsList,expClassName,
       file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_33longitudinal_feature_only.rda")
  ### do the initialization for main+test
  setwd('/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/')
  source('pipelines.R')
  require(Biobase,quietly = T)
  
  load("/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_33longitudinal_feature_only.rda")
  output.directory='/Volumes/Work/Vahid_work/classification_newcode_data/longitudinal_runner/'
  if(!file.exists(output.directory)){
    dir.create(output.directory)
  } else {
    unlink(output.directory, recursive = T, force = T)
    dir.create(output.directory)
  }
  
  
  trainingInput=inputExpData
  trainingLabels=labels
  testInput=testInputExpData
  results=.myInitializer(trainingInput,labels=trainingLabels,testInputData=testInput,testLabels=testLabels,method=initializerFns,prevMethod="")
  
  counter=1
  outputfilesList=""
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



########################################
#### get the features out from main ####

ommandArgs(trailingOnly=TRUE)

rm(list=ls())
# args="/Volumes/Work/Vahid_work/classification_newcode_data/main_runner/"
.myFeatureCollector=function(args){
  # args="/Volumes/Work/Vahid_work/classification_newcode_data/ld_runner/"
  path=gsub("\"", "", args[1])
  dataList=read.table(paste0(path,"dataList.txt"))
  if (length(strsplit(as.character(dataList$V1[1]),"/")[[1]])>2){
    path=''
  }else{
  path=unlist(strsplit(args,"/"))
  path=path[-length(path)]
  path=paste(path,collapse = "/")
  path=paste0(path,"/")}
  
  # path
  res=list()
  error_count=0
  head_count = 0
  total_count=0
  for(i in 1:nrow(dataList)){
    if(file.exists(paste0(path,dataList$V1[i],'resultsPosterior.rda'))){
      load(paste0(path,dataList$V1[i],'resultsPosterior.rda'))
      for(j in 1:length(resultsPosterior)){
        total_count=total_count+1
        if(nrow(resultsPosterior[[j]]$features)==1){
          if(ncol(resultsPosterior[[j]]$features)!=length(resultsPosterior[[j]]$labels)){
            print("error single")
            
          }
          else{error_count=error_count+1}
        }
        else{
          if(nrow(resultsPosterior[[j]]$features)!=length(resultsPosterior[[j]]$labels)){
            print("error")
            print(paste(as.character(i),as.character(j),sep=' '))
            print(nrow(resultsPosterior[[j]]$features))
            print(length(resultsPosterior[[j]]$labels))
            # error_count=error_count+1
          }
          head_count=head_count+1
          res=c(res,list(list(features=colnames(resultsPosterior[[j]]$features),method=resultsPosterior[[j]]$method)))
        }
      }
    }
  }
  print(error_count)
  print(head_count)
  print(total_count)
  save(res,file=paste0(args[1],"feature.rda"))
  
}
.myFeatureCollector(args)

# load feature.rda
load("/Volumes/Work/Vahid_work/classification_newcode_data/main_runner_new/feature.rda")
library(jsonlite)

exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/main_runner_new/feature.json")
list2 <- fromJSON("test.json")
#
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/3_resVis/")
source('myDataCollector.R')
args="/Volumes/Work/Vahid_work/classification_newcode_data/main_runner/"
.myDataCollector(args)

load('/Volumes/Work/Vahid_work/classification_newcode_data/main_runner_new/ResultsArranged.rda')
# name_list = c()
# for(i in 1:length(res)){
#   name_list = c(name_list, res[[i]]$method)
# }
# library(data.table) # install if not installed already
# fwrite(list(name_list), file = "/Volumes/Work/Vahid_work/classification_newcode_data/main_runner/myFile.csv")
# # load feature.rda for main+test data set

write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/main_runner_new/roc_mean.cvs",sep = ",")


args="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/"
args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSetInternalROC/"
.myFeatureCollector(args)
load("/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/feature.rda")
.myFeatureCollector()
#
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/3_resVis/")
source('myDataCollector.R')
args="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/"
.myDataCollector(args)

load('/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/ResultsArranged.rda')
write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/roc_mean.csv",sep = ",")
write.table(ROCres[["resConts"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/roc_counts.csv",sep = ",")

load("/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/feature.json")
# list2 <- fromJSON("test.json")
##############

# resultsPosterior=resultsPosterior
# setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/")
# source('myFoldRunnerFn.R')
# result=.myClassificationBenchmarking(resultsPosterior,ncores=1,expClassName=expClassName,methodsList=c('sqda'))
# 
# 
# i=1
# inputDataList=resultsPosterior
# methodsList=c('sqda')
# experimentClassName=expClassName='ASD'
# trainingInputData=trainFeatures=inputDataList[[i]][["features"]]
# trainingLabels=trainLabels=inputDataList[[i]][["labels"]]
# testInputData=tstFeatures=inputDataList[[i]][["testInputData"]]
# testLabels=tstLabels=inputDataList[[i]][["testLabels"]]

# Redo the LD and long.
# Use the current data 

# Currently just the Gene information.
# But you also need to check the consistency. 
# get high AU-ROC model information from So you check the consistency with the AU-ROC first 
# and then go directly to the gene information. You need gene from the high AU-ROC score method.
# 
# So get the overlap of the two model first.

# library(diffr)
# diffr("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/ClassificationModule.R", 
      # "/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/6_sample scores/ClassificationModule.R")
###################
# for the main ####
###################
rm(list=ls())
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/3_resVis/")
source('myDataCollector.R')


args="/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/"
.myDataCollector(args)
.myFeatureCollector(args)



args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSet/"
# args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/mainDataset_HT12_WG6/classificationSet/"
.myDataCollector(args)
.myFeatureCollector(args)

load('/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/ResultsArranged.rda')
a1=as.data.frame(ROCres[["resMeans"]])
rm(ROCres)
a1[is.na(a1)] <- 0
dim(a1)
# rm(ROCres)
load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSet/ResultsArranged.rda")
# load('/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main_v1/ResultsArranged.rda')
a2=as.data.frame(ROCres[["resMeans"]])
a2[is.na(a2)] <- 0
dim(a2)
a1=a1[match(rownames(a2),rownames(a1)),]
dim(a1)
a2=a2[match(rownames(a2),rownames(a1)),]
dim(a2)
a3 = as.matrix(abs(a2-a1))
View(a3)

rm(list=ls())
load('/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/ResultsArranged.rda')
write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/roc_mean.csv",sep = ",")
write.table(ROCres[["resConts"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/roc_counts.csv",sep = ",")

load("/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/final_classificationSet_main/feature.json")


#################
# check the test dataset
# ###################
# for the test ####
###################
{rm(list=ls())
args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/"
.myDataCollector(args)
.myFeatureCollector(args)

args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/testDataset_HT12/testDataset_HT12/"
.myDataCollector(args)
# .myFeatureCollector(args)
#
rm(list=ls())
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/ResultsArranged.rda")
a1=as.data.frame(ROCres[["resMeans"]])
rm(ROCres)
a1[is.na(a1)] <- 0
dim(a1)
load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/testDataset_HT12/testDataset_HT12/ResultsArranged.rda")
# 
a3 = as.matrix(abs(a2-a1))
View(a3)
rm(list=ls())
load('/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/ResultsArranged.rda')
write.table(ROCres[["resConts"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/roc_counts.csv",sep = ",")
write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/roc_mean.csv",sep = ",")
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test/feature.json")
# list2 <- fromJSON("test.json")
}

{# ###################
  # for the test with main####
  ###################
rm(list=ls())
args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_test_main_added/"
.myDataCollector(args)
.myFeatureCollector(args)

load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_test_main_added/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_test_main_added/feature.json")

# args="/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/testDataset_HT12/testDataset_HT12/"
# .myDataCollector(args)
}

# args="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/"

###################
# for the LD ####
###################

rm(list=ls())
args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/"
.myDataCollector(args)
.myFeatureCollector(args)

args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD_old/"
.myDataCollector(args)
# .myFeatureCollector(args)
#
rm(list=ls())
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/ResultsArranged.rda")
a1=as.data.frame(ROCres[["resMeans"]])
rm(ROCres)
a1[is.na(a1)] <- 0
dim(a1)
# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/testDataset_HT12/testDataset_HT12/ResultsArranged.rda")
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD_old/ResultsArranged.rda")

a2=as.data.frame(ROCres[["resMeans"]])
a2[is.na(a2)] <- 0
dim(a2)
a1=a1[match(rownames(a2),rownames(a1)),]
dim(a1)
a2=a2[match(rownames(a2),rownames(a1)),]
dim(a2)
# a3=as.vector(abs(as.matrix(a2)-as.matrix(a1)))

# 
a3 = as.matrix(abs(a2-a1))
View(a3)
rm(list=ls())
load('/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/ResultsArranged.rda')
write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/roc_mean.csv",sep = ",")
write.table(ROCres[["resConts"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/roc_counts.csv",sep = ",")

load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_LD/feature.json")

###################
# for the longitudinal ####
###################

rm(list=ls())
args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/"
.myDataCollector(args)
.myFeatureCollector(args)

args="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal_old/"
.myDataCollector(args)
# .myFeatureCollector(args)
#
rm(list=ls())
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/ResultsArranged.rda")
a1=as.data.frame(ROCres[["resMeans"]])
rm(ROCres)
a1[is.na(a1)] <- 0
dim(a1)
# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/testDataset_HT12/testDataset_HT12/ResultsArranged.rda")
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal_old/ResultsArranged.rda")

a2=as.data.frame(ROCres[["resMeans"]])
a2[is.na(a2)] <- 0
dim(a2)
a1=a1[match(rownames(a2),rownames(a1)),]
dim(a1)
a2=a2[match(rownames(a2),rownames(a1)),]
dim(a2)
# a3=as.vector(abs(as.matrix(a2)-as.matrix(a1)))

# 
a3 = a1-a2
View(a3)
rm(list=ls())
load('/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/ResultsArranged.rda')
write.table(ROCres[["resConts"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/roc_counts.csv",sep = ",")
write.table(ROCres[["resMeans"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/roc_mean.csv",sep = ",")
load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/feature.rda")
exportJSON <- toJSON(res)
write(exportJSON, "/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_longitudinal/feature.json")


load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/mainDataset_HT12_WG6/inputDataJabba_main.rda")
View(inputExpData@featureData@data)
write.table(inputExpData@featureData@data, file="/Volumes/Work/Vahid_work/classification_newcode_data/final_result_plot/gene_info.csv",sep = "\t")


load("/Volumes/Work/Vahid_work/classification_newcode_data/test_hold/data.rda")
middleFns =list()
middleFns[[1]]='no'
middleFns
posteriorFns=list()
posteriorFns[[1]]='pcr'
posteriorFns[[2]]='plsr'
posteriorFns[[3]]='cppls'
save(classificationMethodsList,data,middleFns,posteriorFns, file = "/Volumes/Work/Vahid_work/classification_newcode_data/test_hold/data.rda")
save(resultsMiddle,posteriorFns, file = "/Volumes/Work/Vahid_work/classification_newcode_data/test_hold/resultsMiddle.rda")



#######################3
inputExpSet=resultsMiddle[[1]]$features
a=as.numeric(rownames(inputExpSet))
inputExpSet <- inputExpSet[order(a),]
dim(inputExpSet)
labels=resultsMiddle[[1]]$labels
testInputData=resultsMiddle[[1]]$testInputData
row.names(testInputData) = as.numeric(rownames(testInputData))
testInputData <- testInputData[ order(row.names(testInputData)),]
testLabels=resultsMiddle[[1]]$testLabels

if(is.null(outputsize)){
  outputsize=0.9
}
if(outputsize>0.99){
  print("output size is supposed to be less than 1")
  print("considering output size of 0.9")
  outputsize=0.9
}
pcaObj=prcomp(t(.matrixExtraction(inputExpSet)),scale. = T)
varExplained=cumsum(pcaObj$sdev^2 / sum(pcaObj$sdev^2))
ncomp=min(which(varExplained>=outputsize))

dataPrepared=list(as.numeric(as.factor(labels)),t(.matrixExtraction(inputExpSet)))
class(dataPrepared)="data.frame"
names(dataPrepared)=c("labels","covariates")
row.names(dataPrepared)=colnames(inputExpSet)
#print(paste0("starting with ",max(ncomp,5)," components"))

#pls.options(parallel = makeCluster(4, type = "PSOCK"))
res_pcr <- pcr(labels ~ covariates, ncomp = 92, data = dataPrepared, validation = "CV",segments=5,scale=TRUE)
explvar(res_pcr)
aload = unclass(res_pcr[["loadings"]])
aload = abs(aload)
aload = sweep(aload, 2, colSums(aload), "/")
aload = sort(aload[,1], decreasing = TRUE)
aload_cum = cumsum(aload)
minweight = min(which(aload_cum>=0.66))
gene_name = names(aload)[1:minweight]
hist(aload[,1])
plot(aload_cum)

res_plsr <- plsr(labels ~ covariates, ncomp = 3, data = dataPrepared, validation = "CV", scale=TRUE)
explvar(res_plsr)
plot(res_plsr, "weights", comps = 1:2, legendpos = "topleft",labels = "numbers", xlab = "nm")
aweights = unclass(res_plsr[["loading.weights"]])
aweights = abs(aweights)
aweights = sweep(aweights, 2, colSums(aweights), "/")
aweights = as.matrix(aweights) %*% as.matrix(explvar(res_plsr))
akk = sort(aweights, decreasing = TRUE)
weightExplained = cumsum(akk)
the_last = min(which(weightExplained>=weightExplained[length(weightExplained)]*0.66))
write(names(akk)[1:the_last], "/Volumes/Work/Vahid_work/classification_newcode_data/name.txt")
print(min(which(weightExplained>=weightExplained[length(weightExplained)]*0.66)))

akk = sort(apply(aweights,1, sum), decreasing = TRUE)
weightExplained = cumsum(akk)
minweight = min(which(weightExplained>=1))
minweight

res_cppls <- cppls(labels ~ covariates, ncomp = 5, data = dataPrepared, validation = "CV",segments=5,scale=TRUE)
ncomp=selectNcomp(res_cppls, method = "randomization", plot = F)
ncomp=max(ncomp,3)
#print(paste0("CPPLS number of components by CV:", ncomp))
res_cppls <- cppls(labels ~ covariates, ncomp = ncomp, data = dataPrepared, validation = "CV",scale=TRUE)
ncomp
aweights = unclass(res_cppls[["loading.weights"]])
aweights = abs(aweights)
aweights = sweep(aweights, 2, colSums(aweights), "/")
View()

# View(sort(aweights[,1], decreasing = TRUE))


# testInputData=predict(res,newdata = t(.matrixExtraction(testInputData)))
#print(paste0("CPPLS number of components by CV:", ncomp))
res_cppls <- cppls(labels ~ covariates, ncomp = ncomp, data = dataPrepared, validation = "CV",scale=TRUE)
ncomp
# colSums(sweep(aload, 2, colSums(aload), "/"))
aweights = unclass(res_cppls[["loading.weights"]])
aweights = abs(aweights)
aweights = sweep(aweight0s, 2, colSums(aweights), "/")
View(sort(apply(aweights,1, sum), decreasing = TRUE))
View(sort(aweights[,1], decreasing = TRUE))
head(aweights[,1:3])
# colSums(sweep(aweights, 2, colSums(aweights), "/"))
sum(aload[,1])

# 
# library(pls)
# data(gasoline)
# gasTrain <- gasoline[1:50,]
# gasTest <- gasoline[51:60,]
# gas1 <- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")
# summary(gas1)
# plot(gas1, ncomp = 2, asp = 1, line = TRUE)
# plot(gas1, plottype = "scores", comps = 1:3)
# plot(gas1, "loadings", comps = 1:2, legendpos = "topleft",labels = "numbers", xlab = "nm")

write.table(inputExpData@featureData@data[["Entrez_Gene_ID"]], file="/Volumes/Work/Vahid_work/classification_newcode_data/final_result_plot/ALL_GENE_FEATURE.txt",row.names = F,col.names = F)
