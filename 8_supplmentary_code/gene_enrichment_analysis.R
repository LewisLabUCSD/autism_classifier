#get the gene information
args = commandArgs(trailingOnly=TRUE)

rm(list=ls())
setwd("/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/6_sample scores/")
args="/Volumes/Work/Vahid_work/classification_bokan/new_runner/microarray_test_data/"

# Making test_data
#load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/testDataset_HT12/testDataset_HT12/out1/data.rda")
{
load("/Volumes/Work/Vahid_work/classification_newcode_data/HT_as_testv2_decon_genes_testDataset_HT12_out1_data.rda")
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
         file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_65_feature_only.rda")
### do the initialization for main+test
setwd('/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/')
source('pipelines.R')
require(Biobase,quietly = T)

load("/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_65_feature_only.rda")
output.directory='/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/'
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
############
# than go to the jabba and do 
# sudo chmod -R a+rwx /data/bokan
# nohup Rscript myFoldRunnerFn.R /data/bokan/classification_mainCode/new_runner/out3/ > /data/bokan/classification_mainCode/new_runner/out3/log.txt &
# 


### making the LD data
{
  load("/Volumes/Work/Vahid_work/classification_newcode_data/data_LD.rda")
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
       file="/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_55LD_feature_only.rda")
  ### do the initialization for main+test
  setwd('/Volumes/Work/Vahid_work/classification_newcode/autism_classifier/2_mainCode/')
  source('pipelines.R')
  require(Biobase,quietly = T)
  
  load("/Volumes/Work/Vahid_work/classification_newcode_data/inputExpData_170_with_55LD_feature_only.rda")
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
  # args="/Volumes/Work/Vahid_work/classification_newcode_data/test_runner/"
  path=gsub("\"", "", args[1])
  dataList=read.table(paste0(path,"dataList.txt"))
  path=unlist(strsplit(args,"/"))
  path=path[-length(path)]
  path=paste(path,collapse = "/")
  path=paste0(path,"/")
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
rm(list=ls())
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
# load("/Volumes/Work/Vahid_work/classification_newcode_data/final_testDataset_HT12_test_old/ResultsArranged.rda")

a2=as.data.frame(ROCres[["resMeans"]])
a2[is.na(a2)] <- 0
dim(a2)
a1=a1[match(rownames(a2),rownames(a1)),]
dim(a1)
a2=a2[match(rownames(a2),rownames(a1)),]
dim(a2)
# a2=a2[match(rownames(a1), rownames(a2)),]

# a3=as.vector(abs(as.matrix(a2)-as.matrix(a1)))
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

