#load result
load("/Volumes/Work/Vahid_work/classification_mainCode/classification_results.rda")
data_main_pdata = dataMain[["data"]][["concensus"]]
data_test_pdata = dataTest[["data"]][["concensus"]]
data_long_pdata = dataLong[["data"]][["concensus"]]
data_ld_pdata = dataLD[["data"]][["concensus"]]

sum(data_main_pdata$subjectId %in%  data_HT_pdata$subjectId )
sum(data_test_pdata$subjectId %in%  data_HT_pdata$subjectId )


# load 170 training data
load("/Volumes/Work/Vahid_work/classification_mainCode/inputExpData_170.rda")
test_subject = inputExpData@phenoData@data$subjectId
train_id = inputExpData@featureData@data[["ID"]]
duplicated(data_main_pdata$subjectId )
duplicated(test_subject)
data_main_pdata$subjectId %in% test_subject
sum(train_id %in% gene_id_test)

# load 65 training data

load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/testDataset_HT12/testDataset_HT12/out1/data.rda")
data_HT12_WG6_test_pdata = data[["testInputData"]]@phenoData@data
# check test are all(which is 64) in the 226 data set
sum(data_HT12_WG6_test_pdata[['subjectId']] %in% data_test_pdata$subjectId)
length(data_test_pdata$subjectId)
gene_id_test = data[["testInputData"]]@featureData@data$ID
# 
# data_ld_pdata$subjectId
# data_main_pdata$subjectId %in% data_test_pdata$subjectId
# load("/Volumes/Work/Vahid_work/classification_bokan/inputDataJabba.rda")
# data_226_pdata = inputExpData@phenoData@data
# 
# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/mainDataset_HT12_WG6/inputDataJabba.rda")
# data_HT12_WG6_pdata = inputExpData@phenoData@data
# # check training have 154 in the 226 data set
# sum(data_HT12_WG6_pdata$subjectId %in% data_226_pdata$subjectId)
# data_test_pdata$subjectId %in%  data_HT_pdata$subjectId 
# 
# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/nondeconvoluted_genes/mainDataset_HT12_WG6/inputDataJabba.rda")
# data_HT12_WG6_nd_gene_pdata = inputExpData@phenoData@data
# # check training have 154 in the 226 data set
# sum(data_HT12_WG6_nd_gene_pdata$subjectId %in% data_226_pdata$subjectId)
# 
# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/testDataset_HT12/testDataset_HT12/out1/data.rda")
# data_HT12_WG6_test_pdata = data[["testInputData"]]@phenoData@data
# # check test are all(which is 64) in the 226 data set
# sum(data_HT12_WG6_test_pdata[['subjectId']] %in% data_HT_pdata$subjectId)
# 
# data[["testInputData"]]@phenoData@data[["blooddrawid"]]
# 
# 
# data_HT_pdata = inputExpData@phenoData@data
# sum(data_main_pdata$subjectId %in% data_HT_pdata$subjectId)
# 
# sum(data_HT_pdata$subjectId %in% data_main_pdata$subjectId)
# 
# data_test_pdata$subjectId %in%  data_HT_pdata$subjectId 
# akk = inputExpData[,which(data_HT_pdata$subjectId %in% data_main_pdata$subjectId)]
# sum(data_HT_pdata$subjectId %in% data_main_pdata$subjectId)
# 
# data_main_pdata$subjectId %in% data_HT_pdata$subjectId 


#####################
load("/Volumes/Work/Vahid_work/classification_bokan/inputDataJabba_ht_as_testv2_perm1.rda")
probe_from_ht_as_testv2_perm1=inputExpData@featureData@data[["ID"]]
exprs_from_ht_as_testv2_perm1 = inputExpData@assayData[["exprs"]]
colnames(exprs_from_ht_as_testv2_perm1) = inputExpData@phenoData@data[["blooddrawid"]]
rownames(exprs_from_ht_as_testv2_perm1) = inputExpData@featureData@data[["ID"]]



# load("/Volumes/Work/Vahid_work/classification_bokan/inputDataJabba_wF_v2.rda")
# probe_from_wF_v2=inputExpData@featureData@data[["PROBE_ID"]]
# exprs_from_wF_v2 = inputExpData@assayData[["exprs"]]

head(probe_from_wF_v2)

head(probe_from_ht_as_testv2_perm1)
sum(probe_from_wF_v2 %in% probe_from_ht_as_testv2_perm1)

load("/Volumes/Work/Vahid_work/classification_bokan/inputDataJabba.rda")
probe_from_jabba=inputExpData@featureData@data[["ID"]]
exprs_jabba = inputExpData@assayData[["exprs"]]
colnames(exprs_jabba) = inputExpData@phenoData@data[["blooddrawid"]]
rownames(exprs_jabba) = inputExpData@featureData@data[["ID"]]


sum(probe_from_jabba %in% probe_from_ht_as_testv2_perm1)
sum(probe_from_jabba %in% union(probe_from_wF_v2,probe_from_ht_as_testv2_perm1))

# compare 2 samples

load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/nondeconvoluted_probes/mainDataset_HT12_WG6/inputDataJabba.rda")

probe_from_HT12_WG6=inputExpData@featureData@data[["ID"]]
exprs_from_HT12_WG6 = inputExpData@assayData[["exprs"]]
colnames(exprs_from_HT12_WG6) = inputExpData@phenoData@data[["blooddrawid"]]


sum(colnames(exprs_jabba) %in% colnames(exprs_from_HT12_WG6))

colnames(exprs_jabba)[1]
"576-01"
rownames(exprs_jabba) %in% rownames(exprs_from_HT12_WG6)
"13774"
rownames(exprs_jabba)[1:10]
"ILMN_1810810"
exprs_jabba[c("ILMN_1667796","ILMN_2331890","ILMN_1695261","ILMN_2094718","ILMN_1792528","ILMN_2100437" ), "576-01"]
exprs_from_HT12_WG6[c("ILMN_1667796","ILMN_2331890","ILMN_1695261","ILMN_2094718","ILMN_1792528","ILMN_2100437" ), "576-01"]
exprs_from_ht_as_testv2_perm1[c("ILMN_1667796","ILMN_2331890","ILMN_1695261","ILMN_2094718","ILMN_1792528","ILMN_2100437" ), "576-01"]
sum(probe_from_jabba %in% probe_from_HT12_WG6)

# sum(probe_from_jabba %in% union(probe_from_wF_v2,probe_from_ht_as_testv2_perm1))


# extract the ROC information

load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSet/out1/results.rda")
classification_1 = result[["PRresults"]]

# load("/Volumes/Work/Vahid_work/classification/classification/HT_as_test_v2/data/deconvoluted_genes/sampleScores/mainDataset_HT12_WG6/classificationSetInternalROC/out1/results.rda")
# classification_2 = result[["PRresults"]]



