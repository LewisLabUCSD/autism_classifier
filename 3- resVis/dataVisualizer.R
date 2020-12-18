library(ggplot2)
library(reshape2)
library(gplots)
library(WGCNA)
library(colorRamps)

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

#AU-ROC and AU-PR values that are zero are related to methods that didn't generate output, so we set them as 0.5
rm(list=ls())
#path="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/classificationDiscoverySet/"
#path="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/server_results/independent_dataSets/WG6/independentSetCorrect_complete/"
#path="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/WG6_scaled/WG6Scale/newWG6testScale/"
#path="~/OneDrive - UC San Diego/vahid_mac/Documents/data/classificationPaper/newData/WG6_kernel/WG6Kernel/newWG6test/"
load(paste0(path,"ResultsArranged.rda"))


perc85=.myResultOrganizerFn(perc85)
perc9=.myResultOrganizerFn(perc9)
perc95=.myResultOrganizerFn(perc95)
res85=.myResultOrganizerFn(res85)
res9=.myResultOrganizerFn(res9)
res95=.myResultOrganizerFn(res95)
ROCres=.myResultOrganizerFn(ROCres,replaceZero = T)
PRres=.myResultOrganizerFn(PRres,replaceZero = T)

tmp=which(res95>0.8&perc95>0.5,arr.ind = T)
tmp=data.frame(feature=row.names(res95)[tmp[,1]],classification=colnames(res95)[tmp[,2]],value=res95[tmp],perc=perc95[tmp])
tmp$name=paste0(tmp$feature,"_",tmp$classification)
tmpROC=ROCres[which(row.names(ROCres) %in% tmp$feature),]
tmpROC=tmpROC[,colnames(tmpROC) %in% unique(tmp$classification)]
x=which(tmpROC>0.8,arr.ind = T)
tmpROC=data.frame(feature=row.names(tmpROC)[x[,1]],classification=colnames(tmpROC)[x[,2]],roc=tmpROC[x])
tmpROC$name=paste0(tmpROC$feature,"_",tmpROC$classification)

tmp=merge(tmpROC,tmp,by="name")

#global heatmap for AUROC

ROCmat=ROCres
ROCmat=ROCmat[,-which(colnames(ROCmat)=="sqda")]
slpallete=c(blue2green(9),colorRampPalette(c("green", "yellow"))(9)[-1],colorRampPalette(c("yellow", "red"))(6)[-1])
hm=heatmap.2(as.matrix(ROCmat), margins = c(5,10),trace="none",density.info=c("none"),col=slpallete[max(round(min(ROCmat),1)*25,2):(round(max(ROCmat),2)*25)],labRow = F)
#finding high performing cluster
##set k in a way that separates the high performing cluster

k=5
n=cutree(as.hclust(hm$rowDendrogram),k)
nPalette=WGCNA:::labels2colors(seq(1:k))
#nPalette[4]="black"

nPalette=nPalette[n]

#nPalette[nPalette=="brown"]="orange"

heatmap.2(as.matrix(ROCmat), margins = c(5,10),trace="none",density.info=c("none"),col=slpallete[max(round(min(ROCmat),1)*25,2):(round(max(ROCmat),2)*25)],labRow = F,RowSideColors=nPalette)
#correspondence of colors with cluster numbers
data.frame(Number=1:k,color=WGCNA:::labels2colors(seq(1:k)))

#which color?
slColor=4
slClusters=ROCmat[n==slColor,]
pathList=row.names(slClusters)

save(pathList,file=paste0(path,"selectedPathList.rda"))


#detailed heatmaps of precision-recall and ROC

PRresArranged=PRres
PRresArranged=PRresArranged[,-which(colnames(PRresArranged)=="sqda")]
PRresArranged=melt(PRresArranged,id.vars='filtrationRoute')
p <- ggplot(PRresArranged, aes(variable,filtrationRoute)) 
p=p + geom_tile(aes(fill = value), colour = "white")+ scale_fill_gradient(name="colorKey",low = "black",high = "yellow",limits=c(PRsmRes[5], PRsmRes[6]))
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y = element_text(size=3),axis.text.x = element_text(size=4))+ggtitle("Precison-Recall results")


ROCsmRes=summary(as.vector(as.matrix(ROCresArranged)))
x=apply(ROCresArranged,1,function(x){length(which(x>ROCsmRes[5]))})
ROCresArranged=ROCresArranged[x>0,]

ROCresArranged=cbind(filtrationRoute=row.names(ROCresArranged),ROCresArranged)
ROCresArranged=melt(ROCresArranged,id.vars='filtrationRoute')
p <- ggplot(ROCresArranged, aes(variable,filtrationRoute)) 
p=p + geom_tile(aes(fill = value), colour = "white")+ scale_fill_gradient(name="colorKey",low = "black",high = "yellow",limits=c(ROCsmRes[5], ROCsmRes[6]))
p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.y = element_text(size=3),axis.text.x = element_text(size=4))+ggtitle("ROC results")

