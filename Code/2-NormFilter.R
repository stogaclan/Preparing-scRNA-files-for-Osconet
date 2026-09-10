rm(list=ls())
library(Oscope)
library(SingleCellExperiment)
customlib='C:/Users/SotoLabPC/Desktop/Osconet_stefano/OscoNet-master/R/CustomLibrary.R'
source(customlib)

####Setting Input Parameters

options(stringsAsFactors = FALSE)
OsconetInputdir="C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/OsconetInput"
meancutlow=1
alpha=0.2 #when alpha=alphacut=0.95 zero filtering does not filter out any further gene
case='her4_naive' #change this to something specific to your dataset or cell type


####################################################################
##                  Do not modify after this line                   ##
####################################################################

if (file.exists(OsconetInputdir)==FALSE){
  dir.create(OsconetInputdir)
}

# Read the TPM data and the corresponding metadata  containing the lables (cluster id) identified externally to this analysis by our collaborators
filename=paste0('C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/her4_naive_counts.csv')
metafile=paste0('C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/her4_naive_clusters.csv')
ALLDATA=read.csv(filename)
metadata=read.csv(metafile)
clusterId=sort(unique(metadata$x))
nc=length(clusterId)
ss=paste0('Number of Clusters in the input file ',nc)
print(ss)
for (i in c(1:nc)){
  label=clusterId[i]
  
  # #doing this because some clusters seemed to be exclusively NAs and was giving an issue
  if (label == 0 | label == 3 | label == 4 | label == 9) {
    next
  }
  
  
  ss=paste0('Normalizing Cluster ',label)
  print('###########')
  print(ss)
  D=Cluster_filtering(label,metadata,ALLDATA,alpha,meancutlow)
  fD=paste0('C:/Users/SotoLabPC/Desktop/Osconet_stefano/larval_test/OsconetInput/filter',case,'D',label,'_',alpha,'.csv')
  write.csv(D, file=fD,row.names = TRUE)
  
}





