# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)

split_dir <- '/LOCATION OF SPLIT FCS FILES'
export_dir <- '/LOCATION WHERE TO SAVE KNN IMPUTED FILES'

ids <- list.dirs(split_dir, recursive=FALSE, full.names=FALSE)
for (id in ids){
  print(id)
  # Read in the split and ground-truth flowframes
  FCS1 <- read.FCS(paste0(split_dir, id, '/', id, '_ff1.fcs')) 
  FCS2 <- read.FCS(paste0(split_dir, id, '/', id, '_ff2.fcs'))
  
  RelevantMarkers1 = seq(1, ncol(exprs(FCS1))-1)
  RelevantMarkers2 = seq(1, ncol(exprs(FCS2))-1)
  
  FCS1.data = as.data.frame(FCS1@exprs)[,RelevantMarkers1]
  VarNames1 = colnames(FCS1.data)
  
  FCS2.data = as.data.frame(FCS2@exprs)[,RelevantMarkers2]
  VarNames2 = colnames(FCS2.data)
  
  # Find shared Markers
  Matches = as.vector(matrix(0,ncol = length(VarNames1)))
  for (i in c(1:length(VarNames1))){
    if(any(VarNames1[i]==VarNames2))
      Matches[i] = which(VarNames1[i]==VarNames2)
  }
  
  # Reorder data
  Shared_Index = which(Matches>0)
  Data1.nonshared = FCS1.data[,-Shared_Index]
  VarNames1.nonshared = VarNames1[-Shared_Index]
  Data2.nonshared = FCS2.data[,-Matches[Shared_Index]]
  VarNames2.nonshared = VarNames2[-Matches[Shared_Index]]
  
  # This is the order data {Shared  Non_Shared}
  FCS1.data = cbind(FCS1.data[,Shared_Index], Data1.nonshared)
  VarNames1 = c(VarNames1[Shared_Index], VarNames1.nonshared)
  FCS2.data = cbind(FCS2.data[,Matches[Shared_Index]], Data2.nonshared)
  VarNames2 = c(VarNames2[Matches[Shared_Index]], VarNames2.nonshared)
  
  # Combine data
  m = length(Shared_Index)           # Number of shared markers
  
  IDX1 = FNN::get.knnx(FCS2.data[,1:m],FCS1.data[,1:m], k = 1, algorithm = "kd_tree")
  IDX1 = IDX1$nn.index
  IDX2 = FNN::get.knnx(FCS1.data[,1:m],FCS2.data[,1:m], k = 1, algorithm = "kd_tree")
  IDX2 = IDX2$nn.index
  
  Data.combine.1 = matrix(0,nrow = dim(FCS1.data)[1],ncol = dim(FCS1.data)[2]+dim(FCS2.data)[2]-m)
  Data.combine.2 = matrix(0,nrow = dim(FCS2.data)[1],ncol = dim(FCS1.data)[2]+dim(FCS2.data)[2]-m)
  
  for (i in c(1:dim(FCS1.data)[1])){
    Data.combine.1[i,1:m] = as.matrix(FCS1.data[i,1:m])
    Data.combine.1[i,(m+1):dim(FCS1.data)[2]] = as.matrix(FCS1.data[i,(m+1):dim(FCS1.data)[2]])
    Data.combine.1[i,(dim(FCS1.data)[2]+1):dim(Data.combine.1)[2]] = as.matrix(apply(FCS2.data[IDX1[i,],(m+1):dim(FCS2.data)[2]],2,median))
  }
  
  for (i in c(1:dim(FCS2.data)[1])){
    Data.combine.2[i,1:m] = as.matrix(FCS2.data[i,1:m])
    Data.combine.2[i,(m+1):dim(FCS1.data)[2]] = as.matrix(apply(FCS1.data[IDX2[i,],(m+1):dim(FCS1.data)[2]],2,median))
    Data.combine.2[i,(dim(FCS1.data)[2]+1):dim(Data.combine.2)[2]] = as.matrix(FCS2.data[i,(m+1):dim(FCS2.data)[2]])
  }
  
  Data.combine = rbind(Data.combine.1,Data.combine.2)
  VarNames.combine = c(VarNames1,VarNames2[(m+1):length(VarNames2)])
  
  Data.combine = as.data.frame(Data.combine)
  colnames(Data.combine) = VarNames.combine
  
  write.csv(Data.combine, paste0(export_dir, id, '_R_kNN.csv'))
}