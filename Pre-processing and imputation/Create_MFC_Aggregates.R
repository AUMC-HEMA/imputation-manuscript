# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(FlowSOM)

dir <- '/PATH WHERE IMPUTED DATA IS SAVED'

for (method in c('CyTOFmerge', 'CytoBackBone', 'cyCombine', 'Infinicyt')){
  file.paths <- list.files(dir, pattern=paste0('_', method, '.fcs'), full.names = T, recursive = T)
  ff <- AggregateFlowFrames(file.paths, cTotal=999999999999999999, keepOrder=TRUE, silent=FALSE)
  ff@description$FILENAME <- paste0('Flow_', method, '_agg.fcs')
  ff@description$FIL <- paste0('Flow_', method, '_agg.fcs')
  ff@description$`$FIL` <- paste0('Flow_', method, '_agg.fcs')
  write.FCS(ff, paste0('Flow_', method, '_agg.fcs'))
}