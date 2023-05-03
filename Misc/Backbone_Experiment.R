# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(dplyr)
source('CyTOFmerge_MOD.R')

# Location of pre-gated files
fcs_dir <- 'PATH OF FCS FILES'
export_dir <- 'EXPORT PATH'

# Backbones in different orders depening on correlation
# with CD45RA
backbones <- list('TIM-3' = c('BV711-A'),
                  'TIM-3|TIGIT' = c('BV711-A', 'PC7-A'),
                  'TIM-3|TIGIT|CD27' = c('BV711-A', 'PC7-A', 'BV786-A'),
                  'TIM-3|TIGIT|CD27|CCR7' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1|KLRG1' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A', 'APC-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1|KLRG1|CD8' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A', 'APC-A', 'BUV496-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1|KLRG1|CD8|CD95' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A', 'APC-A', 'BUV496-A', 'PE-CF594-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1|KLRG1|CD8|CD95|CD4' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A', 'APC-A', 'BUV496-A', 'PE-CF594-A', 'BUV737-A'),
                  'TIM-3|TIGIT|CD27|CCR7|CD45|CD57|PD-1|KLRG1|CD8|CD95|CD4|CD28' = c('BV711-A', 'PC7-A', 'BV786-A', 'BV421-A', 'HV500c-A', 'FITC-A', 'BV605-A', 'APC-A', 'BUV496-A', 'PE-CF594-A', 'BUV737-A', 'PE-A'),
                  'CD28' = c('PE-A'),
                  'CD28|CD4' = c('PE-A', 'BUV737-A'),
                  'CD28|CD4|CD95' = c('PE-A', 'BUV737-A', 'PE-CF594-A'),
                  'CD28|CD4|CD95|CD8' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A'),
                  'CD28|CD4|CD95|CD8|KLRG1' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57|CD45' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A', 'HV500c-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57|CD45|CCR7' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A', 'HV500c-A', 'BV421-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57|CD45|CCR7|CD27' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A', 'HV500c-A', 'BV421-A', 'BV786-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57|CD45|CCR7|CD27|TIGIT' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A', 'HV500c-A', 'BV421-A', 'BV786-A', 'PC7-A'),
                  'CD28|CD4|CD95|CD8|KLRG1|PD-1|CD57|CD45|CCR7|CD27|TIGIT|TIM-3' = c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 'FITC-A', 'HV500c-A', 'BV421-A', 'BV786-A', 'PC7-A', 'BV711-A'))

# Incrementally build up 20 different backbones
full_backbone <- c('PE-A', 'BUV737-A', 'PE-CF594-A', 'BUV496-A', 'APC-A', 'BV605-A', 
                   'FITC-A', 'HV500c-A', 'BV421-A', 'BV786-A', 'PC7-A', 'BV711-A')
for (seed in 1:20){
  # Shuffle the markers
  set.seed(seed)
  # Sample until backbone size
  shuffled <- sample(full_backbone, 7)
  for (i in 1:7){
    markers <- shuffled[1:i]
    print(markers)
    name <- paste0('Seed', as.character(seed), '|', 'Sampled', as.character(i))
    backbones[[name]] <- markers
  }
}

# Variable markers
markers1 <- c('BUV395-A', 'PerCP-Cy5-5-A')
markers2 <- c('APC-R700-A')
variable <- c(markers1, markers2)

file.names <- list.files(fcs_dir, pattern = "fcs$", full.names = F, recursive = T)

for (file in file.names){
  print(file)
  id <- strsplit(file, ' ')[[1]][5]
  id <- substr(id, 1, 7)
  print(id)
  
  for (backbone in names(backbones)){
    # Read file and perform basic pre-processing
    ff <- read.FCS(paste0(fcs_dir, file))
    # Transform the data
    trans_list <- list('100'=c('BV605-A'),
                       '150'=c('FITC-A', 'PerCP-Cy5-5-A', 'APC-R700-A', 'APC-H7-A',
                               'BV421-A', 'HV500c-A', 'BV711-A', 'BV786-A', 'BUV395-A',
                               'BUV496-A', 'BUV737-A', 'PE-A', 'PE-CF594-A', 'PC7-A'),
                       '300'=c('APC-A'))
    for (cofactor in names(trans_list)){
      ff <- transform(ff,transformList(trans_list[[cofactor]],
                                       arcsinhTransform(a = 0,
                                                        b = 1/as.integer(cofactor),
                                                        c = 0)))}
    # Sample flowframe to equal number of cells
    # This prevents some issues with non-even dataframe manipulations later on
    n = nrow(exprs(ff))%/%1000*1000
    idx <- sample(seq(1, nrow(exprs(ff))), n, replace=FALSE)
    ff@exprs <- ff@exprs[idx,]
    
    # Export real data
    exprs <- exprs(ff)
    
    # Where to save results
    export <- paste0(export_dir, id)
    dir.create(export)
    prefix <- paste0(export, '/' , backbone, '_', id)
    print(prefix)
    
    # # Add column designating original ordering
    original_ID <- as.matrix(seq(1, nrow(exprs(ff))))
    colnames(original_ID) <- 'original_ID'
    ff <- fr_append_cols(ff, original_ID)
    
    # Split into marker sets, each containing half of the data
    sample_size <- nrow(exprs(ff)) / 2
    ff1_idx <- sample(seq(1, nrow(exprs(ff))), sample_size, replace=FALSE)
    ff1_idx <- sort(ff1_idx)
    ff2_idx <- setdiff(seq(1, nrow(exprs(ff))), ff1_idx)
    
    ff1 <- ff[ff1_idx,c(backbones[[backbone]], markers1, 'original_ID')]
    ff2 <- ff[ff2_idx,c(backbones[[backbone]], markers2, 'original_ID')]

    # Create a marker to keep track of imputed marker
    dataset <- as.matrix(seq(1, nrow(exprs(ff))))
    colnames(dataset) <- 'dataset'
    dataset[ff1_idx,] <- 1
    dataset[ff2_idx,] <- 2
    ff <- fr_append_cols(ff, dataset)

    # Export datasets
    write.FCS(ff, paste0(prefix, '_gt.fcs'))
    write.FCS(ff1, paste0(prefix, '_ff1.fcs'))
    write.FCS(ff2, paste0(prefix, '_ff2.fcs'))
    
    # Impute
    gt <- read.FCS(paste0(prefix, '_gt.fcs'))
    p <- pData(gt@parameters)
    
    # Impute using CyTOFmerge
    # Exclude original_ID columns using -1 in RelevantMarkers
    CyTOFmerge_exprs <- CombineFCS_MOD(FCSfile1=paste0(prefix, '_ff1.fcs'),
                                       RelevantMarkers1 = seq(1, ncol(exprs(ff1))-1),
                                       FCSfile2=paste0(prefix, '_ff2.fcs'), 
                                       RelevantMarkers2 = seq(1, ncol(exprs(ff2))-1), 
                                       arcsinhTrans = FALSE)
    CyTOFmerge_exprs <- as.matrix(CyTOFmerge_exprs)
    # Rename columns to original channel names
    channels <- colnames(CyTOFmerge_exprs)
    for (channel in channels){
      name <- subset(p, desc==channel)$name[[1]]
      colnames(CyTOFmerge_exprs)[which(colnames(CyTOFmerge_exprs) == channel)] <- name
    }
    # CyTOFmerge reshuffles the data so that ff1 is followed by ff2 instead of 
    # the original order in the ground-truth flowframe. We retrieve this original
    # ordering by using the original_ID variable.
    sample_size <- nrow(exprs(gt)) / 2
    imputed_ff1 <- CyTOFmerge_exprs[seq(1, sample_size),]
    imputed_ff2 <- CyTOFmerge_exprs[seq(sample_size+1, nrow(CyTOFmerge_exprs)),]
    # Retrieve original ordering of events
    ff1_idx <- ff1@exprs[,'original_ID']
    ff2_idx <- ff2@exprs[,'original_ID']
    ff_CyTOFmerge <- gt
    ff_CyTOFmerge@exprs[ff1_idx,c(variable)] <- imputed_ff1[,c(variable)]
    ff_CyTOFmerge@exprs[ff2_idx,c(variable)] <- imputed_ff2[,c(variable)]
    
    write.FCS(ff_CyTOFmerge, paste0(prefix, '_CyTOFmerge_exprs.fcs'))
  }
}
