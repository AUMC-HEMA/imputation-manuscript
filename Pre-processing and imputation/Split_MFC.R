# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(dplyr)

# Location of pre-gated files
fcs_dir <- '/PATH OF PRE-GATED FILES'
export_dir <- '/LOCATION WHERE TO SAVE FILES'

# Define different marker subsets
scatters <- c('FSC-A', 'FSC-H', 'FSC-W', 'SSC-A', 'SSC-H', 'SSC-W')
leukocytes <- c('HV500c-A')
tcells <- c('BUV395-A')
monocytes <- c('PerCP-Cy5-5-A')
tsub <- c('BUV737-A', 'BUV496-A')
tdiff <- c('BV421-A', 'APC-R700-A')
backbones <- list('Leukocytes_T-cells_Monocytes_T-sub_T-diff' = c(leukocytes, tcells, monocytes, tsub, tdiff))

# Split variable markers
# Pseudo-tube 1: CD57, KLRG1, PD1, CD27
markers1 <- c('FITC-A', 'APC-A', 'BV605-A', 'BV786-A')
# Pseudo-tube 2: CD28, CD95, TIM-3, TIGIT
markers2 <- c('PE-A', 'PE-CF594-A', 'BV711-A', 'PC7-A')
variable <- c(markers1, markers2)

file.names <- list.files(fcs_dir, pattern = "fcs$", full.names = F, recursive = T)
for (file in file.names){
  print(file)
  id <- strsplit(file, ' ')[[1]][5]
  id <- substr(id, 1, 7)
  print(id)
  # Read file and perform basic pre-processing
  ff <- read.FCS(paste0(dir, file))
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
                                                      c = 0)))
  }
  # Sample flowframe to equal number of cells
  # This prevents some issues with non-even dataframe manipulations later on
  n = nrow(exprs(ff))%/%1000*1000
  idx <- sample(seq(1, nrow(exprs(ff))), n, replace=FALSE)
  ff@exprs <- ff@exprs[idx,]
  
  for (backbone in names(backbones)){
    # Export real data
    exprs <- exprs(ff)
    
    # Where to save results
    export <- paste0(export_dir, id)
    dir.create(export)
    prefix <- paste0(export, '/' ,id)

    # Add column designating original ordering
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
  }
}