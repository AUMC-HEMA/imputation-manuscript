# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(CytoBackBone)
library(FlowSOM)
library(cyCombine)
library(dplyr)
devtools::source_url("https://github.com/tabdelaal/CyTOFmerge/blob/master/CombineFCS.R?raw=TRUE")
source('CytoBackBone_MOD.R')

# Location of split files
split_dir <- '/PATH OF PSEUDO-TUBE FILES'
export_dir <- '/PATH WHERE TO SAVE RESULTS'
ids <- list.dirs(split_dir, recursive=FALSE, full.names=FALSE)

# Define different marker subsets
leukocytes <- c('HV500c-A')
tcells <- c('BUV395-A')
monocytes <- c('PerCP-Cy5-5-A')
tsub <- c('BUV737-A', 'BUV496-A')
tdiff <- c('BV421-A', 'APC-R700-A')
backbone <- c(leukocytes, tcells, monocytes, tsub, tdiff)
# Pseudo-tube 1: CD57, KLRG1, PD1, CD27
markers1 <- c('FITC-A', 'APC-A', 'BV605-A', 'BV786-A')
# Pseudo-tube 2: CD28, CD95, TIM-3, TIGIT
markers2 <- c('PE-A', 'PE-CF594-A', 'BV711-A', 'PC7-A')
variable <- c(markers1, markers2)

for (id in ids){
  print(id)
  # Read in the split and ground-truth flowframes
  ff1 <- read.FCS(paste0(split_dir, id, '/', id, '_ff1.fcs')) 
  ff2 <- read.FCS(paste0(split_dir, id, '/', id, '_ff2.fcs'))
  gt <- read.FCS(paste0(split_dir, id, '/', id, '_gt.fcs'))

  prefix <- paste0(export_dir, id, '/')
  dir.create(prefix)

  ##############################################################################
  # Impute using CytoBackBone
  # Requires FCS object instead of FlowFrame
  write.FCS(ff1, 'ff1.fcs')
  ff1_FCS <- import.FCS('ff1.fcs', trans='none')
  write.FCS(ff2, 'ff2.fcs')
  ff2_FCS <- import.FCS('ff2.fcs', trans='none')
  
  # Add (missing) scatter descriptions
  p <- pData(gt@parameters)
  p$desc[c(1:6)] <- c('FSC.A', 'FSC.H', 'FSC.W', 'SSC.A', 'SSC.H', 'SSC.W')
  backbone_names <- p[p$name %in% backbone,]$desc
  
  # Set CytoBackBone distance threshold for NN algorithm
  th <- length(backbone_names)*0.2
  
  # Prevent using original IDs
  ff1_FCS <- ff1_FCS[,head(ff1_FCS@markers, -1)]
  ff2_FCS <- ff2_FCS[,head(ff2_FCS@markers, -1)]
  
  # Run modified version of algorithm which keeps track of iterations
  merged <- merge_mod(ff1_FCS, ff2_FCS, BBmarkers=backbone_names, normalize=FALSE, th=th)
  
  # Export iteration data for cell-cell matching in Python
  write.csv(merged$metadata, paste0(prefix, id, '_CBB_meta.csv'))
  write.csv(merged$FCS1_excluded, paste0(prefix, id, '_FCS1_exclusions.csv'))
  write.csv(merged$FCS2_excluded, paste0(prefix, id, '_FCS2_exclusions.csv'))
  
  # Convert merged FCS into FlowFrame
  export(merged$FCS, paste0(prefix, id, '_merged.fcs'), transform=FALSE)
  ff_CBB <- read.FCS(paste0(prefix, id, '_merged.fcs'))
  # Counter-act some possible artefacts from exporting/reading FCS
  merged_exprs <- merged$FCS@intensities
  colnames(merged_exprs) <- merged$FCS@markers
  ff_CBB@exprs <- merged_exprs
  
  # Re-insert parameter naming from ground-truth flowframe
  new_p <- pData(ff_CBB@parameters)
  ff_CBB_exprs <- exprs(ff_CBB)
  channels <- colnames(exprs(ff_CBB))
  for (channel in channels){
    name <- subset(p, desc==channel)$name[[1]]
    new_p[new_p$desc==channel,]$name <- name
    colnames(ff_CBB_exprs)[which(colnames(ff_CBB_exprs) == channel)] <- name
  }
  ff_CBB@parameters@data <- new_p
  exprs(ff_CBB) <- ff_CBB_exprs
  
  ##############################################################################
  # Impute using CyTOFmerge
  # Exclude original_ID columns using -1 in RelevantMarkers
  CyTOFmerge_exprs <- CombineFCS(FCSfile1='ff1.fcs', 
                                 RelevantMarkers1 = seq(1, ncol(exprs(ff1))-1),
                                 FCSfile2='ff2.fcs', 
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
  
  ############################################################################
  # Impute using cyCombine
  dataset1 <- data.frame(exprs(ff1))
  dataset2 <- data.frame(exprs(ff2))
  names1 <- head(colnames(dataset1), -1)
  names2 <- head(colnames(dataset2), -1)
  overlap <- intersect(names1, names2)
  missing_1 <- names1[!(names1 %in% overlap)]
  missing_2 <- names2[!(names2 %in% overlap)]
  
  cyCombine <- impute_across_panels(dataset1 = dataset1, 
                                    dataset2 = dataset2,
                                    overlap_channels = overlap,
                                    impute_channels1 = missing_2,
                                    impute_channels2 = missing_1,
                                    xdim=15,
                                    ydim=15)
  cyCombine_exprs <- as.matrix(rbind(cyCombine$dataset1, cyCombine$dataset2))
  colnames(cyCombine_exprs) <- gsub(".", "-", colnames(cyCombine_exprs), fixed=TRUE)
  
  # Insert imputed data in the original flowframe
  imputed_ff1 <- cyCombine_exprs[seq(1, sample_size),]
  imputed_ff2 <- cyCombine_exprs[seq(sample_size+1, nrow(cyCombine_exprs)),]
  ff_cyCombine <- gt
  ff_cyCombine@exprs[ff1_idx,c(variable)] <- imputed_ff1[,c(variable)]
  ff_cyCombine@exprs[ff2_idx,c(variable)] <- imputed_ff2[,c(variable)]
  
  # Mark values with NAs, we add those back after clustering
  cyCombine_NAs <- ff_cyCombine@exprs[!complete.cases(ff_cyCombine@exprs), ]
  ff_cyCombine@exprs <- ff_cyCombine@exprs[complete.cases(ff_cyCombine@exprs), ]
  
  ##############################################################################
  # Import Infinicyt imputed values
  infinicyt_exprs <- read.csv(paste0(split_dir, id, '/', id, '_infinicyt.txt'), sep=';')
  # Rename columns
  colnames(infinicyt_exprs)[colnames(infinicyt_exprs) == "PD.1.BV605.A"] <- 'PD-1.BV605-A'
  colnames(infinicyt_exprs)[colnames(infinicyt_exprs) == "TIM.3.BV711.A"] <- 'TIM-3.BV711-A'
  colnames(infinicyt_exprs) <- gsub("^[^.]*.", "", colnames(infinicyt_exprs))
  colnames(infinicyt_exprs) <- gsub(".", "-", colnames(infinicyt_exprs), fixed=TRUE)
  infinicyt_exprs <- as.matrix(infinicyt_exprs)
  # Recreate flowframe with correct ordering
  imputed_ff1 <- infinicyt_exprs[seq(1, sample_size),]
  imputed_ff2 <- infinicyt_exprs[seq(sample_size+1, nrow(infinicyt_exprs)),]
  ff_infinicyt <- gt
  ff_infinicyt@exprs[ff1_idx,c(variable)] <- imputed_ff1[,c(variable)]
  ff_infinicyt@exprs[ff2_idx,c(variable)] <- imputed_ff2[,c(variable)]
  
  ############################################################################
  imputed_ffs <- list('CytoBackBone'=ff_CBB,
                      'CyTOFmerge'=ff_CyTOFmerge,
                      'cyCombine'=ff_cyCombine,
                      'Infinicyt'=ff_infinicyt)
  
  for (method in names(imputed_ffs)){
    # Concatenate ground-truth expression to imputed data in a new flowframe
    combined <- gt
    # Mark ground-truth data with a 0
    imp_state <- as.matrix(rep(0, nrow(exprs(combined))))
    colnames(imp_state) <- 'imp_state'
    combined <- fr_append_cols(combined, imp_state)
    
    # Mark imputed data with a 1
    imp_state <- as.matrix(rep(1, nrow(exprs(imputed_ffs[[method]]))))
    colnames(imp_state) <- 'imp_state'
    imputed_ffs[[method]] <- fr_append_cols(imputed_ffs[[method]], imp_state)
    
    # Concatenate
    cols <- colnames(imputed_ffs[[method]]@exprs)
    combined <- combined[,cols]
    combined@exprs <- rbind(combined[,cols]@exprs, imputed_ffs[[method]]@exprs)
    
    print('Building FlowSOM')
    fSOM <- FlowSOM(combined,
                    compensate = FALSE,
                    transform = FALSE,
                    colsToUse = c(backbone, variable), 
                    xdim = 15, 
                    ydim = 15,
                    nClus=40)

    # Export combined flowframe of ground-truth and imputed data with cluster labels
    clusters <- as.matrix(GetClusters(fSOM))
    clusters <- cbind(clusters, as.matrix(as.numeric(GetMetaclusters(fSOM))))
    colnames(clusters) <- c('fSOM_cluster', 'fSOM_metacluster')
    combined <- fr_append_cols(combined, clusters)
    
    # Re-add the NAs found in cyCombine
    # We add a false cluster of -1 to mark false-negatives
    if (method == 'cyCombine' && (dim(cyCombine_NAs)[1] != 0)){
      imp_state <- rep(1, nrow(cyCombine_NAs))
      cyCombine_NAs <- cbind(cyCombine_NAs, imp_state)
      fSOM_cluster <- rep(-1, nrow(cyCombine_NAs))
      cyCombine_NAs <- cbind(cyCombine_NAs, fSOM_cluster)
      fSOM_metacluster <- rep(-1, nrow(cyCombine_NAs))
      cyCombine_NAs <- cbind(cyCombine_NAs, fSOM_metacluster)
      combined@exprs <- rbind(combined@exprs, cyCombine_NAs)
      
      # Replace the NAs with an arbitrary high number to solve FCS format issues
      df <- data.frame(combined@exprs)
      df <- df %>% mutate_all(~ifelse(is.na(.), 999999, .))
      exprs <- as.matrix(df)
      colnames(exprs) <- gsub(".", "-", colnames(exprs), fixed=TRUE)
      combined@exprs <- exprs
    }

    write.FCS(combined, paste0(prefix, id, '_', method, '_exprs.fcs'))
  }
}
