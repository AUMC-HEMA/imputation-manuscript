library(HDCytoData)
library(flowCore)

ff <- Mosmann_rare_flowSet(metadata = FALSE)@frames$V1

export_dir <- 'EXPORT PATH'

# Backbone
# Live/Dead, CD3, CD4, CD45RA, CD8a, CD14, CD69
backbone <- c('Violet G 550/40-A', 'Violet D 605/40-A', 'Green D 610/20-A',
              'Violet C 660/40-A', 'Violet B 705/70-A', 'Violet A 780/60-A',
              'Red A 780/60-A')
# Tube 1
# TNFa, IL17a, CXCR5, IL2
markers1 <- c('Violet H 450/50-A', 'Blue A 710/50-A', 'Blue B 515/20-A',
              'Red B 710/50-A')
# Tube 2
# IFNg, GZB-SA, CCL4, IL5
markers2 <- c('Green A 780/40-A', 'Violet E 585/42-A', 'Green E 575/25-A',
              'Red C 660/20-A')
variable <- c(markers1, markers2)

# Transform
ff <- transform(ff,transformList(c(backbone, variable), 
                                 arcsinhTransform(a = 0, b = 1/150, c = 0)))

# Sample flowframe to equal number of cells
# This prevents some issues with non-even dataframe manipulations later on
n = nrow(exprs(ff))%/%1000*1000
idx <- sample(seq(1, nrow(exprs(ff))), n, replace=FALSE)
ff@exprs <- ff@exprs[idx,]

# Add column designating original ordering
original_ID <- as.matrix(seq(1, nrow(exprs(ff))))
colnames(original_ID) <- 'original_ID'
ff <- fr_append_cols(ff, original_ID)

# Split into marker sets, each containing half of the data
sample_size <- nrow(exprs(ff)) / 2
ff1_idx <- sample(seq(1, nrow(exprs(ff))), sample_size, replace=FALSE)
ff1_idx <- sort(ff1_idx)
ff2_idx <- setdiff(seq(1, nrow(exprs(ff))), ff1_idx)
ff1 <- ff[ff1_idx,c(backbone, markers1, 'original_ID')]
ff2 <- ff[ff2_idx,c(backbone, markers2, 'original_ID')]

# Create a marker to keep track of imputed marker
dataset <- as.matrix(seq(1, nrow(exprs(ff))))
colnames(dataset) <- 'dataset'
dataset[ff1_idx,] <- 1
dataset[ff2_idx,] <- 2
ff <- fr_append_cols(ff, dataset)

# Export datasets
write.FCS(ff, paste0(export_dir, 'Mosmann_rare', '_gt.fcs'))
write.FCS(ff1, paste0(export_dir, 'Mosmann_rare', '_ff1.fcs'))
write.FCS(ff2, paste0(export_dir, 'Mosmann_rare', '_ff2.fcs'))