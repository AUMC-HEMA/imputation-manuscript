library(HDCytoData)
library(flowCore)

ff <- Nilsson_rare_flowSet(metadata = FALSE)@frames$V1

export_dir <- 'EXPORT PATH'

# Backbone
# PI, CD45, CD11b, CD34, CD19, CD3, CD38
backbone <- c('PE-Texas Red-A', 'Alexa Fluor 700-A', 'Qdot 605-A', 'PerCP-Cy5-5-A',
              'PE-Cy5-A', 'APC-Cy7-A', 'FITC-A')
# Tube 1: Lymphoid markers
# CD4, CD10, CD45RA
markers1 <- c('Qdot 655-A', 'PE-Cy7-A', 'DAPI-A')
# Tube 2: Myeloid markers
# CD123, CD90, CD49f, CD110
markers2 <- c('PE-A', 'QDot 800-A', 'QDot705-A', 'APC-A')
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
write.FCS(ff, paste0(export_dir, 'Nilsson_rare', '_gt.fcs'))
write.FCS(ff1, paste0(export_dir, 'Nilsson_rare', '_ff1.fcs'))
write.FCS(ff2, paste0(export_dir, 'Nilsson_rare', '_ff2.fcs'))
