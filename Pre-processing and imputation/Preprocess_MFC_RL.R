# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(dplyr)

# Location of pre-gated files
fcs_dir <- ''
export_dir <- ''

gt <- read.FCS(paste0(fcs_dir, 'GROUND TRUTH FILE'))
ff1 <- read.FCS(paste0(fcs_dir, 'FF1'))
ff2 <- read.FCS(paste0(fcs_dir, 'FF2'))

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

# Compensate files
gt <- compensate(gt, gt@description$SPILL)
ff1 <- compensate(ff1, ff1@description$SPILL)
ff2 <- compensate(ff2, ff2@description$SPILL)

# Transform
# Transform the data
trans_list <- list('100'=c('BV605-A'),
                   '150'=c('FITC-A', 'PerCP-Cy5-5-A', 'APC-R700-A', 'APC-H7-A', 
                           'BV421-A', 'HV500c-A', 'BV711-A', 'BV786-A', 'BUV395-A',
                           'BUV496-A', 'BUV737-A', 'PE-A', 'PE-CF594-A', 'PC7-A'),
                   '300'=c('APC-A'))

for (cofactor in names(trans_list)){
  gt <- transform(gt,transformList(trans_list[[cofactor]],
                                   arcsinhTransform(a = 0, 
                                                    b = 1/as.integer(cofactor),
                                                    c = 0)))
  ff1 <- transform(ff1,transformList(trans_list[[cofactor]],
                                     arcsinhTransform(a = 0, 
                                                      b = 1/as.integer(cofactor),
                                                      c = 0)))
  ff2 <- transform(ff2,transformList(trans_list[[cofactor]],
                                     arcsinhTransform(a = 0, 
                                                      b = 1/as.integer(cofactor),
                                                      c = 0)))}
# Subset the actual markers
ff1 <- ff1[,c(backbones$`Leukocytes_T-cells_Monocytes_T-sub_T-diff`, markers1)]
ff2 <- ff2[,c(backbones$`Leukocytes_T-cells_Monocytes_T-sub_T-diff`, markers2)]

write.FCS(gt, paste0(export_dir, '001_gt.fcs')
write.FCS(ff1, paste0(export_dir, '001_ff1.fcs')
write.FCS(ff2, paste0(export_dir, '001_ff2.fcs')
