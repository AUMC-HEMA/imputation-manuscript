# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(CytoML)
library(flowWorkspace)

wsp_dir <- '/LOCATION OF FLOWJO WORKSPACE FILE'
agg_dir <- '/LOCATION OF AGGREGATED FCS FILE'
export_dir <- '/WHERE TO SAVE .CSV OF CELL LABELS'

for (name in c('Flow_cyCombine', 'Flow_CyTOFmerge', 'Flow_CytoBackBone', 'Flow_Infinicyt')){
  print(name)
  wsfile <- paste0(wsp_dir, name, '_SingleGate.wsp')
  ws <- open_flowjo_xml(wsfile, sample_names_from = "keyword")
  gs <- flowjo_to_gatingset(ws, name=1, execute=FALSE, transform=FALSE)
  gh <- gs[[1]]
  ff <- read.FCS(paste0(agg_dir, name, '_agg.fcs'))
  gates <- list()
  pops <- gs_get_leaf_nodes(gs)
  for (pop in pops){
    # Add as last column the leaf gate
    gate <- gs_pop_get_gate(gs, pop)[[1]]
    filter <- filter(ff, gate)
    bools <- as.matrix(as.integer(filter@subSet))
    colnames(bools) <- pop
    gates[[pop]] <- bools
  }
  gates <- data.frame(gates)
  colnames(gates) <- pops
  write.csv(gates, paste0(export_dir, name, 'SingleGate_labels.csv'))
}
