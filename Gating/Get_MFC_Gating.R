# Author: Tim Mocking
# Contact: t.r.mocking@amsterdamumc.nl
library(flowCore)
library(CytoML)
library(flowWorkspace)

wsp_dir <- '/LOCATION OF FLOWJO WORKSPACE FILE'
agg_dir <- '/LOCATION OF AGGREGATED FCS FILE'
export_dir <- '/WHERE TO SAVE .CSV OF CELL LABELS'

for (name in c('Flow_cyCombine', 'Flow_CyTOFmerge', 'Flow_CytoBackBone', 'Flow_Infinicyt')){
  wsfile <- paste0(wsp_dir, name, '.wsp')
  ws <- open_flowjo_xml(wsfile, sample_names_from = "keyword")
  gs <- flowjo_to_gatingset(ws, name=1, execute=FALSE, transform=FALSE)
  gh <- gs[[1]]
  
  ff <- read.FCS(paste0(agg_dir, name, '_agg.fcs'))
  indices <- list()
  pops <- gs_get_leaf_nodes(gs)
  for (pop in pops){
    # Remove last leaf in case of PD-1 pos
    if (pop == '/CD4 positive/NAIVE_SCM/NAIVE/PD1 positive NAIVE'){
      pop <- "/CD4 positive/NAIVE_SCM/NAIVE"
    }
    if (pop == '/CD4 positive/NAIVE_SCM/SCM/PD1 positive stem cell memories'){
      pop <- "/CD4 positive/NAIVE_SCM/SCM"
    }
    parents <- list()
    parent_nodes <- strsplit(pop, '/')[[1]]
    parent_nodes <- parent_nodes[-length(parent_nodes)]
    # From every leaf-node add all parent nodes in 0,1 format to a matrix
    for (parent in parent_nodes){
      if (parent != ''){
        print(parent)
        gate <- gs_pop_get_gate(gs, parent)[[1]]
        filter <- filter(ff, gate)
        bools <- as.matrix(as.integer(filter@subSet))
        colnames(bools) <- parent
        parents[[parent]] <- bools
      }
    }
    parents <- do.call(cbind, parents)
    # Add as last column the leaf gate
    gate <- gs_pop_get_gate(gs, pop)[[1]]
    filter <- filter(ff, gate)
    bools <- as.matrix(as.integer(filter@subSet))
    colnames(bools) <- pop
    gates <- cbind(parents, bools)
    n_gates <- ncol(gates)
    # Identify indices where all gates are TRUE including leaf
    gates <- cbind(gates, all_total = rowSums(gates[,1:n_gates]))
    # Identify indices where all gates are TRUE except leaf
    gates <- cbind(gates, parents_total = rowSums(gates[,1:n_gates-1]))
    
    # Find indices positive for all gates (including leaf gate)
    pos_idx <- which(gates[,'all_total']==n_gates, arr.ind = T)
    # Find all indices which are only negative for leaf gate
    neg_idx = which((gates[,'parents_total']==n_gates-1 & gates[,'all_total']!=n_gates), arr.ind=T)
    # Add to list
    indices[[paste0(pop, '+')]] <- pos_idx
    indices[[paste0(pop, '-')]] <- neg_idx
    # Recursively add all negative gates
    if (pop == '/CD8 positive/CCR7, CD45RA subset/CD57, CD28 subset/KLRG1, TIGIT subset/CD27-, TIGIT subset'){
      gates <- cbind(gates, level1_total = rowSums(gates[,1:3]))
      gates <- cbind(gates, level1_parents = rowSums(gates[,1:2]))
      neg_idx = which((gates[,'level1_parents'] == 2 & gates[,'level1_total'] != 3), arr.ind=T)
      indices[[paste0('/CD8 positive/CCR7, CD45RA subset/CD57, CD28 subset-')]] <- neg_idx
      gates <- cbind(gates, level2_total = rowSums(gates[,1:4]))
      gates <- cbind(gates, level2_parents = rowSums(gates[,1:3]))
      neg_idx = which((gates[,'level2_parents'] == 3 & gates[,'level2_total'] != 4), arr.ind=T)
      indices[[paste0('/CD8 positive/CCR7, CD45RA subset/CD57, CD28 subset/KLRG1, TIGIT subset-')]] <- neg_idx
    }
  }
  # Create final labeling
  labels <- as.matrix(rep('Other', nrow(ff@exprs)))
  for (pop in names(indices)){
    labels[indices[[pop]],] <- pop
  }
  write.csv(labels, paste0(export_dir, name, '_labels.csv'))
}
