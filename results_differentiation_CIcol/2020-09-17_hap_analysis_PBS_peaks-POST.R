#### Input ####

library("ape")
library("stringr")
library("phangorn")


fasta_list <- dir("results_peaks/", pattern =".fasta")

pdf("pbs_peak_phylogenies_nj.pdf", width = 3, height = 5)
for (fas in fasta_list) {

  # plot  
  fai = ape::read.FASTA(file = sprintf("results_peaks/%s", fas) ,type = "DNA")
  dis = ape::dist.dna(fai)
  tre = ape::nj(dis)
  tre = phangorn::midpoint(tre)
  fan = stringr::str_replace(fas, pattern = "alignment.fasta", "")
  fan = stringr::str_replace(fan, pattern = "^sel_", "")
  fan = stringr::str_replace_all(fan, pattern = "_", " ")
  
  # dataframe of edges
  tre_edge           = as.data.frame(tre$edge)
  colnames(tre_edge) = c("edge_start","edge_end")
  tre_edge$ix_edges = as.numeric(rownames(tre_edge))
  tre_edge$ends_in_tip = tre_edge$edge_end <= length(tre$tip.label)
  # dataframe of nodes
  tre_nods = data.frame(taxa = c(tre$tip.label, tre$node.label))
  tre_nods$edge_end = as.numeric(rownames(tre_nods))
  tre_nods$is_tip   = tre_nods$edge_end <= length(tre$tip.label)
  
  # merge them
  tre_edge = merge(tre_edge, tre_nods, all.x = T, all.y = T, by.x = "edge_end", by.y = "edge_end")
  
  # find phenotype info
  tre_phe = data.frame(row.names = tre$tip.label)
  tre_phe$taxa = rownames(tre_phe)
  tre_phe$phenotype = grepl("alive", tre_phe$taxa)
  tre_phe$color = "slategray"
  tre_phe[tre_phe$phenotype,"color"] = "springgreen3"
  tre_phe[!tre_phe$phenotype,"color"] = "deeppink3"

  # # add info to plot
  # tre_data = merge(tre_edge, tre_phe, by.x = "taxa", by.y = "taxa",all.x = T)
  # tre_data = tre_data[order(tre_data$ix_edges),]
  # tre_data["phenotype"][is.na(tre_data["phenotype"])] <- FALSE
  # 
  # # add colors
  # tre_data$color = "slategray"
  # tre_data[tre_data$phenotype  & tre_data$ends_in_tip,"color"] = "springgreen3"
  # tre_data[!tre_data$phenotype & tre_data$ends_in_tip,"color"] = "deeppink3"
  # tre_data$width = 1
  # tre_data[tre_data$ends_in_tip,"width"] = 2
  
  # tre$edge.length[tre$edge.length<1e-3] = 1e-3
  
  plot.phylo(tre, font=1, type="p", show.tip.label = F, edge.color = "slategray", root.edge = T)
  tiplabels(pch = 19, col = tre_phe$color, height = 4, cex = 0.4)
  title(sprintf("%s",fan), cex.main=0.8)
  
}
dev.off()