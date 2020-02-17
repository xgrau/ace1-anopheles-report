### Define input ####

# input files
library(ape)
library(phytools)
library(stringr)

prefix    = "phylo"

phy_list  = c("hapalignment_upstream.iqt.treefile",
              "hapalignment_duplication.iqt.treefile",
              "hapalignment_downstream.iqt.treefile")
phy_name  = c("upstream","duplication","downstream")

dups = read.table("haps_with_dups.csv", header = T)
hap_has_dup = as.vector(dups[dups$has_dup,"hap"])

#  plot phylogenies
pdf(file=paste(prefix,"circol.pdf",sep="."),height=6,width=6)
cou      = 0
for (phy_fn in phy_list) {
  
  # tree to distance matrix
  cou    = cou+1
  phy    = read.tree(phy_fn)
  
  # label
  phy$tip.label_sep  = gsub("_"," ",phy$tip.label)
  phy$tip.species    = gsub("^[A-Z][A-Z]","",word(phy$tip.label_sep,2))
  phy$tip.species[phy$tip.species == ""] = "gamcol"
  phy$tip.population = word(phy$tip.label_sep,2)
  phy$tip.hasdup     = phy$tip.label %in% hap_has_dup
  phy$tip.colfactor  = as.factor(paste(phy$tip.species,phy$tip.hasdup,sep="_"))
  phy$tip.colors     = c("blue","turquoise3","orangered3","orange")[phy$tip.colfactor]

  # rename tips for nice plot
  phy$tip.label      = rep("·", length(phy$tip.label))
  
  # resize branches for nice plot
  # phy$edge.length[phy$edge.length >  0.005] = 0.005
  # phy$edge.length[phy$edge.length == 0]    = 5e-5
  
  # plot
  plot.phylo(phy, type="unr",
             use.edge.length=T, show.tip.label=T, show.node.label=F,
             tip.color = phy$tip.colors, lab4ut = "axial",
             edge.color = "slategray3",
             font = 1, edge.width = 0.5, node.depth=1, cex=2,
             main=phy_name[cou])
  ape::add.scale.bar(x=0, y=0, lcol="slategray", cex=0.8)
  
}

dev.off()


pdf(file=paste(prefix,"llarg.pdf",sep="."),height=150,width=30)
cou      = 0
for (phy_fn in phy_list) {
  
  # tree to distance matrix
  cou    = cou+1
  phy    = read.tree(phy_fn)
  
  # label
  phy$tip.label_sep  = gsub("_"," ",phy$tip.label)
  phy$tip.species    = gsub("^[A-Z][A-Z]","",word(phy$tip.label_sep,2))
  phy$tip.species[phy$tip.species == ""] = "gamcol"
  phy$tip.population = word(phy$tip.label_sep,2)
  phy$tip.hasdup     = phy$tip.label %in% hap_has_dup
  phy$tip.colfactor  = as.factor(paste(phy$tip.species,phy$tip.hasdup,sep="_"))
  phy$tip.colors     = c("blue","turquoise3","orangered3","orange")[phy$tip.colfactor]
  
  
  # plot
  plot.phylo(phy, type="phy",
             use.edge.length=T, show.tip.label=T, show.node.label=T,
             tip.color = phy$tip.colors, lab4ut = "axial",
             edge.color = "slategray3",
             font = 1, edge.width = 0.5, node.depth=1, cex=1.2, underscore = T,
             main=phy_name[cou])
  ape::add.scale.bar(x=0, y=0, lcol="slategray", cex=0.8)
  
}

dev.off()


stop("ara")

