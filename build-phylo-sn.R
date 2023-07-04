### Building a phylogeny for the species involved in pollen size/number analysis
# Claire Smith
# Last updated: 4 July 2023

library(V.PhyloMaker2)
library(ape)
library(dplyr)

# Read in list of species, genera, families...
sn.sp.list <- read.csv("raw-data/species_list_sizenum.csv")

# Create phylogeny using V.Phylomaker2
sn.phylo <- V.PhyloMaker2::phylo.maker(sn.sp.list, tree = GBOTB.extended.TPL)

### plot the phylogenies with node ages displayed.
# par(mfrow = c(1, 3))
plot.phylo(sn.phylo$scenario.3, cex = 1.5, main = "scenario.3")
# nodelabels(round(branching.times(sn.phylo$scenario.3), 1), cex = 1)
# dev.off()

# Default is scenario 3 -- not sure what difference it makes? The diff
# scenarios just seem to rearrange exact order of sedges, and change branch lengths a little

sn.tree <- sn.phylo$scenario.3
sn.tree$node.label<-NULL # without this caper:pgls throws an error

write.tree(sn.tree, file = "processed-data/sizenum_tree.nwk")

