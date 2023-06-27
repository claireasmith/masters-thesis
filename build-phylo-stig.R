### Building a phylogeny for the species involved in pollen capture
# Claire Smith
# Last updated: 27 June 2023

library(V.PhyloMaker2)
library(ape)
library(dplyr)

# Read in list of species, genera, families...
capture.sp.list <- read.csv("raw-data/species_list_capture.csv")

# Create phylogeny using V.Phylomaker2
capture.phylo <- V.PhyloMaker2::phylo.maker(capture.sp.list, tree = GBOTB.extended.TPL)

### plot the phylogenies with node ages displayed.
# par(mfrow = c(1, 3))
plot.phylo(capture.phylo$scenario.3, cex = 1.5, main = "scenario.3")
# nodelabels(round(branching.times(capture.phylo$scenario.3), 1), cex = 1)
# dev.off()

# Default is scenario 3 -- not sure what difference it makes? The diff
# scenarios just seem to rearrange exact order of sedges, and change branch lengths a little

capture.tree <- capture.phylo$scenario.3
capture.tree$node.label<-NULL # without this caper:pgls throws an error

write.tree(capture.tree, file = "processed-data/capture_tree.nwk")

