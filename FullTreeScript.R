#load packages
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(rgl)
library(caTools)
library(magick)

#class(all_sixteen_tips) <- c("multiPhylo", class(all_sixteen_tips)) #add object class

all_filepath_trees <- list.files("~/Summer Scholarship 2019/EpichloeAll/EpiAllTrees", full.names = TRUE) #filepath to best trees
all_read_trees <- lapply(all_filepath_trees, read.tree)
all_full_tree <- sapply(all_read_trees, Ntip) == 23 #all tips, required for ASTRAL
all_23_trees <- all_read_trees[all_full_tree]
class(all_23_trees) <- c("multiPhylo", class(all_23_trees))
class(all_read_trees) <- c("multiPhylo", class(all_read_trees))
write.tree(all_sixteen_tips, file = "all_full_trees.phy")
write.tree(all_read_trees, file = "all_trees.phy") #newick format for ASTRAL (use write.nexus for nexus format)

#reinput ASTRAL tree to root + show node confidences
species_tree <- read.tree("~/Summer Scholarship 2019/EpichloeAll/allspeciestree")
rooted_tree <- root(species_tree, "gansuensis")
plot.phylo(rooted_tree, show.tip.label = TRUE, show.node.label = TRUE)

#make consensus net
cnet <- consensusNet(all_sixteen_tips, .2) #change for different relationships
plot(cnet, "2D", show.edge.label = TRUE)#can show.edge.label =TRUE for clarity
plot(cnet)
play3d(spin3d(axis=c(0,1,0), rpm=6), duration=10)
# create animated gif file - messy at current
movie3d(spin3d(axis=c(0,1,0), rpm=6), duration=10, dir = "~")

#produce dist_topo, see any obvious discrepancies / unpublishable
dist_topo <- dist.topo(all_sixteen_tips, method = "PH85")
tree <- upgma(dist_topo)
plot(tree)






#orginal graph using trees exactly same as species tree, with no examination with outgroup

unroot_all <- unroot(all_23_trees)
unroot_species <- unroot(species_tree)
dist_from_spp_tree <- function(tr, spp_tree = unroot_species){dist.topo(tr, spp_tree)}

D <- sapply(unroot_all, dist_from_spp_tree)
any(D == 0)
same_as_spp_best <- D==0
num_same_best <- which(same_as_spp_best == TRUE)

trees_same_best <- all_23_trees[num_same_best]
plot(trees_same_best[1])


num_full_best <- which(all_full_tree_best == TRUE) #list of number of those trees which are full
num_tree_best <- num_full_best[num_same_best]
filepaths_best <- all_filepath_trees[num_tree_best]
library(stringr)
og_num_best <- str_extract(string = filepaths_best, pattern = "og_\\d+")
#only four gene trees match the species tree! og_2010, 2263, 8538, 8671

#A useful pattern!
#og_num_only <- str_match(string = filepaths_best, pattern = "og_(\\d+)")


#now comparing those which are same as gene tree to ones which are different
#need to compare gene trees which are monophyletic for xxx + KMK to those which are not

bothlabels <- sapply(all_read_trees, function(tr) "C2857" %in% tr$tip.label & "Efe2368" %in% tr$tip.label)
inclusive <- which(bothlabels == TRUE) #number on gene tree with both labels
toanalyse <- all_read_trees[inclusive]

monophyly <- sapply(toanalyse, is.monophyletic, tips=c("C2857", "Efe2368"))
length(which(monophyly==TRUE))
includedfilepaths <- all_filepath_trees[inclusive]
og_num_both <- str_extract(string = includedfilepaths, pattern = "og_\\d+")
negatives <- which(monophyly==FALSE)
dif_between_festucae <- og_num_both[negatives]


#### isolate secodn column from protein ortho long file to use later on in festucae data??
#just kidding need geneid - second column is gene id ahh

protein_ortho_long <- read.delim("~/Summer Scholarship 2019/EpichloeAll/ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)

strainC2857 <- subset.data.frame(protein_ortho_long, V3=="C2857")
includedfilepaths <- all_filepath_trees[inclusive]
og_num_both <- str_extract(string = includedfilepaths, pattern = "og_\\d+")

inframe_both <- which(strainC2857$V1 %in% og_num_both)
C2857genes_both <- strainC2857[c(inframe_both),2]
write.table(C2857genes_both, quote = FALSE, row.names = FALSE, col.names = FALSE , file = "C2857genes_both.txt")

EALboth <- read.delim("~/Summer Scholarship 2019/Kakapo/EALboth.txt", stringsAsFactors=FALSE)
indatabase <- which(EALgenes_both %in% EALboth$Protein.stable.ID)
mono <- monophyly[indatabase]

EALboth$monophy = mono
plot <- ggplot(EALboth, aes(Gene.start..bp., Chromosome.scaffold.name))
points <- plot + geom_jitter(aes(colour = monophy), width = 0, height = 0.1)
fumigatusgraph <- points + labs(x = "Gene bp start", y = "Chromosome number", title = "A. fumigatus gene locations", colour = "Same as Species Tree")
scalechange <- fumigatusgraph + scale_x_continuous("Position in chromosome (Mbp)", labels = function(x) x/1e6)


