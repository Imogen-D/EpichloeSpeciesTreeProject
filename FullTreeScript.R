library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(rgl)
library(caTools)
library(stringr)
library(tidyr)
library(XML)
library("rentrez")
library(RColorBrewer)
library(magick)


#get filepaths to all RAxML bestTree genetrees
all_filepath_trees <- list.files("RAxMLTrees/", full.names = TRUE)
all_read_trees <- lapply(all_filepath_trees, read.tree)
all_full_tree <- sapply(all_read_trees, Ntip) == 24 #all tips, required for ASTRAL
all_twentyfour_tips <- all_read_trees[all_full_tree]
class(all_twentyfour_tips) <- c("multiPhylo", class(all_twentyfour_tips))
class(all_read_trees) <- c("multiPhylo", class(all_read_trees))
write.tree(all_twentyfour_tips, file = "all_full_trees.phy")
write.tree(all_read_trees, file = "all_trees.phy") #newick format for ASTRAL

#for this code iteration for others
all_twentyfour_tips <- read.tree("FULLGENETREES") #multiphylo object supplied in repository
all_full_trees <- read.tree("ALLGENETREES") #multiphylo object supplied in repository

#all_full_trees.phy is utilised in ASTRAL, windows command prompt. Species map file was used for species tree and outlog produced.

#reinput ASTRAL tree to root + show node confidences
strain_tree <- read.tree("straintree") #supplied in repository
strain_rooted_tree <- root(strain_tree, "CCE27021")
plot.phylo(strain_rooted_tree, show.tip.label = TRUE, show.node.label = TRUE) #is exported and saved

species_tree <- read.tree("speciestree") #species tree used species map file, supplied in repository
rooted_tree <- root(species_tree, "C.purpurea")
plot.phylo(rooted_tree, show.tip.label = TRUE, show.node.label = TRUE) #is exported and saved

#make consensus net
cnet <- consensusNet(all_twentyfour_tips, .2)
plot(cnet, "2D", show.edge.label = FALSE)
plot(cnet) #exported

#produce dist_topo, see any obvious discrepancies / unpublishable
dist_topo <- dist.topo(all_twentyfour_tips, method = "PH85")
dist <- upgma(dist_topo)
plot(dist)



protein_ortho_long <- read.delim("ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)

#monophyly studies
#FESTUCAE
festlabels <- sapply(all_read_trees, function(tr) "C2857" %in% tr$tip.label & "Efe2368" %in% tr$tip.label)
festinclusive <- which(festlabels == TRUE) #number on gene tree with both labels
festtoanalyse <- all_read_trees[festinclusive]

festmonophyly <- sapply(festtoanalyse, is.monophyletic, tips=c("C2857", "Efe2368"))
length(which(festmonophyly==TRUE))
festincludedfilepaths <- all_filepath_trees[festinclusive]
festog_num_both <- str_extract(string = festincludedfilepaths, pattern = "og_\\d+")
festnegatives <- which(festmonophyly==FALSE)
dif_between_festucae <- festog_num_both[festnegatives]

festdif_filepaths <- festincludedfilepaths[festnegatives]
festdifreadtrees <- lapply(festdif_filepaths, read.tree)
class(festdifreadtrees) <- c("multiPhylo", class(festdifreadtrees))
TFfullfesttrees <- sapply(festdifreadtrees, Ntip) == 24 
all_full_fest_trees <- festdifreadtrees[TFfullfesttrees]
write.tree(festdifreadtrees, file = "festASTRALinput.phy")

#ASTRAL, saved as festucaedifferences

fest_tree <- read.tree("festucaedifferences")
fest_rooted_tree <- root(fest_tree, "CCE27021")
plot.phylo(fest_rooted_tree, show.tip.label = TRUE, show.node.label = TRUE)

festcnet <- consensusNet(all_full_fest_trees, .2)
plot(festcnet, "2D", show.edge.label = FALSE)




#BROMICOLA
bromlabels <- sapply(all_read_trees, function(tr) "EbroNfe1" %in% tr$tip.label & "EbroAL0434" %in% tr$tip.label & "EbroAL0426" %in% tr$tip.label)
brominclusive <- which(bromlabels == TRUE) #number on gene tree with both labels
bromtoanalyse <- all_read_trees[brominclusive]

brommonophyly <- sapply(bromtoanalyse, is.monophyletic, tips=c("EbroNfe1", "EbroAL0434", "EbroAL0426"))
bromincludedfilepaths <- all_filepath_trees[brominclusive]
bromog_num <- str_extract(string = bromincludedfilepaths, pattern = "og_\\d+")
bromnegatives <- which(brommonophyly==FALSE)
dif_between_bromicola <- bromog_num[bromnegatives]

bromdif_filepaths <- bromincludedfilepaths[bromnegatives]
bromdifreadtrees <- lapply(bromdif_filepaths, read.tree)
class(bromdifreadtrees) <- c("multiPhylo", class(bromdifreadtrees))
TFfullbromtrees <- sapply(bromdifreadtrees, Ntip) == 24 
all_full_brom_trees <- bromdifreadtrees[TFfullbromtrees]
write.tree(bromdifreadtrees, file = "bromASTRALinput.phy")

#ASTRAL, saved as bromicoladifferences
                     
brom_tree <- read.tree("bromicoladifferences")
brom_rooted_tree <- root(brom_tree, "CCE27021")
plot.phylo(brom_rooted_tree, show.tip.label = TRUE, show.node.label = TRUE)

bromcnet <- consensusNet(all_full_brom_trees, .2)
plot(bromcnet, "2D", show.edge.label = FALSE)




#TYPHINA
typhlabels <- sapply(all_read_trees, function(tr) "E8Q19" %in% tr$tip.label & "E8Q16" %in% tr$tip.label & "Ety8" %in% tr$tip.label)
typhinclusive <- which(typhlabels == TRUE) #number on gene tree with both labels
typhtoanalyse <- all_read_trees[typhinclusive]

typhmonophyly <- sapply(typhtoanalyse, is.monophyletic, tips=c("E8Q19", "E8Q16", "Ety8"))
typhincludedfilepaths <- all_filepath_trees[typhinclusive]
typhog_num <- str_extract(string = typhincludedfilepaths, pattern = "og_\\d+")
typhnegatives <- which(typhmonophyly==FALSE)
dif_between_typhina <- typhog_num[typhnegatives]

typhdif_filepaths <- typhincludedfilepaths[typhnegatives]
typhdifreadtrees <- lapply(typhdif_filepaths, read.tree)
class(typhdifreadtrees) <- c("multiPhylo", class(typhdifreadtrees))
TFfulltyphtrees <- sapply(typhdifreadtrees, Ntip) == 24 
all_full_typh_trees <- typhdifreadtrees[TFfulltyphtrees]
write.tree(typhdifreadtrees, file = "typhASTRALinput.phy")

#ASTRAL, saved as bromicoladifferences

typh_tree <- read.tree("typhinadifferences")
typh_rooted_tree <- root(typh_tree, "CCE27021")
plot.phylo(typh_rooted_tree, show.tip.label = TRUE, show.node.label = TRUE)

typhcnet <- consensusNet(all_full_typh_trees, .2)
plot(typhcnet, "2D", show.edge.label = FALSE)


#producing metadata heatmaps
                     
#SEXUAL REPRODUCTION, SUBTRIBE, GEOGRAPHY

epi_info <- read.csv("EpichloeTaxaComma.csv") #where csv is knowledge of host species
names(epi_info)[names(epi_info) == "ï..name"] <- "Species"

unroot_tree <- drop.tip(rooted_tree, "C.purpurea")
x <- plot(unroot_tree)
ggtree <- ggtree(unroot_tree)
joined <- ggtree %<+% epi_info
x <- joined +
  geom_tiplab(aes(color = SexualReproduction)) +
  theme(legend.position = "right")
  scale_color_brewer(palette = "Set1") +
  labs(title="Sexual Reproduction in Epichloe")
par(mar = c(5, 5, 5, 5))
x
ggsave("SexualRepro.pdf", width = 50, height = 30, units = "cm", limitsize = FALSE)






geo_mat <- epi_info[,12:17] #geography
row.names(geo_mat) <- epi_info$Species

cols= c("darkgrey", "forestgreen")
gheatmap(x, geo_mat, offset = 1, colnames_angle = 90, colnames_offset_y = -0.5, width = 0.3) + scale_fill_manual(values=cols) + theme(legend.title = element_text()) + labs(fill = "Distribution") + labs(title = "Sexual reproduction and global distribtion of Epichloe")  
#heatmap of geography
ggsave("Geo.pdf", width = 50, height = 30, units = "cm", limitsize = FALSE)


HostNumber <- str_count(epi_info[,"KnownHostRange"], ",") + 1  #number of hosts

hosts <- read.csv("EpiHosts.csv")

WideHostGenus <- read.csv("WideHostGenus.csv", row.names=1, stringsAsFactors=FALSE)


class(hosts$Host[[1]])
taxon_names <- sapply(hosts$Host, rotl::tnrs_match_names)
planttaxainfo <- as.data.frame(taxon_names)

removed_26flags <- taxon_names[,-26]
removed_23flags <- removed_26flags[,-23] #ISSUES WITH ROTL TAXON FINDING REMOVED
id <- as.integer(removed_23flags[4,]) 
first_tree <- rotl::tol_induced_subtree(ott_ids = id)
as.character(planttaxainfo[4,])


#making heatmap with host Genus
hostplot <- gheatmap(x, WideHostGenus, offset = 1.5, colnames_angle = 90, colnames_offset_y = -2, color = "black", width = 0.5) + 
  scale_fill_brewer(palette = "Paired") + 
  theme(legend.title = element_text()) + 
  labs(fill = "Host Presence") + 
  labs(title = "Sexual reproduction and host genera of Epichloe") +
  geom_tiplab(align=TRUE)  

ggsave("HostsGenusHeatmap.pdf", plot = last_plot(), height = 12, width = 28)


Genera <- colnames(WideHostGenus)

get_tribe <- function(Genu){
  taxon_search <- entrez_search(db="taxonomy", term = Genu)
  taxon_rec <- entrez_fetch(db="taxonomy", id=taxon_search$ids, rettype="xml", parsed=TRUE)
  get_taxon(taxon_rec, "tribe")
}

Genus <- hosts$Host
tribes <- sapply(Genus, get_tribe)

hosts$tribes = tribes

WideTribe <- read.csv("WideTribe.csv", row.names=1, stringsAsFactors=FALSE)


tribeheatmap <- gheatmap(x, WideTribe, offset = 1.5, colnames_angle = 90, colnames_offset_y = -2, color = "black", width = 0.5) +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.title = element_text()) +
  labs(fill = "Host Presence") +
  labs(title = "Sexual reproduction and hosts of Epichloe") + 
  geom_tiplab(align=TRUE)
ggsave("HostTribeHeatmap.pdf", plot = last_plot(), height = 12, width = 28)


#CYCLIC PEPTIDES
                     
pa <- read.delim("pa.csv", header=FALSE, stringsAsFactors=FALSE)
#CSV describing presence of orthologs
                     
unroot_strain_tree <- drop.tip(strain_rooted_tree, "Cpur")
x <- plot(unroot_strain_tree)
ggstraintree <- ggtree(unroot_strain_tree)
joinedpa <- ggstraintree %<+% pa
y <- joinedpa +
  geom_tiplab(aes()) +
  theme(legend.position = "right")
  scale_color_brewer(palette = "Set1")


paspread <- spread(pa, V1, V3)
rownames(paspread) <- paspread$V2
paspread$V2 <- NULL

cols= c("darkgrey", "forestgreen")
gheatmap(y, paspread, offset = 4, colnames_angle = 90, colnames_offset_y = -0.5, width = 1.5) + scale_fill_manual(values=cols) + theme(legend.title = element_text()) + labs(fill = "Distribution")  

