
args = commandArgs(trailingOnly=TRUE)

library(ape)


input_phylogenetic_tree <- args[1]
input_occurence_data    <- args[2]

phylogenetic_tree <- read.csv(input_phylogenetic_tree, dec = ".", sep = ";", header = T)
occurence_data    <- read.tree(input_occurence)

ajout des taxons et branches abscentes sur l''arbre


# create community matrix for phyloregions (special format for phyloreg)
long parse après ajouts des taxons sur l''arbre
transforme les données de communautées en matrice condencées .(long2sparse/pez)
#create comm matrix and phylogeny with only present species from community
compare les taxons et la phylogénie avec une matrice de communauté, tente d''ajouter les taxons manquants .(match_phylo_comm/phylodiversity)

This function computes the standard effect size of PD by correcting for changes in species richness.
The novelty of this function is its ability to utilize sparse community matrix making it possible to
efficiently randomize very large community matrices spanning thousands of taxa and sites.(PD_ses/phyloregion)

clacul MPD ET MNTD (picante)

phy.dist_pycnos  <- cophenetic(com.data.pycno$phy)
pycno_mpd        <- as.data.frame(picante::mpd(comm_pycnos, phy.dist_pycnos))

Computes the cophenetic distances for a hierarchical clustering .(cophenetic/stats)


library(phyloregion)

chosefile()

readfile()

ophiL2P <- long2sparse(pycno_grid, grid = "grids", species = "newscientificname")
com.data.pycno <- match_phylo_comm(pycno_tree_merge, sparse_grid3_pycnos, delete_empty_rows = F)
comm_pycnos    <- com.data.pycno$comm
ses_PD_pycnos <- PD_ses(comm_pycnos, com.data.pycno$phy, model = "tipshuffle", reps = 999)








