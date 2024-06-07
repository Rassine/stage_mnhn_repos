library(rgeos) 
library(rgdal)
library(sp)
library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(picante)  
library(phyloregion)
library(tidyr) 
library(pez)
library(vegan)


# Import you phylogeny and occurrence data table
actino_tree <- file.choose()# choose your file
actino_tree <- read.tree(actino_tree)
plot(actino_tree)

actinos_gridPF_1 <- file.choose()  # choose your file actinos_gridPF_1_2024
actinos_gridPF_1 <- read.csv2(actinos_gridPF_1, dec = ".", sep = ";", header = T, na.strings = "")  #read your file


# Import the 3x3 degrees grid cell file cutted by Polar Front (PF zone), 
#in this study one grid cell is considered as a species community
#Lambert Azimuthal Equal Area projection 'laea'

grid_aq_3_PF <- read_sf("grid_aq_3_PF.shp", layer = "grid_aq_3_PF")
prj_laea     <- "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
grid_aq_3_PF <- st_transform(grid_aq_3_PF, crs = prj_laea)

plot(st_geometry(grid_aq_3_PF))


#create a community matrix ----
# create a sparce, more light community matrix than the classic one

sparse_grid3_actinos <- long2sparse(actinos_gridPF_1, grid = "grids", species = "newscientificname")
class(sparse_grid3_actinos)
dim(sparse_grid3_actinos)  # check on matrix dimension (species and cells number)
rownames(sparse_grid3_actinos)

# Calculate PD measure and its z-value
set.seed(34)  # random seed for SES calculations and zscore

# verify the correspondence btw species on a tree tips and in community matrix
com.data.actinos <- match_phylo_comm(actino_tree, sparse_grid3_actinos, delete_empty_rows = F)
dim(com.data.actinos$comm)  # check on how many species or cells are gone

# put the community matrix to a data frame format
comm_actinos <- com.data.actinos$comm
comm_actinos <- as.data.frame(as.matrix(comm_actinos))

# check on phylogeny
com.data.actinos$phy

# Calculate PD measure and its z-value 
ses_PD_actinos <- PD_ses(com.data.actinos$comm, com.data.actinos$phy, model = "tipshuffle", reps = 999)

# Merge each grid cell with produced PD values 
PD_actinos <- merge(grid_aq_3_PF, ses_PD_actinos, by = "grids", all.x = T, all.y = F)

# delete NA cells
PD_actinos <- PD_actinos %>% filter(!is.na(PD_obs))  #no NA allowed for following models

## richness column created by PD_ses is not species richness 
#but total abundance in each grid cell !!! 
# calculate Species richness differently, with picante
pd_actinos_picante       <- as.data.frame(picante::pd(comm_actinos, com.data.actinos$phy, include.root = F))
pd_actinos_picante$grids <- rownames(pd_actinos_picante)
PD_actinos               <- merge(PD_actinos, pd_actinos_picante, by = "grids", all.x = T, all.y = F)

# delete cells with one species
PD_actinos <- PD_actinos %>% filter(!is.na(PD))

# Calculate endemism
PE <- phylo_endemism(com.data.actinos$comm, com.data.actinos$phy)
WE <- weighted_endemism(com.data.actinos$comm)
PD_actinos <- merge(PD_actinos, data.frame(grids = names(PE), PE = PE), by = "grids")
PD_actinos <- merge(PD_actinos, data.frame(grids = names(WE), WE = WE), by = "grids")

# standardize phylogeneitc endemism by weighted occurrence based endemesm with loess regression
PD_actinos$std_PE <- as.vector(scale(resid(loess(PE ~ WE, data=PD_actinos))))

# stndardize PD with species richness by loess non-linear regression ----
PD_actinos$std_resid_loess <- as.vector(scale(resid(loess(PD_obs ~ SR, data = PD_actinos))))

# stndardize PD with species richness by GAM non-linear regression ----
library(mgcv)
PD_actinos$std_resid_gam <- as.vector(scale(resid(gam(PD_obs ~ s(SR), data = PD_actinos))))

## Calculate MPD and MNTD indices ----
# use the same comm_pycnos community matrix
# calculate phenetic distances between species from phylogeny needed for these two indices calculation
phy.dist_actinos <- cophenetic(com.data.actinos$phy)




# MPD calculation
actino_mpd        <- as.data.frame(picante::mpd(comm_actinos, phy.dist_actinos))
names(actino_mpd) <- "mpd"
actino_mpd$grids  <- PD_actinos$grids       # add grid cells column to merge with global PD_actinos table


# randomize MPD with Standard Effect Size (SES) - and get the zscore (corrected version of MPD)
ses_mpd_actino       <- ses.mpd(comm_actinos, phy.dist_actinos, null.model = "taxa.labels", runs = 999)
ses_mpd_actino$grids <- row.names(ses_mpd_actino)
ses_actinos          <- merge(grid_aq_3_PF, ses_mpd_actino, by = "grids", all.x = T, all.y = F)



# same for MNTD index
actino_mntd        <- as.data.frame(picante::mntd(comm_actinos, phy.dist_actinos))
names(actino_mntd) <- "mntd"
actino_mntd$grids  <- PD_actinos$grids


ses_mntd_actino       <- ses.mntd(comm_actinos, phy.dist_actinos, null.model = "taxa.labels", runs = 999)
ses_mntd_actino$grids <- row.names(ses_mntd_actino)
ses_actinos           <- merge(ses_actinos, ses_mntd_actino, by = "grids", all.x = T, all.y = F)



ses_actinos <- ses_actinos %>% filter(!is.na(ses_actinos$ntaxa.x))

PD_actinos <- merge(PD_actinos, actino_mpd,  by = "grids", all.x = T, all.y = F)
PD_actinos <- merge(PD_actinos, actino_mntd, by = "grids", all.x = T, all.y = F)
PD_actinos <- merge(PD_actinos, st_drop_geometry(ses_actinos), by = "grids", all.x = T, all.y = F)

# Create pairwise matrix of scatterplots between all indices calculated
PD_actinos_noXY <- PD_actinos %>%
                   st_drop_geometry(PD_actinos) %>%
                   select( grids 
                         , PD 
                         , zscore 
                         , std_resid_loess 
                         , std_resid_gam 
                         , SR 
                         , mpd 
                         , mntd 
                         , mpd.obs.z 
                         , mntd.obs.z 
                         , PE 
                         , WE 
                         , std_PE,)

str(PD_actinos_noXY)
pairs( PD_actinos_noXY
     , lower.panel = panel.smooth
     , upper.panel = panel.cor
     , cex = 1.5
     , pch = 1
     , cex.labels = 1
     , gap = 1/4)

## CREATE PHYLOREGIONS ----
library(phyloregion)
# check for species community matrix and phylogeny correspondence
com.data.actinos <- match_phylo_comm(actino_tree, sparse_grid3_actinos, delete_empty_rows = F)

# calculate phylogenetic Beta diversity - a phylogenetic distance matrix between grid cells
phylo_beta_actinos <- phylobeta(com.data.actinos$comm, com.data.actinos$phy, index.family = "sorensen")
class(phylo_beta_actinos[[1]])
names(phylo_beta_actinos[[1]])

#select the less distorting clustering method, best fitting between phylogenetic distances in phylobeta matrix
# and raw distances from branch lengths of the tree
select_linkage(phylo_beta_actinos[[1]])
# select the most correlated method, probably would be UPGMA

#select optimal number of clusters with selected method
optim <- optimal_phyloregion(phylo_beta_actinos[[3]], k = 30)
plot(optim$df$k, optim$df$ev, pch = 20)  # k - nbr of clusters VS explained variance given k
# k has to be selected by a user
# pass the grid cell to spatial format 
grid_aq_3_sp <- as_Spatial(grid_aq_3_PF)

# calculate phyloregions clusters
y <- phyloregion(phylo_beta_actinos[[3]], shp = grid_aq_3_sp, k = 8, method = "average")
summary(y)
phylo_nmds <- y$NMDS

# take an shp spatial file for phyloregions and put it to sf format
phyloreg_sf <- y$shp
st_crs(phyloreg_sf)
plot(phyloreg_sf)
actinos_phyloreg_sf <- st_as_sf(phyloreg_sf)
# plot the evolutionaty distincivness of phyloregions, just serves to visualise 
plot(actinos_phyloreg_sf["ED"])

PD_actinos <- merge(PD_actinos, y$region.df, by = "grids", all.x = T, all.y = F)
names(PD_actinos)

names(y$region.df)
table(y$region.df$COLOURS)
#ed_palette <- c(unique(y$region.df$COLOURS))
PD_actinos <- PD_actinos %>% select(!c(cluster, ED, COLOURS))
PD_actinos <- merge(PD_actinos, y$region.df, by = "grids", all.x = T, all.y = F)

# calculate mean index values of all grid cells for each phylogerion
# Applying group_by & summarise
PD_actinos %>% group_by(cluster, .add = T) %>% summarise(mean_stdPD = mean(std_resid_loess), mean_stdPE = mean(std_PE), mean_MPD = mean(mpd), mean_ED = mean(ED))
# check visually with boxplots
library(ggpubr)
ggboxplot(PD_actinos, x = "cluster", y = "std_resid_loess")
ggboxplot(PD_actinos, x = "cluster", y = "std_PE")
ggboxplot(PD_actinos, x = "cluster", y = "mpd")

#is there any significat differnce btw clusters ?
library(rstatix)
res.kruskal <- PD_actinos %>% kruskal_test(std_resid_loess ~ cluster)
# or
kruskal.test(PD_actinos$std_resid_loess, g = PD_actinos$cluster)

# get the NMDS between phylogerions and plot it
phylo_nmds <- y$NMDS
phylo_nmds <- y$NMDS$points
class(phylo_nmds)
phylo_nmds <- as.data.frame(phylo_nmds)
phylo_nmds$phyloreg <- rownames(phylo_nmds)

# create a tiff figure
tiff("nmds_phyloregions_actinos.tiff", width = 6, height = 5, res = 600, units = "in", compression = "lzw")

#if you want the figure out with tiff delete "actino_nmds <- " wtich you'll need for global plot creation
actino_nmds <- ggplot()                                                                                        +
               geom_point(data = phylo_nmds, aes(x = MDS1, y = MDS2, fill = phyloreg), pch = 21, size = 6)     +
               scale_fill_brewer(palette = "Set3")                                                             +   # customized color palette
               labs(fill = "")                                                                                 +
               geom_text(data = phylo_nmds, aes(x = MDS1, y = MDS2, label = phyloreg), vjust = 0.3, size = 4)  +
               theme_minimal()                                                                                 +
               theme(legend.position = "none")
dev.off()


