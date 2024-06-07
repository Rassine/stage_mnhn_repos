library(ape)
library(phylobase)
library(phytools)
#library(treebase)
#library(geiger)
#library(adiv)
library(dplyr)
library(sf)
library(adephylo)
library(phylocomr)
library(pez)
library(phyloregion)
library(picante)
library(vegan)
  # Pycno_tree  ###### all fossils correlated_Anna (un de 3 phylos)
# in hundreds of Millions of years
# Ballesteros et al, 2020 doi:10.1093/molbev/msaa228
tree<-file.choose()
pycno_tree<-read.nexus(tree)
peycno_tree_OTL<-read.newick(tree)
plot.phylo(peycno_tree_OTL, show.node.label = TRUE, type="phylogram",
           cex=0.35, label.offset=3)
## Ballestros et al, 110 species
plot.phylo(pycno_tree, show.node.label = TRUE, type="phylogram",
           cex=0.35, label.offset=3)
axisPhylo(cex.axis=0.75)

plotTree(pycno_tree, node.numbers = T)
is.ultrametric(pycno_tree)

## make tree ultrametric ----
#dist_pycno<-cophenetic.phylo(pycno_tree)
tip.heights<-adephylo::distRoot(pycno_tree)
heights.summary<-table(tip.heights)# there are several one, should be one !
options(digits=22) # set to maximum allowed digits
real.tree.height<-as.numeric(names(which.max(heights.summary)))
over.under<-tip.heights-real.tree.height
class(over.under)
## we can now paint the branches that were problematic using phytools()
painted.tree<-paintBranches(pycno_tree[which(round(over.under,5)!=0)],"2") 
# here I am rounding to 5 decimal places... pretty arbitrary choice
#plotSimmap(painted.tree,lwd=4)
## extract all terminal edges for tips that do not have the final height we want:
tip.ids <- pycno_tree$edge[, 2] <= Ntip(pycno_tree)
terminal.edges <- pycno_tree$edge.length[tip.ids]
## add/subtract the extra length from the terminal branches
corrected.terminal.edges<-terminal.edges-over.under
## change the termnial edges in the phylo object
pycno_tree$edge.length[tip.ids]<-corrected.terminal.edges
tip.heights2<-adephylo::distRoot(pycno_tree)
table(tip.heights2)
is.ultrametric(pycno_tree)

to_drop<-c("Peripatopsis","Scutigera","Strigamia","Craterostigmus","Parhyale","Drosophila","Tribolium",
           "Eremobates","Siro","Pachylicus","Trogulus","Phalangium","Leiobunum","Ricinoides","Limulus","Centruroides",
           "Bothriurus","Liphistius","Leucauge","Damon","Mastigoproctus")
pycno_tree<-drop.tip(pycno_tree, to_drop, trim.internal = T, subtree=F, root.edge = 0)
plot.phylo(pycno_tree, show.tip.label = T, type="phylogram",
           cex=0.4, label.offset=3)
is.ultrametric(pycno_tree)

#just do it with phytools ;)
ultra_pycno_tree<-force.ultrametric(pycno_tree, method="extend")
table(distRoot(pycno_tree))

min(pycno_tree$edge.length)
is.binary(pycno_tree)

#### occurrence data
arthro<-file.choose()
arthropoda<-read_sf("arthropoda_land.shp",
                    layer="arthropoda_land")
arthropoda <- cbind(arthropoda,st_coordinates(arthropoda))#x and Y only
arthropoda_df<-st_drop_geometry(arthropoda)
arthropoda$species<-gsub(" ", "_", arthropoda$species)

#arthropoda<-read_sf(dsn="D:\\Post-doc Antarctica\\R AQ\\pts_phylum\\arthropoda_land.shp",
 #                   layer="arthropoda_land")
#arthropoda<-st_drop_geometry(arthropoda)
arthropoda_worms<-file.choose()
arthropoda_worms<-read.csv(arthropoda_worms, dec=".", sep=";", header=T)
arthropoda_worms<-arthropoda_worms[,c(2,7)]
arthropoda_worms<-plyr::rename(arthropoda_worms, c("scientificname"="species", "scientificname.1"="newscientificname"))
arthropoda_worms$newscientificname<-gsub(" ", "_", arthropoda_worms$newscientificname)
arthropoda_worms$species<-gsub(" ", "_", arthropoda_worms$species)

setdiff(arthropoda$species, arthropoda_worms$species)
length((unique(arthropoda$species)))
plot_sf(st_geometry(arthropoda))
length(unique(arthropoda_worms$species))
arthropoda<-subset(arthropoda, !is.na(species))
#test to merge species new names (after synonymy)
#merge new names with data without points on land
arthropoda<-merge(arthropoda,arthropoda_worms,by="species",all.x=TRUE, all.y=FALSE)

#pycnos_baq<-merge(inside_pycnos_baq,arthropoda_worms,by="species",all.x=TRUE, all.y=FALSE)
#length(unique(pycnos_baq$newscientificname))

#arthro_names_phy<-unique(arthropoda_worms$newscientificname)

pycnos<-subset(arthropoda, class=="Pycnogonida")
length(unique(pycnos$newscientificname))

#dupl<-pycno_tree$tip.label[duplicated(pycno_tree$tip.label)]
pycnogonida_phylo_names<-as.data.frame(pycno_tree$tip.label)
names(pycnogonida_phylo_names)<-"species"
pycno_phy_names<-tidyr::separate(pycnogonida_phylo_names, species, c("genus","species","other"), 
                                 sep = "_", extra = "merge", fill = "right")
pycno_phy_names$spp<-paste(pycno_phy_names$genus,pycno_phy_names$species, sep="_")
pycnogonida_phylo_names<-cbind(pycnogonida_phylo_names, pycno_phy_names)
write.csv2(pycnogonida_phylo_names, file="pycnos_phylo_names_change.csv")

#rename tip labels with corrected names (delete some Gen_sp, will limit some polytomies)
t<-file.choose()
pycnogonida_phylo_names<-read.table(t, dec=".", sep=";", header=T)
pycno_tree$tip.label<-pycnogonida_phylo_names$spp
#or
#### PHYLOGENETIC IMPUTATIONS ###########
set.seed(42)
library(pez)
setdiff(pycnos$newscientificname, pycno_tree$tip.label)
setdiff(pycno_tree$tip.label,pycnos$newscientificname)

pycnos_to_insert<-setdiff(pycnos$newscientificname, pycno_tree$tip.label)
pycno_tree_merge<-congeneric.merge(pycno_tree, pycnos_to_insert, split="_")
plotTree(pycno_tree_merge)#added 231 species
plot.phylo(pycno_tree_merge, show.tip.label = T, type="phylogram",
                    cex=1, label.offset=1, show.node.label = T)

count_pycnos<-pycnos3031%>%
  count(newscientificname)

not_inserted_pycnos<-setdiff(pycnos_to_insert, pycno_tree_merge$tip.label)
setdiff(pycno_grid$newscientificname, pycno_tree_merge$tip.label)
# on laisse tomber - verifier apres filtrage des occ=1 et à l'interieur de PF
write.tree(pycno_tree_merge, file="pycno_tree_merge.nex")

#pycnos_baq <- pycnos_baq %>% mutate_all(na_if,"")
#pycnos_bas_insert<-subset(pycnos_baq, !newscientificname %in% pycno_tree_merge$tip.label)
#pycnos_bas_insert<-unique(pycnos_bas_insert$species)
### change some species names 
#pycno_tree_merge2<-congeneric.merge(pycno_tree_merge, pycnos_bas_insert)
#setdiff(pycnos_bas_insert$species, pycno_tree_merge2$tip.label)

diff<-pycnos_bas_insert%>%
  select(species, newscientificname)%>%
  st_drop_geometry()
diff$newscientificname<-diff$species
pycnos_baq<-merge(pycnos_baq, diff, by="species", all.x=T, all.y=F)
pycnos_baq<-pycnos_baq %>% 
  mutate(newscientificname.x = coalesce(newscientificname.x,newscientificname.y))%>%
  select(!newscientificname.y)
names(pycnos_baq)[names(pycnos_baq) == 'newscientificname.x'] <- 'newscientificname'
### new species names mixed in newscientificname

# a corriger les noms trouvés par le 
#write.csv2(pycnos, file="pycnos_spp_names.csv")

## check for duplicata in new names
pycnos_df<-st_drop_geometry(pycnos)
pycnos_nodup<-pycnos_df%>%
  distinct(X,Y,newscientificname,.keep_all = TRUE)%>%
  filter(!is.na(newscientificname))%>%
  filter(!newscientificname=="")
length(unique(pycnos_nodup$newscientificname))#280
## OR ##
#pycnos<-pycnos[!duplicated(pycnos[,c("newscientificname")]),]

pycnos3031 = st_as_sf(pycnos_nodup, coords = c("X", "Y"), crs = 3031)
st_crs(pycnos3031)
plot(st_geometry(pycnos3031))
pycnos3031$newscientificname<-gsub(" ", "_", pycnos3031$newscientificname)
pycnos3031<-subset(pycnos3031, !is.na(newscientificname))
length(unique(pycnos3031$newscientificname))#

100 * sum(is.na(pycnos3031$year), na.rm = TRUE) / length(pycnos3031$year)


st_write(pycnos3031, "pycnos3031.shp", layer="pycnos3031", driver="ESRI Shapefile")
pycnos3031<-st_read(dsn="D:\\as_proj\\pycnos3031.shp",
                    layer="pycnos3031")
plot(st_geometry(pycnos3031))
pycnos3031<-subset(pycnos3031, !is.na(pycnos3031$newscientificname))

prj_laea <- "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 

pycnos_laea<-st_transform(pycnos3031, crs=prj_laea)

pycno_grid<-st_intersection(pycnos_laea, grid_aq_3_PF)## 
plot(st_geometry(pycno_grid))
plot(st_geometry(grid_aq_3_PF), add=T)
length(unique(pycno_grid$newscientificname))#
length(unique(pycno_grid$grids))#

pycnogonida<-as.data.frame(unique(pycno_grid$newscientificname))## 
names(pycnogonida)<-"species"
pycnogonida$species<-gsub(" ",  "_", pycnogonida$species)

pycno_grid <- cbind(pycno_grid,st_coordinates(pycno_grid))#x and Y only
pycno_grid<-pycno_grid%>%
  distinct(X,Y,newscientificname,.keep_all = TRUE)

pycno_sampl_eff<-pycno_grid%>%
  distinct(X,Y,.keep_all = TRUE)
count_stations_pycno<-pycno_sampl_eff%>%
  group_by(grids, .add=T) %>%
  summarise(count = dplyr::n_distinct(X,Y))

#count species frequency in total
count_pycnos<-pycno_grid%>%
  count(newscientificname)
count_pycnos<-st_drop_geometry(count_pycnos)
pycno_grid<-merge(pycno_grid, count_pycnos, by="newscientificname", all.x=T, all.y=F)
rm(count_pycnos)
table(pycno_grid$n) # 29
table(pycno_grid$family)
pycnos_occ1<-count_pycnos%>%filter(n==1)
write.csv2(pycnos_occ1, file="pycnos_occ1.csv")

pycno_grid_no_occ_1<-pycno_grid%>%
  filter(!n ==1)
length(unique(pycno_grid_no_occ_1$newscientificname))#
length(unique(pycno_grid$newscientificname))#

#clculate SR in grid cells and delete SR =1
SR_pycnos<-pycno_grid%>%
  group_by(grids, .add=T) %>%
  st_drop_geometry()%>%
  summarise(count = dplyr::n_distinct(newscientificname))
pycno_grid<-merge(pycno_grid, SR_pycnos, by = "grids", all.x=T, all.y=F)
pycno_grid<-pycno_grid%>%
  filter(!count == 1)
length(unique(pycno_grid$newscientificname))# 223
length(unique(pycno_grid$grids))

100 * sum(is.na(pycno_grid$year), na.rm = TRUE) / length(pycno_grid$year)

## create community matrix with phyloregions (special format for phyloreg)
library(phyloregion)
set.seed(34)
# faster way to create a community matrix (then has to pass to classic data frame for picante)
#phyloregions::PD function does not calculate SR in the same time
sparse_grid3_pycnos <- long2sparse(pycno_grid, grid="grids", species="newscientificname")
dim(sparse_grid3_pycnos)

diff<-setdiff(colnames(sparse_grid3_pycnos),pycno_tree_merge$tip.label)
diff

## add new genera on a tree
#tree_test <- makeNodeLabel(pycno_tree_merge, "u", 
 #                          nodeList = list(Pallenopsis = "Pallenopsis"))

plot.phylo(pycno_tree_merge, show.tip.label = T, type="phylogram",
           cex=0.4, label.offset=3,show.node.label =T)
plotTree(pycno_tree_merge,node.numbers=T, fsize=0.5)

pycno_tree_merge<-bind.tip(pycno_tree_merge,"Bathypallenopsis_longiseta",where=350, 
                           position=0.5*pycno_tree_merge$edge.length[which(pycno_tree_merge$edge[,2]==350)])
# use to dd species to a specific node number
tiff("pycno_tree_add.tiff",width=10,height=8,res=600, units="in", compression ="lzw")
plotTree(pycno_tree_merge,node.numbers=T, fsize=0.2, ftype="i")
dev.off()

pycno_tree_merge<-bind.tip(pycno_tree_merge,"Cilunculus_cactoides",where=329, 
                           position=0.5*pycno_tree_merge$edge.length[which(pycno_tree_merge$edge[,2]==329)])
pycno_tree_merge<-congeneric.merge(pycno_tree_merge, c("Cilunculus_kravcovi","Cilunculus_acanthus"), split="_")

pycno_grid_no_occ_1%>%
  filter(newscientificname == "Sexanymphon_mirabilis")%>%
  count(newscientificname)
#Bathypallenopsis_longiseta #O
#Cheilopallene_gigantea#2 NO!
#Cilunculus_acanthus#occ=9
#Cilunculus_cactoides#occ=54
#Cilunculus_kravcovi occ=4
#Dodecolopoda_mawsoni occ=23
#Dromedopycnon_acanthus#occ=2 NO!
#Heteronymphon_exiguum occ=25 
# Oropallene_dimorpha occ=3 NO!
# Sexanymphon_mirabilis occ=8
pycno_tree_merge<-bind.tip(pycno_tree_merge,"Dodecolopoda_mawsoni",where=345, 
                           position=0.5*pycno_tree_merge$edge.length[which(pycno_tree_merge$edge[,2]==345)])
# chek the spatial distribution of specific species
ggplot()+geom_sf(data=subset(pycno_grid_no_occ_1,newscientificname == "Sexanymphon_mirabilis"),
                             shape = 16, size=2, alpha = 0.6, color = "purple4" )+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red2", size=6,linetype = "dashed")

pycno_tree_merge<-bind.tip(pycno_tree_merge,"Heteronymphon_exiguum",where=389, 
                           position=0.5*pycno_tree_merge$edge.length[which(pycno_tree_merge$edge[,2]==389)])

pycno_tree_merge<-bind.tip(pycno_tree_merge,"Sexanymphon_mirabilis",where=390, 
                           position=0.5*pycno_tree_merge$edge.length[which(pycno_tree_merge$edge[,2]==390)])


# check if all species from a community data are present on a tree
com.data.pycno<-match_phylo_comm(pycno_tree_merge, sparse_grid3_pycnos, delete_empty_rows = F)
dim(com.data.pycno$comm)# some species gone
com.data.pycno$phy

diff<-setdiff(colnames(sparse_grid3_pycnos),com.data.pycno$phy$tip.label)
diff
pycnos_occ1<-read.csv2("pycnos_occ1.csv")
pycnos_occ1<-pycnos_occ1%>%filter(keep == "no")

pycno_grid<-pycno_grid%>%
  filter(!newscientificname %in% pycnos_occ1$newscientificname)
length(unique(pycno_grid$newscientificname))# 214

write.tree(pycno_tree_merge, file="pycno_tree_merge_300723.nex")
## recent tree version ----
pycno_tree_merge<-read.tree("pycno_tree_merge_300723.nex")

write.csv2(pycno_grid, file="pycno_grid_10823.csv")

# calculate phylodiversity ----
sparse_grid3_pycnos <- long2sparse(pycno_grid, grid="grids", species="newscientificname")
dim(sparse_grid3_pycnos)
diff<-setdiff(colnames(sparse_grid3_pycnos), pycno_tree_merge$tip.label)
diff # 4 species can not be added

# check if all species from a community data are present on a tree

com.data.pycno<-match_phylo_comm(pycno_tree_merge, sparse_grid3_pycnos, delete_empty_rows = F)
dim(com.data.pycno$comm)# some species gone
com.data.pycno$phy$tip.label
comm_pycnos<-com.data.pycno$comm
#comm_pycnos<-comm_pycnos[!row.names(comm_pycnos) %in% "cell_3.1337",]
dim(comm_pycnos) # 
#comm_pycnos<-as.data.frame(as.matrix(comm_pycnos))

ses_PD_pycnos<-PD_ses(comm_pycnos, com.data.pycno$phy, model="tipshuffle", reps=999)
PD_pycnos<-merge(grid_aq_3_PF, ses_PD_pycnos, by="grids", all.x =T, all.y=F)
PD_pycnos <- merge(PD_pycnos, data.frame(grids=names(phylo_endemism_pycnos), PEW=phylo_endemism_pycnos), by="grids")

PD_pycnos<-PD_pycnos%>%
  filter(!is.na(PD_obs))#enleve les cellules de grid_as=q non utilisées
setdiff(rownames(comm_pycnos),PD_pycnos$grids)
## delete two celles with strangely zero, so we have 140 grid cells and not 142: "cell_3.1337" "cell_3.1355"         
# phylogenetic endemism
PE<-phylo_endemism(com.data.pycno$comm, com.data.pycno$phy)
WE<-weighted_endemism(com.data.pycno$comm)
PD_pycnos <- merge(PD_pycnos, data.frame(grids=names(PE), PE=PE), by="grids")
PD_pycnos <- merge(PD_pycnos, data.frame(grids=names(WE), WE=WE), by="grids")
PD_pycnos$std_PE <- as.vector(scale(resid(loess(PE ~ WE, data=PD_pycnos))))
## richness in PD_ses is actually a frequency of species !!! cause comm data is not 1/0

############stop
library(mgcv)
is.rooted(com.data.pycno$phy)
is.rooted(prunedTree)
plot(com.data.pycno$phy)

pd_pycno_pez<-as.data.frame(picante::pd(as.matrix(comm_pycnos), com.data.pycno$phy, include.root =F))#single-sp comm will not be cal
#prunedTree <- prune.sample(comm_pycnos,com.data.pycno$phy)
#pd_pycno_pez<-as.data.frame(picante::pd(as.matrix(comm_pycnos), prunedTree, include.root = T))#single-sp comm will not be cal

pd_pycno_pez$grids<-rownames(pd_pycno_pez)
PD_pycnos<-merge(PD_pycnos, pd_pycno_pez, by="grids", all.x =T, all.y=F)
# delete celles with one species ! 33 here
#PD_pycnos_withSR1<-PD_pycnos
PD_pycnos<-PD_pycnos%>%
  filter(!is.na(PD))
# delete these  cells from comm data
#comm_pycnos<-comm_pycnos[(rownames(comm_pycnos) %in% PD_pycnos$grids),]
#comm_pycnos<-comm_pycnos[,!(colSums(comm_pycnos) == 0)]

PD_pycnos$std_resid_loess<-as.vector(scale(resid(loess(PD~SR,data=PD_pycnos))))
#PD_pycnos$std_resid_gam <- as.vector(scale(resid(gam(PD ~ s(SR), data=PD_pycnos))))
#PD_pycnos$std_zscore<-as.vector(scale(PD_pycnos$zscore))
library(ggplot2)
PD_pycnos%>%ggplot(aes(x=SR, y=richness))+geom_point()

## calc MPD and MNTD ----
# use of comm_pycnos lready without 1 species cells
phy.dist_pycnos<-cophenetic(com.data.pycno$phy)
pycno_mpd<-as.data.frame(picante::mpd(comm_pycnos, phy.dist_pycnos))## ?? comm_comm 
names(pycno_mpd)<-"mpd"
pycno_mpd$grids<-PD_pycnos$grids
head (pycno_mpd)
pycno_mntd<-as.data.frame(picante::mntd(comm_pycnos, phy.dist_pycnos))
names(pycno_mntd)<-"mntd"
pycno_mntd$grids<-PD_pycnos$grids
head(pycno_mntd)
PD_pycnos<-merge(PD_pycnos, pycno_mpd, by="grids", all.x=T, all.y=F)
PD_pycnos<-merge(PD_pycnos, pycno_mntd, by="grids", all.x=T, all.y=F)

PD_pycnos$std_mpd <- as.vector(scale(resid(loess(mpd ~ SR, data=PD_pycnos))))
PD_pycnos$std_mntd <- as.vector(scale(resid(loess(mntd ~ SR, data=PD_pycnos))))


ses_mpd_pycnos<-ses.mpd(comm_pycnos, phy.dist_pycnos, null.model="taxa.labels", runs=999)
head(ses_mpd_pycnos)
ses_mpd_pycnos$grids<-row.names(ses_mpd_pycnos)
ses_pycnos<-merge(grid_aq_3_PF, ses_mpd_pycnos, by="grids", all.x =T, all.y=F)

ses_mntd_pycnos<-ses.mntd(comm_pycnos, phy.dist_pycnos, null.model="taxa.labels", runs=999)
head(ses_mntd_pycnos)
summary(ses_mntd_pycnos)
ses_mntd_pycnos$grids<-row.names(ses_mntd_pycnos)
ses_pycnos<-merge(ses_pycnos, ses_mntd_pycnos, by="grids", all.x =T, all.y=F)

PD_pycnos<-merge(PD_pycnos, st_drop_geometry(ses_pycnos), by="grids", all.x=T, all.y=F)

PD_pycnos_noXY<-PD_pycnos%>%
  st_drop_geometry(PD_pycnos)%>%
  select(grids, PD, std_resid_loess, SR, mpd, mntd, 
         mpd.obs.z,mntd.obs.z,std_mpd, std_mntd,
         PE, WE, std_PE,)
pairs(PD_pycnos_noXY[,c(2:13)], lower.panel = panel.smooth,  upper.panel = panel.cor, cex = 1, cex.labels = 1)

ED_pycnos<-picante::evol.distinct(com.data.pycno$phy)
names(ED_pycnos)<-c("newscientificname", "ED_tot")
ED_pycnos<-merge(ED_pycnos, count_pycnos, by = "newscientificname", all.x=T, all.y=F)
#glm(ED_tot~n, data=ED_pycnos)
#ggplot(ED_pycnos, aes(ED_tot, n))+geom_point()+geom_smooth(method='lm')
# ps de relation !

pycno_grid<-merge(pycno_grid, ED_pycnos, by="newscientificname", all.x=T, all.y=F)
cor.test(pycno_grid$n, pycno_grid$ED_tot)## 
ggplot(pycno_grid, aes(ED_tot, n))+geom_point()+geom_smooth(method='lm')
# attention - there are still species not added to the phylogeny 
# species with highest ED are those i've added...cause alone to their genera
data <- ED_pycnos %>% arrange(desc(ED_pycnos$ED_tot))
obs <- nrow(data) 
topED5_pycno5<-data %>% filter(row_number() < obs * 0.05)
names(topED5_pycno5)<-c("newscientificname", "topED_5")
pycno_grid<-merge(pycno_grid, topED5_pycno5, by="newscientificname", all.x=T, all.y=F)

## spatial congruence and overlap between PD SR horspots and MPAs ----
centroids_sf<-st_centroid(grid_aq_3_PF)
centoids_pd_pycno<-merge(centroids_sf, PD_pycnos_noXY, by="grids", all.x =T, all.y=F)
centoids_pd_pycno<-centoids_pd_pycno%>%filter(!is.na(PD))
st_crs(centoids_pd_pycno)
plot(st_geometry(centoids_pd_pycno))

summary(mpa_eez)
mpa_eez$mpa<-gsub(" ",  "_", mpa_eez$mpa)
mpa_eez<-st_transform(mpa_eez, crs=prj_laea)
centroids_mpa_eez<-st_intersection(centroids_sf, mpa_eez)
centroids_mpa_proposed<-st_intersection(centroids_sf, mpa_proposed)
centroids_mpa_proposed<-centroids_mpa_proposed%>%select(grids, name, geometry)
plot(st_geometry(centroids_mpa_proposed))  
plot(st_geometry(mpa_proposed), add=T)
plot(st_geometry(centroids_mpa_eez))
plot(st_geometry(mpa_eez), add=T)

centroids_mpa_pd_pycno<-merge(centoids_pd_pycno, st_drop_geometry(centroids_mpa_eez), 
                              by="grids", all.x =T, all.y=F )
centroids_mpa_pd_pycno<-merge(centroids_mpa_pd_pycno, st_drop_geometry(centroids_mpa_proposed), 
                              by="grids", all.x =T, all.y=F )

# define hotspots of PD and SR

## 10 percents
data <- centroids_mpa_pd_pycno %>% 
  select(std_resid_loess, grids)%>%
  arrange(desc(std_resid_loess))
obs <- nrow(data) 
hotspot10_resid_loess_pycno<-data %>% filter(row_number() < obs * 0.1)# or *0.05
names(hotspot10_resid_loess_pycno)<-c("hotspot10_resid_loess", "grids", "geometry")
plot(hotspot10_resid_loess_pycno["hotspot10_resid_loess"])

centroids_mpa_pd_pycno<-merge(centroids_mpa_pd_pycno, st_drop_geometry(hotspot10_resid_loess_pycno), 
                              by="grids", all.x =T, all.y=F )

data <- centroids_mpa_pd_pycno %>% 
  select(std_resid_loess, grids)%>%
  arrange(std_resid_loess)
obs <- nrow(data) 
coldspot10_resid_loess_pycno<-data %>% filter(row_number() < obs * 0.1)# or *0.05
names(coldspot10_resid_loess_pycno)<-c("coldspot10_resid_loess", "grids", "geometry")
plot(coldspot10_resid_loess_pycno["coldspot10_resid_loess"])

centroids_mpa_pd_pycno<-merge(centroids_mpa_pd_pycno, st_drop_geometry(coldspot10_resid_loess_pycno), 
                              by="grids", all.x =T, all.y=F )

data <- centroids_mpa_pd_pycno %>% 
  select(SR, grids)%>%
  arrange(desc(SR))
obs <- nrow(data) 
hotspot10_SR<-data %>% filter(row_number() < obs * 0.1)
names(hotspot10_SR)<-c("hotspot10_SR", "grids", "geometry")
plot(hotspot10_SR)

centroids_mpa_pd_pycno<-merge(centroids_mpa_pd_pycno, st_drop_geometry(hotspot10_SR), 
                               by="grids", all.x =T, all.y=F )

hotspots_pycno<-centroids_mpa_pd_pycno%>%select(grids, hotspot10_SR, hotspot10_resid_loess,coldspot10_resid_loess)
PD_pycnos<-merge(PD_pycnos, st_drop_geometry(hotspots_pycno), by="grids", all.x=T, all.y=F)

PD_pycnos_no1sp<-subset(PD_pycnos, !is.na(PD_pycnos$mpd))# or SR == 1
PD_pycnos_no1sp$std_resid_loess_mpd<-as.vector(scale(resid(loess(mpd~SR,data=PD_pycnos_no1sp))))
PD_pycnos_no1sp$std_resid_loess_mntd<-as.vector(scale(resid(loess(mntd~SR,data=PD_pycnos_no1sp))))


## creates polytomies at each genus node, do not use
pycno_tree<-force.ultrametric(pycno_tree, method="extend")
is.ultrametric(pycno_tree)
## number of species added varies at each run

ses.mpds.pycnos<-list()
for(i in 1:10){
  print (i)
  imp.tree <- congeneric.impute(pycno_tree, pycnos_to_insert)
  c.data <- comparative.comm(imp.tree, as.matrix(comm_pycnos))
  ses.mpds.pycnos[[i]] <- .ses.mpd(c.data)$mpd.obs.z
  names(ses.mpds.pycnos[[i]])<-rownames(comm_pycnos)
  
}
## ses.mpds mean per grid cell
ses.mpds.pycnos.X<-lapply(ses.mpds.pycnos, as.matrix)
ses.mpds.pycnos.X<-do.call(cbind, ses.mpds.pycnos.X)
mean_ses_pd_pycnos<-as.data.frame(rowMeans(ses.mpds.pycnos.X, na.rm=T))
rownames(mean_ses_pd_pycnos)<-rownames(comm_pycnos)
mean_ses_pd_pycnos<-as.data.frame(mean_ses_pd_pycnos)
mean_ses_pd_pycnos$grids<-rownames(mean_ses_pd_pycnos)
names(mean_ses_pd_pycnos)<-c("mean.ses.pd","grids")
write.csv(mean_ses_pd_pycnos, file="mean_ses_pd_pycnos.csv")
## loess with MPD 
loess.pycno<-list()
for(i in 1:2){
  print (i)
  imp.tree <- congeneric.impute(pycno_tree, pycnos_to_insert)
  c.data <- comparative.comm(imp.tree, comm_pycnos)
  MPD_obs <- .ses.mpd(c.data)
  residuals_loess<-loess(mpd.obs~ntaxa, data=MPD_obs, na.action=na.exclude)$residuals
  loess.pycno[[i]]<-scale(residuals_loess)
}
loess.pycno[[1]]
loess.pycno.table<-lapply(loess.pycno, rbind)
class(loess.pycno.table[[1]])
loess.pycno.table<-lapply(loess.pycno.table, as.data.frame)
library(data.table)
loess_PD_pycno<-rbindlist(loess.pycno.table, use.names=T, fill=T)
mean_loess_PD_pycno<-as.data.frame(colMeans(loess_PD_pycno, na.rm=T))
### when data have the same length in grid cell numbers
#mean_loess_PD_ophi<-as.data.frame(rowMeans(loess_PD_ophi, na.rm=T))
mean_loess_PD_pycno$grids<-rownames(mean_loess_PD_pycno)
names(mean_loess_PD_pycno)<-c("mean.loess","grids")
write.csv(mean_loess_PD_pycno, file="mean_loess_PD_pycno.csv")

### Loess regression with PD vlues
set.seed(45)
pd.pycno.loess<-list()
for(i in 1:1000){
  print (i)
  imp.tree <- congeneric.impute(pycno_tree, pycnos_to_insert)
  c.data<-match_phylo_comm(imp.tree, sparse_grid3_pycnos, delete_empty_rows = F)
  X<-pd(as.matrix(c.data$comm), c.data$phy, include.root = T)
  X$grids<-rownames(c.data$comm)
  Y<-subset(X, !PD==0)
  resids_loess <- resid(loess(Y$PD ~ Y$SR))
  pd.pycno.loess[[i]]<-cbind(resids_loess,Y)
}

class(pd.pycno.loess[[1]])
setdiff(pd.pycno.loess[[1]]$grids, pd.pycno.loess[[4]]$grids)

class(loess.pycno.table[[1]])
library(data.table)
library(dplyr)
pd.pycno.loess<-lapply(pd.pycno.loess,as.data.table)
pd.pycno.loess<-rbindlist(pd.pycno.loess, use.names=T, fill=T)
mean_loess_PD_pycno<- pd.pycno.loess%>%group_by(grids)%>%summarise_all("mean")
mean_loess_PD_pycno$scaled_loess_resid<-scale(mean_loess_PD_pycno$resids_loess)
write.csv(mean_loess_PD_pycno, file="mean_loess_PD_pycno.csv")
t<-file.choose()
mean_loess_PD_pycno<-read.table(t, dec=".", sep=",", header=T)
##ATTENTION  IL FAUT AUSSI LA RICHESSE EN MOYENNE PAR PLOT IMPUTE!!!
mean_loess_PD_pycno<-merge(grid_aq_3, mean_loess_PD_pycno, by="grids", all.x=T, all.y=F)

## some tests
h <- pd(c.data$comm, c.data$phy, include.root=T)
plot(h$PD, h$SR)
pd.ivs.loess <- resid(loess(h$PD ~ rowSums(c.data$comm)))
pd.ivs.lm <- unname(resid(lm(h$PD ~ rowSums(c.data$comm))))
cor.test(scale(pd.ivs.loess), scale(pd.ivs.lm), method="pearson")
par(mfrow=c(2,1))

plot(h$PD, predict(loess(h$PD ~ rowSums(c.data$comm))))
plot(h$PD, predict(lm(h$PD ~ rowSums(c.data$comm))))

#########################################
# compute phylogenetic beta diversity btw cells
set.seed(42)
phylo.beta.pycno<-list()
for(i in 1:1000){
  print (i)
  imp.tree <- congeneric.impute(pycno_tree, pycnos_to_insert)
  imp.tree<-as.phylo(imp.tree)
  c.data<-match_phylo_comm(imp.tree, sparse_grid3_pycnos, delete_empty_rows = F)
  phylo.beta.pycno[[i]] <- phylobeta(c.data$comm, c.data$phy)$phylo.beta.sor
}

dim(phylo.beta.pycno[[1]])
sparse_grid3_pycnos<-as.matrix(sparse_grid3_pycnos)

## some dist matices have NA (kept empty rows of unadded species)
## have to manage it to get mean of these dist list
#mean_phylo_beta_pycno<-mult_dist_average(phylo.beta.pycno)
## do the same
mean_phylo_beta_pycno<-mean_dist(phylo.beta.pycno)

#Y <- do.call(cbind, phylo.beta.pycno)
#Y <- array(Y, dim=c(dim(phylo.beta.pycno[[1]]), length(phylo.beta.pycno)))

#mean_phylo_beta_pycno<-Reduce(`+`, phylo.beta.pycno) / length(phylo.beta.pycno)

#names(mean_phylo_beta_pycno)<-rownames(comm_pycnos)
#mean_phylo_beta_pycno<-as.dist(mean_phylo_beta_pycno)
#dim(mean_phylo_beta_pycno)


mean_phylo_beta_pycno = as.matrix(mean_phylo_beta_pycno)
mean_phylo_beta_pycno = mean_phylo_beta_pycno[rowSums(is.na(mean_phylo_beta_pycno)) == 0,
                                              colSums(is.na(mean_phylo_beta_pycno)) == 0, drop = FALSE]
mean_phylo_beta_pycno<-as.dist(mean_phylo_beta_pycno)
dim(mean_phylo_beta_pycno)## 

names(mean_phylo_beta_antho)
### le nombre de cellules change car probablement dans les ajouts d'especes random,
## certaines cellules ne contiennent que des especes qui n'ont pas été ajouté à un "run" 
## du coup les cellules sont vides et elimnées de sparse gris des pycnos
## faire la moyenne que des cell grids qui ont le meme nom
## calculer la moyenne de phylobeta dans toutes les cellules en indiquant le max de cellules
## comme dans sparse phylo comm
## !!! in phyloregions::mean_dist all matrices should be of the same dimension

sparse_grid3_pycnos <- long2sparse(pycno_grid, grid="grids", species="newscientificname")
com.data.pycnos<-match_phylo_comm(pycno_tree_merge, sparse_grid3_pycnos, delete_empty_rows = F)
setdiff(com.data.pycnos$phy$tip.label,colnames(sparse_grid3_pycnos))

phylo_beta_pycnos <- phylobeta(comm_pycnos, com.data.pycno$phy)
#select the less distorting clustering method
select_linkage(phylo_beta_pycnos[[3]])
#select optimal number of clusters with selected method
optim<-optimal_phyloregion(phylo_beta_pycnos[[3]], k=15, method="average")
plot(optim$df$k, optim$df$ev)
#
##phyloregionalisation
y <- phyloregion(phylo_beta_pycnos[[3]], shp=grid_aq_3_sp, k=7, method="average")
summary(y)
phylo_nmds<-y$NMDS
pycno_phyloreg_sf<-y$shp
plot(pycno_phyloreg_sf)
pycno_phyloreg_sf<-st_as_sf(pycno_phyloreg_sf)

class(pycno_phyloreg_sf["ED"])

tiff("ED_K7_phylobeta_pycnos_209sp_bis.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+geom_sf(data=pycno_phyloreg_sf, aes(fill = as.factor(as.character(COLOURS))))
plot(pycno_phyloreg_sf["ED"])
plot(Coastline,  col="white", add=T)
dev.off()

phylo_nmds<-y$NMDS
phylo_nmds<-y$NMDS$points
class(phylo_nmds)
phylo_nmds<-as.data.frame(phylo_nmds)
phylo_nmds$phyloreg<-rownames(phylo_nmds)

tiff("nmds_phyloregions_pycnos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
pycno_nmds<-ggplot()+geom_point(data=phylo_nmds, aes(x=MDS1, y=MDS2, fill=phyloreg), pch=21, size=6)+
  scale_fill_brewer(palette = "Dark2") +   # customized color palette
  labs(fill = "")+
  geom_text(data=phylo_nmds, aes(x=MDS1, y=MDS2, label=phyloreg), vjust=0.3, size=4)+
  theme_minimal()+
  theme(legend.position="none")
dev.off()


y$membership
y$region.df
PD_pycnos<-merge(PD_pycnos, y$region.df, by="grids", all.x=T, all.y=F)
names(PD_pycnos)
#subset(pycno_grid, grids == "cell_3.975")# cluster 10

#####################################
#some phylogenetic trees stuff ######
is.binary(pycno_tree)
plot(pycno_tree_impute, cex=0.5)
min(pycno_tree$edge.length)
max(pycno_tree$edge.length)
## created bigger branch lengths
pycno_tree_congener_di<-multi2di(pycno_tree_congener, random=T, equiprob=T)
is.binary(pycno_tree_congener_di)
plot(pycno_tree_congener_di, cex=0.8)
axisPhylo() # Shows the scale of the branch lengths in millions of years.

## with add species to genus
library(phytools)
pycno_tree_added<-add.species.to.genus(pycno_tree,"Tanystylum_styligerum" , where="random")#with where=root created polytomy
plot(pycno_tree_added, cex=0.5)

pycno_tree_added<-add.species.to.genus(ultra_pycno_tree,pycnogonida)

for(i in 1:length(pycnogonida_test)){
  ultra_pycno_tree<-add.species.to.genus(ultra_pycno_tree,pycnogonida[i],where="random")
}
#bugs...if too many species
pycno_tree_added<-add.species.to.genus(pycno_tree_added,"Nymphon_articulare")

# change genera to ancient Pallenopsis and add as polytomy 
#"Bathypallenopsis_antipoda","Bathypallenopsis_longirostris","Bathypallenopsis_longiseta"
#"Boehmia_dubia" - nomen dubium
# Family of Callipallenidae, add as a polytomy/random family addition
#to one of existing genera of this family:
# Callipallene, Austropallene, Stylopallene (1 esp), Anoropallene (1 esp)
#"Cheilopallene_gigantea", "Cheilopallene_trappa","Oropallene_dimorpha",
#"Oropallene_dolichodera", "Oropallene_metacaula", "Seguapallene_insignatus"
# Ammotheidae
#"Cilunculus_acanthus","Cilunculus_cactoides",         
#"Cilunculus_kravcovi","Cilunculus_sewelli","Cilunculus_spinicristus",
#"Dromedopycnon_acanthus"
#Colossendeidae 
#"Dodecolopoda_mawsoni","Hedgpethia_eleommata"
# Nymphonidae
#"Heteronymphon_exiguum", "Sexanymphon_mirabilis"

#change in occurrence data
#"Meridionale_ambigua" to "Pallenella_ambigua"
pycnos3031$newscientificname<-sub("Meridionale_ambigua", "Pallenella_ambigua", pycnos3031$newscientificname)

## adding single species to a phylogeny ----
plot(ultra_pycno_tree,no.margin=TRUE,edge.width=2,cex=0.7, node.number=T,label.offset=0.1)
nodelabels(frame="none",adj=c(1.1,-0.4))
tiplabels(frame="none", cex=0.6)
#nodelabels(text=1:pycno_tree$Nnode,node=1:pycno_tree$Nnode+Ntip(pycno_tree),frame="none",adj=c(1.1,-0.4))
## same but ugly
plotTree(pycno_tree,no.margin=TRUE,edge.width=1,cex=0.5, node.numbers=T)
# with ape::bind.tree needs a "phylo" object with edge lengths to add
branching.times(pycno_tree)
pycno_tree$edge.length[which(pycno_tree$edge[,2]==88)]
1.901771
pycno_tree$edge.length[which(pycno_tree$edge[,2]==198)]
0.235502
pycno_tree$edge.length[which(pycno_tree$edge[,2]==198)]/2
0.117751
1.901771+0.203616+0.117751
2.223138
tip<-list(edge=matrix(c(2,1),1,2),
          tip.label="Heteronymphon_exiguum",
          edge.length=2.223138,
          Nnode=1)
class(tip)<-"phylo"
## OR
tip <- rtree(1, tip.label="new_tip");
# for Nymphonidae
tree_test<-bind.tree(pycno_tree, tip,where=198,
                    position=0.5* pycno_tree$edge.length[which(pycno_tree$edge[,2]==198)])
## or with phytools
tree_test<-bind.tip(pycno_tree, "tip",where=198,
                     position=0.5* pycno_tree$edge.length[which(pycno_tree$edge[,2]==198)])

## or for one by one additions
node<-198
tt<-splitTree(pycno_tree,split=list(node=node,
                              bp=pycno_tree$edge.length[which(pycno_tree$edge[,2]==node)]))
tt[[2]]<-add.random(tt[[2]],tips=c("spp1", "spp2", "spp3"))
new.tree<-paste.tree(tt[[1]],tt[[2]])
plotTree(new.tree)
is.ultrametric(new.tree)# yes, as pycno tree is 

plotTree(tree_test,no.margin=TRUE,edge.width=2,cex=0.7, node.numbers=F)
nodelabels(frame="none",adj=c(1.1,-0.4))
#mrca_pycno<-as.data.frame(mrca(pycno_tree, full=T))

tree<-file.choose()
pycno_tree<-read.nexus(tree)
pycno_tree_anna2<-read.newick(tree)
plot(pycno_tree_anna, cex=0.8)

intersect(pycnogonida$species, pycno_tree$tip.label)# 29/281!!! quasi nothing is in this tree!
pycnogonida<-pycnogonida$species

#########################
## from open tree of life ##----
##########################
pycno_ott_id <- rotl::tnrs_match_names("Arthropoda")$ott_id
#OR
#pycnos2<-pycnos[c(1:10),c(2,3,4,5,14)]
pycno_spp_otol<-tnrs_match_names(names = pycnos3$newscientificname, context_name = "All life")
pycno_spp_otol2 = filter(pycno_spp_otol, !is.na(unique_name))
pycno_spp_otol2 = unique(pycno_spp_otol2)

tree_otl = try(tol_induced_subtree(ott_ids = pycno_spp_otol2$ott_id, label_format = "name"))
no_id = c (3582562, 3583148,4963791)
tree_otl = tol_induced_subtree(ott_ids = filter(pycno_spp_otol2, !ott_id %in% no_id)$ott_id, label_format = "name")
plot.phylo(tree_otl, cex=0.7, show.node.label = TRUE, show.tip.label = FALSE)
is.rooted(tree_otl)# yes!
tree_otl$node.label
tree_otl$tip.label
tree_otl$edge

library(phangorn)
plot(tree_otl, show.tip.label=F)
nodelabels(frame="none")#only for a tree
tiplabels()
Descendants(tree_otl, 286, "tips")# get all descendant tips of internal node

# for pycnogonida, out of 275 spp timetree has only 3 
timetree_pair(tax1 = "Colossendeis macerrima", tax2 = "Achelia bituberculata")

itree_otl = timetree_phylo(tree_otl)

## get node names and id - edge table

node_labels_in_edge <- tree_otl$node.label[tree_otl$edge[,1]-Ntip(tree_otl)]
tips_nodes <- tree_otl$edge[,2]

select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}

edge_table <- data.frame(
  "parent" = tree_otl$edge[,1],
  "par.name" = sapply(tree_otl$edge[,1], select.tip.or.node, tree = tree_otl),
  "child" = tree_otl$edge[,2],
  "chi.name" = sapply(tree_otl$edge[,2], select.tip.or.node, tree = tree_otl)
)

write.tree(itree_otl$phy, file = "data/tree_otl_nodes.tre")
itree_otl$info %>% filter(!is.na(median)) %>% 
  select(node, estimated) %>% 
  write.table(file = "data/ages", quote = F, sep = "\t", row.names = F, col.names = F)

# error cause only 3 species in time tree !!!
# get newick file from time tree directly
ages_df <- data.frame(
  a = c('Arthropoda','Pantopoda'),
  b = c(583, 296)
)
library(phylocomr)
dated_pycno_test<-phylocomr::ph_bladj(ages=ages_df, phylo=tree_otl)
tips<-c("Nymphon_gracile","Colossendeis_macerrima")
tips2<-setdiff(tree_otl$tip.label,tips)

pycno_pruned<-drop.tip(tree_otl, tree_otl$tip.label[-match(tips, tree_otl$tip.label)])
pruned.tree<-drop.tip(tree_otl, setdiff(tree_otl$tip.label, tips ))

## get nodes numbers for specific species
node <- c(getMRCA(itree_otl, tip = c("Nymphon gracile","Achelia bituberculata")))# Ach. bit. n'existe plus dans otol
is.phylo(tree_otl)
library(datelife)
node_info<-tree_get_node_data(tree_otl)# does not work without branch lenghts!!


# ott_id or node_id (not both) specify the node which will be the root of the tree
pycno_otol<-tol_subtree(ott_id=pycno_ott_id, label_format = "name")#gets a huge polytomy of pycnos (at least we can see the ssp)
summary(pycno_otol)
is.rooted(pycno_otol)
plot(pycno_otol,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:pycno_otol$Nnode,node=1:tree$pycno_otol+Ntip(pycno_otol))


