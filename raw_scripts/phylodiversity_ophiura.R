library(ape)
library(phylobase)
library(phytools)
library(geiger)
#library(adiv)
library(dplyr)
library(sf)
library(adephylo)
library(phyloregion)
library(pez)
library(abind)
# Echinodermata ----
echinos<-read_sf(dsn="F:\\Post-doc Antarctica\\R AQ\\pts_phylum\\echinodermata_land.shp",
                 layer="echinodermata_land")
#echinos<-st_drop_geometry(echinos)#delete geometry column with point coordinates
echinos_worms<-file.choose()# F:\Post-doc Antarctica\R AQ\taxo_syno 
echinos_worms<-read.csv(echinos_worms, dec=".", sep=";", header=T)
echinos_worms<-echinos_worms[,c(2,7)]
echinos_worms<-plyr::rename(echinos_worms, c("scientificname"="species"))
#test to merge species new names (after synonymy)
#merge new names with data without points on land = occurences
length(unique(echinodermata$species))
length(unique(echinos_worms$species))

echinodermata<-merge(echinos,echinos_worms,by="species",all.x=TRUE, all.y=FALSE)
length(unique(echinodermata$species))
length(unique(echinodermata$newscientificname))
table(echinodermata$class)
ophiures<-subset(echinodermata, class=="Ophiuroidea")
length(unique(ophiures$newscientificname))#263
## check for duplicata in new names
ophiures <- cbind(ophiures,st_coordinates(ophiures))#x and Y only
rm(ophiures)
ophiures_df<-st_drop_geometry(ophiures)
ophiures_nodup<-ophiures_df%>%
  distinct(newscientificname,X,Y,.keep_all = TRUE)%>%
  filter(!is.na(newscientificname))
length(unique(ophiures_nodup$newscientificname))#262
## 
ophiures3031 <- st_as_sf(ophiures_nodup, coords = c("X", "Y"), crs = 3031)
st_crs(ophiures3031)
plot(st_geometry(ophiures3031))
plot(ophiures3031)
ophiures3031$newscientificname<-gsub(" ", "_", ophiures3031$newscientificname)
ophiures3031<-subset(ophiures3031, !is.na(newscientificname))
st_write(ophiures3031, "ophiures3031.shp", layer="ophiures3031", driver="ESRI Shapefile")
ophiures3031<-st_read(dsn="D:\\as_proj\\ophiures3031.shp",
                     layer="ophiures3031")

ophiura_laea<-st_transform(ophiures3031, crs=prj_laea)
plot(st_geometry(ophiura_laea))

100 * sum(is.na(ophiures3031$year), na.rm = TRUE) / length(ophiures3031$year)
# or
mean(is.na(ophiures3031$year)) * 100

## ophiures tree O'Hara et al, 2019 https://www.nature.com/articles/s41586-019-0886-z
tree<-file.choose()
ophiure_tree<-read.tree(tree)
plot.phylo(ophiure_tree, cex=0.4, show.tip.label = TRUE, type="phylogram")
is.ultrametric(ophiure_tree)
ophiures_phylo_names<-as.data.frame(ophiure_tree$tip.label)
names(ophiures_phylo_names)<-"species"##u total 1185 spp
########
#ophiura_species_aq<-as.data.frame(ophiures_phylo_names$species[grep("_AN_", 
#                                                                   as.character(ophiures_phylo_names$species))])
#names(ophiura_species_aq)<-"species"
#ophiure_tree_aq<-drop.tip(ophiure_tree, ophiure_tree$tip.label[-match(ophiura_species_aq$species, ophiure_tree$tip.label)]) # sub-tree des esp du plot en question
#plot.phylo(ophiure_tree_aq, cex=0.7, show.node.label = TRUE)
ophiura_phy_names<-tidyr::separate(ophiures_phylo_names, species, c("genus","species","other"), 
                                   sep = "_", extra = "merge", fill = "right")
ophiura_phy_names<-cbind(ophiura_phy_names,ophiures_phylo_names)
ophiura_phy_names$spp<-paste(ophiura_phy_names$genus,ophiura_phy_names$species, sep=" ")
#ophiures_phylo_names<-cbind(ophiures_phylo_names, ophiura_phy_names)
write.csv2(ophiura_phy_names, file="ophiura_phy_names.csv")
ophiura_phy_names2<-read.table("F:\\as_proj\\ophiura_phy_names.csv", sep=";", dec=".", header=T)
ophiura_phy_names2$spp<-gsub(" ", "_", ophiura_phy_names2$spp)
#rename tip labels with new names (CORRECTED by hand)
ophiure_tree$tip.label<-ophiura_phy_names2$spp
plot.phylo(ophiure_tree, cex=0.7, show.node.label = TRUE)

spp_to_delete<-subset(ophiura_phy_names2, delete=="y")
spp_to_delete<-spp_to_delete$spp

ophiure_tree<-drop.tip(ophiure_tree, spp_to_delete)#753
plot.phylo(ophiure_tree, cex=0.4, show.tip.label = TRUE, type="phylogram")
## delete duplicated species tips
## https://stackoverflow.com/questions/39403443/collapse-a-clade-by-tip-labels-while-maintaining-phylogenetic-position
## get home-made function 
## what spp are duplicated

##make changes on the tree :
ophiure_tree$tip.label<-gsub("Amphioplus_daleus", "Silax_daleus", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Amphioplus_verrilli", "Silax_verrilli", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiomusium_australe", "Ophiomusa_australe", ophiure_tree$tip.label)
# add "Ophiomusa_biporica" to the genera
ophiure_tree$tip.label<-gsub("Ophiomusium_constrictum", "Ophiomusa_constricta", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiomusium_lymani", "Ophiomusa_lymany", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiomusium_scalare",  "Ophiomusa_scalare", ophiure_tree$tip.label)
# ophioparva genera existes on the tree -> add "Ophioparva_blochi"
ophiure_tree$tip.label<-gsub("Ophiura_ambigua", "Ophiuroglypha_ambigua", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_brevispinosa", "Ophiuroglypha_brevispinosa", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_carinifera","Ophiuroglypha_carinifera", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_irrorata","Ophiuroglypha_irrorata" , ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_lymani","Ophiuroglypha_lymani" , ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_ossiculata","Ophiuroglypha_ossiculata" , ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_brevispina","Ophiuroglypha_brevispina", ophiure_tree$tip.label)
ophiure_tree$tip.label<-gsub("Ophiura_rugosa","Ophiuroglypha_rugosa" , ophiure_tree$tip.label)

write.tree(ophiure_tree, file="ophiure_tree_new.nex")

#write.tree(new_ophi_tree, file="new_ophi_tree")
#plot.phylo(new_ophi_tree, cex=0.7, show.node.label = TRUE)
### prepare data for PD analysis
#subset(ophiures3031, newscientificname %in% ophi_tree_merge$tip.label)

plot(st_geometry(ophiures3031),  key.pos = NULL, axes = TRUE)
setdiff(ophiures3031$newscientificname, new_ophi_tree$tip.label)
intersect(new_ophi_tree$tip.label,ophiures3031$newscientificname)
length(unique(ophiures3031$newscientificname))#262

## 3x3 degree antarctic ocean grid manip'
ophi_grid<-st_intersection(ophiura_laea, grid_aq_3_PF)
plot(st_geometry(ophi_grid), add=T)
plot(st_geometry(grid_aq_3_PF), add=T)
length(unique(ophi_grid$newscientificname))
length(unique(ophi_grid$grids))#262

## attention duplication de certains points d'occ qui tombent entre deux grid cells !
# les enlever:
ophi_grid <- cbind(ophi_grid,st_coordinates(ophi_grid))#x and Y only
ophi_grid<-ophi_grid%>%
  distinct(X,Y,newscientificname,.keep_all = TRUE)
length(unique(ophi_grid$newscientificname))
length(unique(ophi_grid$grids))

subset(ophi_grid, grids == "cell_3.1547")
# dans cette cellule ya la même espece deux fois avec une mini diff dans les coordonnées entre gbif et baq
## on verra plus tard, là je supprime cete cellule
ophi_grid<-ophi_grid%>%
  filter(!grids == "cell_3.1547")
# espece disparu en quesiton avec la cellule:
ophi_grid%>%filter(newscientificname == "Anophiura banzarei")

## on travaillera pas avec abondance! 
ophi_abund<-ophi_grid%>%
  dplyr::count(newscientificname)
ophi_abund<-st_drop_geometry(ophi_abund)
ophi_grid<-merge(ophi_grid, ophi_abund, by="newscientificname", all.x=T, all.y=F)
rm(ophi_abund)

# delete cells with SR = 1
SR_ophi<-ophi_grid%>%
  group_by(grids, .add=T) %>%
  st_drop_geometry()%>%
  summarise(count = dplyr::n_distinct(newscientificname))
ophi_gridPF_1<-merge(ophi_grid, SR_ophi, by = "grids", all.x=T, all.y=F)
ophi_gridPF_1<-ophi_gridPF_1%>%
  filter(!count == 1)
length(unique(ophi_gridPF_1$newscientificname))
length(unique(ophi_gridPF_1$grids))
table(ophi_gridPF_1$count)

tiff("ophiura_occ.tiff",width=6,height=6,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=ophi_grid_1, shape = 1, size=2)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red2", size=6,linetype = "dashed")+
  theme_bw()
dev.off()

spp_ophiura_to_insert<-setdiff(unique(ophi_grid$newscientificname), ophiure_tree$tip.label)#
spp_ophiura_to_insert<-gsub(" ", "_", spp_ophiura_to_insert)
ophiure_tree$tip.label<-gsub(" ", "_", ophiure_tree$tip.label)
spp_ophiura_to_insert<-subset(spp_ophiura_to_insert, !is.na(spp_ophiura_to_insert))
spp_ophiura_to_insert
write.csv2(as.data.frame(spp_ophiura_to_insert), file="spp_ophiura_to_insert.csv")

library(pez)
ophi_tree_merge <- congeneric.merge(ophiure_tree, spp_ophiura_to_insert, split="_")
setdiff(unique(ophi_grid$newscientificname), ophi_tree_merge$tip.label)#
plotTree(ophi_tree_merge,no.margin=TRUE,edge.width=2,cex=0.4, node.numbers=F)
axisPhylo()
is.ultrametric(ophi_tree_merge)
write.tree(ophi_tree_merge, file="ophiure_tree_merge.nex")

# possibilité d'ajouter comme des pycnos par numero donné et classif connue
setdiff(unique(ophi_grid$newscientificname), ophi_tree_merge$tip.label)#
# use to dd species to a specific node number
# count missing species abundance
ophi_grid%>%
  filter(newscientificname == "Ophiodaces_inanis")%>%
  count(newscientificname)
# check visually species locations
ggplot()+geom_sf(data=subset(ophi_grid,newscientificname == "Ophiodaces_inanis"),
                 shape = 16, size=2, alpha = 0.6, color = "purple4" )+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red2", size=6,linetype = "dashed")
#check the commun node according to other gengera from the same family
tiff("ophi_tree_add2.tiff",width=10,height=8,res=600, units="in", compression ="lzw")
plotTree(ophi_tree_merge, node.numbers=T, fsize=0.1, ftype="i", lwd = 0.1)
dev.off()

## create community matrix for phyloregions (special format for phyloreg)
library(phyloregion)
sparse_grid3_ophi <- long2sparse(ophi_gridPF_1, grid="grids", species="newscientificname")
dim(sparse_grid3_ophi)
class(sparse_grid3_ophi)
comm_ophi<-as.matrix(sparse_grid3_ophi)
comm_ophi<-as.data.frame(comm_ophi)
dim(comm_ophi)

# create comm matrix and phylogeny with only present species from community
com.data.ophi<-match_phylo_comm(ophi_tree_merge, sparse_grid3_ophi, delete_empty_rows = F)
com.data.ophi$phy# 168 sur l'arbre au lieu de 173 - 5 manquent par leur genre...
dim(com.data.ophi$comm)

ses_PD_ophi<-PD_ses(sparse_grid3_ophi, com.data.ophi$phy, model="tipshuffle", reps=999)
PD_ophi<-merge(grid_aq_3_PF, ses_PD_ophi, by="grids", all.x =T, all.y=F)
PD_ophi<-PD_ophi%>%
  filter(!is.na(richness))
com.data.ophi$phy#168

library(mgcv)
is.rooted(ophi_tree_merge)
c.data.ophi<-comparative.comm(ophi_tree_merge, as.matrix(comm_ophi))#130
ophi_tree_merge$node.label<-NULL
dim(c.data.ophi$comm)
### error of duplicate tips in the phylogny

comm_ophi<-as.data.frame(as.matrix(com.data.ophi$comm))
dim(comm_ophi)

pd_ophi_pez<-as.data.frame(picante::pd(comm_ophi, com.data.ophi$phy, include.root = F))
pd_ophi_pez$grids<-rownames(pd_ophi_pez)
PD_ophi<-merge(PD_ophi, pd_ophi_pez, by="grids", all.x =T, all.y=F)
PD_ophi$std_resid_loess<-as.vector(scale(resid(loess(PD_obs~SR,data=PD_ophi))))

# phylogenetic endemism
PE<-phylo_endemism(com.data.ophi$comm, com.data.ophi$phy)
WE<-weighted_endemism(com.data.ophi$comm)
PD_ophi <- merge(PD_ophi, data.frame(grids=names(PE), PE=PE), by="grids")
PD_ophi <- merge(PD_ophi, data.frame(grids=names(WE), WE=WE), by="grids")
PD_ophi$std_PE <- as.vector(scale(resid(loess(PE ~ WE, data=PD_ophi))))
## calc MPD and MNTD ----
# use of comm_pycnos lready without 1 species cells
phy.dist_ophi<-cophenetic(com.data.ophi$phy)
ophi_mpd<-as.data.frame(picante::mpd(comm_ophi, phy.dist_ophi))## ?? comm_comm 
names(ophi_mpd)<-"mpd"
ophi_mpd$grids<-PD_ophi$grids
head (ophi_mpd)

ophi_mntd<-as.data.frame(picante::mntd(comm_ophi, phy.dist_ophi))
names(ophi_mntd)<-"mntd"
ophi_mntd$grids<-PD_ophi$grids
head(ophi_mntd)
PD_ophi<-merge(PD_ophi, ophi_mpd, by="grids", all.x=T, all.y=F)
PD_ophi<-merge(PD_ophi, ophi_mntd, by="grids", all.x=T, all.y=F)

ses_mpd_ophi<-ses.mpd(comm_ophi, phy.dist_ophi, null.model="taxa.labels", runs=999)
head(ses_mpd_ophi)
ses_mpd_ophi$grids<-row.names(ses_mpd_ophi)
ses_ophi<-merge(grid_aq_3_PF, ses_mpd_ophi, by="grids", all.x =T, all.y=F)

ses_mntd_ophi<-ses.mntd(comm_ophi, phy.dist_ophi, null.model="taxa.labels", runs=999)
head(ses_mntd_ophi)
ses_mntd_ophi$grids<-row.names(ses_mntd_ophi)
ses_ophi<-merge(ses_ophi, ses_mntd_ophi, by="grids", all.x =T, all.y=F)

#pas besoin
PD_ophi$loess_mpd <- as.vector(scale(resid(loess(mpd ~ SR, data=PD_ophi))))
PD_ophi$loess_mntd <- as.vector(scale(resid(loess(mntd ~ SR, data=PD_ophi))))

PD_ophi<-merge(PD_ophi, st_drop_geometry(ses_ophi), by="grids", all.x=T, all.y=F)

write.csv(PD_ophi, file="PD_ophi_measures.csv")

PD_ophi_noXY<-PD_ophi%>%
  st_drop_geometry(PD_ophi)%>%
  select(grids, PD, std_resid_loess, SR, mpd, mntd, 
         mpd.obs.z,mntd.obs.z,loess_mpd, loess_mntd,
         PE, WE, std_PE,)
pairs(PD_ophi_noXY[,c(2:13)], lower.panel = panel.smooth,  upper.panel = panel.cor, cex = 1, cex.labels = 1)

## will use std_resid_loess
ED_ophi<-picante::evol.distinct(com.data.ophi$phy)
names(ED_ophi)<-c("newscientificname", "ED_tot")
ophi_grid<-merge(ophi_grid, ED_ophi, by="newscientificname", all.x=T, all.y=F)

data <- ED_ophi %>% arrange(desc(ED_ophi$ED_tot))
obs <- nrow(data) 
topED5_ophi2.5<-data %>% filter(row_number() < obs * 0.025)
names(topED5_ophi2.5)<-c("newscientificname", "topED_2.5")
ophi_grid<-merge(ophi_grid, topED5_ophi2.5, by="newscientificname", all.x=T, all.y=F)

tiff("PD_ophi_loess_hotspot10_ED2.5.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=PD_ophi, aes(fill=std_resid_loess), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=mpa_eez, fill=NA, color="deeppink", linewidth=0.6)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  geom_sf(data=mpa_proposed, fill=NA, color="orange", linewidth=0.6)+
  geom_sf(data=subset(PD_ophi,!is.na(hotspot10_resid_loess)), color="red", fill=NA,linewidth=1)+
  geom_sf(data=subset(ophi_grid, !is.na(topED_2.5.y)),shape = 19, size=2)+ 
   theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  coord_sf(datum = NA)
dev.off()

tiff("SR_ophi_hot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=PD_ophi, aes(fill=SR), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="D", direction = 1, name="Species\nRichness")+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, fill=NA, color="deeppink", linewidth=0.6)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  #geom_sf(data=subset(ophi_grid, !is.na(topED_2.5)),shape = 19, size=2)+ 
  geom_sf(data=mpa_proposed, fill=NA, color="orange", linewidth=0.6)+
  geom_sf(data=subset(PD_ophi,!is.na(hotspot10_SR)), color="red", fill=NA,linewidth=1)+
  theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  coord_sf(datum = NA)
dev.off()

# spp not added are probably synonims already on the tree
#synonims changed by worms while checking accepted names
# disparities btw synonimezed names on tree O'Hatra and worms
#########################################"
## phylo imputation with pez (Will Pearse!) 1000 time to get a mean branch lengths
##########################################

## MPD ###
set.seed(42)
ses.mpds.ophi2<-list()
for(i in 1:2){
  print (i)
  imp.tree <- congeneric.impute(ophiure_tree, ophiura_to_insert)
  c.data <- comparative.comm(imp.tree, comm_ophi)
  ses.mpds.ophi2[[i]] <- .ses.mpd(c.data)$mpd.obs.z
  names(ses.mpds.ophi2[[i]])<-rownames(comm_ophi)
}

## ses.mpds mean per grid cell
ses.mpds.ophi.X<-lapply(ses.mpds.ophi, as.matrix)
ses.mpds.ophi.X<-do.call(cbind, ses.mpds.ophi.X)
mean_ses_pd_ophi<-as.data.frame(rowMeans(ses.mpds.ophi.X, na.rm=T))
rownames(mean_ses_pd_ophi)<-rownames(comm_ophi)
mean_ses_pd_ophi<-as.data.frame(mean_ses_pd_ophi)
mean_ses_pd_ophi$grids<-rownames(mean_ses_pd_ophi)
names(mean_ses_pd_ophi)<-c("mean.ses.pd","grids")
write.csv(mean_ses_pd_ophi, file="mean_ses_pd_ophi.csv")

### PD ###
set.seed(45)
pd.ophi.loess2<-list()
for(i in 1:10){
  print (i)
  imp.tree <- congeneric.impute(ophiure_tree, ophiura_to_insert)
  c.data<-match_phylo_comm(imp.tree, sparse_grid3_ophi, delete_empty_rows = F)
  X<-pd(as.matrix(c.data$comm), c.data$phy, include.root = T)
  X$grids<-rownames(c.data$comm)
  Y<-subset(X, !PD==0)
  resids_loess <- resid(gam(Y$PD ~ s(Y$SR)))
  pd.ophi.loess2[[i]]<-cbind(resids_loess, Y)
}

library(data.table)
library(dplyr)
pd.ophi.loess2<-lapply(pd.ophi.loess2,as.data.table)
pd.ophi.loess2<-rbindlist(pd.ophi.loess2, use.names=T, fill=T)
mean_loess_PD_ophi2<- pd.ophi.loess2%>%group_by(grids)%>%summarise_all("mean")
mean_loess_PD_ophi2$scaled_loess_resid<-scale(mean_loess_PD_ophi2$resids_loess)
write.csv(mean_loess_PD_ophi, file="mean_loess_PD_ophi.csv")

ggplot(mean_loess_PD_ophi2, aes(PD,resids_loess)) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

st<-file.choose()
mean_loess_PD_ophi<-read.table(t, dec=".", sep=",", header=T)
##ATTENTION  IL FAUT AUSSI LA RICHESSE EN MOYENNE PAR PLOT IMPUTE!!!
mean_loess_PD_ophi<-merge(grid_aq_3, mean_loess_PD_ophi, by="grids", all.x=T, all.y=F)
###################
tiff("loess_PD_ophi.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(mean_loess_PD_ophi,!is.na(PD)), aes(fill=scaled_loess_resid))+
  #scale_fill_viridis(option="A", direction = 1)+
  scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=Coastline_sf, fill="white")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", size=0.7)+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="black", size=0.7,linetype = "dashed")+
  theme_minimal()
dev.off()  
###################
tiff("loess_SR_ophi.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(mean_loess_PD_ophi,!is.na(PD)), aes(fill=SR))+
  #scale_fill_viridis(option="A", direction = 1)+
  scale_fill_viridis(option="C", direction = 1, name="Species\nRichness")+
  geom_sf(data=Coastline_sf, fill="white")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", size=0.7)+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="black", size=0.7,linetype = "dashed")+
  theme_minimal()
dev.off()  
### .ses.mpd mets des NA quand une cellule comporte qu'une seule espece !!

# compute phylogenetic beta diversity btw cells
# make a mean with random add of missing species
set.seed(42)
phylo.beta.ophi<-list()
for(i in 1:1000){
  print (i)
  imp.tree <- congeneric.impute(ophiure_tree, ophiura_to_insert)
  imp.tree<-as.phylo(imp.tree)
  #c.data <- comparative.comm(imp.tree, comm)
  c.data<-match_phylo_comm(imp.tree, sparse_grid3_ophi, delete_empty_rows = F)
  phylo.beta.ophi[[i]] <- phylobeta(c.data$comm, c.data$phy)$phylo.beta.sor
}

class(phylo.beta.ophi[[1]])
str(phylo.beta.ophi)

#there are NA NaN in mean distance matrix
Y <- do.call(cbind, phylo.beta.ophi)
Y <- array(Y, dim=c(dim(phylo.beta.ophi[[1]]), length(phylo.beta.ophi)))

mean_phylo_beta_ophi<-colMeans(aperm(Y, c(3, 1, 2)), na.rm = TRUE)
names(mean_phylo_beta_ophi)<-rownames(comm_ophi)
mean_phylo_beta_ophi<-as.dist(mean_phylo_beta_ophi)
dim(mean_phylo_beta_ophi)

## phyloregions ----
dim(com.data.ophi$comm)
phylo_beta_ophi <- phylobeta(com.data.ophi$comm, com.data.ophi$phy)
#select the less distorting clustering method
select_linkage(phylo_beta_ophi[[1]])
#select optimal number of clusters with selected method
optim<-optimal_phyloregion(phylo_beta_ophi[[1]], k=20, method="centroid")
plot(optim$df$k, optim$df$ev)
grid_aq3_sf<-st_as_sf(grid_aq_3_pr)
#names(grid_aq3_sf)<-c("cell_name3", "geometry")
#ophiura_ses_pd<-merge(grid_aq3_sf, mean_ses_pd_ophi, by="grids", all.x =T, all.y=F)
X11()
##phyloregionalisation
y <- phyloregion(phylo_beta_ophi[[1]], shp=grid_aq_3_sp, k=8, method="average")
summary(y)
y$membership
y$region.df
ophi_phyloreg_sf<-y$shp
plot(ophi_phyloreg_sf)
ophi_phyloreg_sf<-st_as_sf(ophi_phyloreg_sf)
plot(ophi_phyloreg_sf["ED"])

PD_ophi<-merge(PD_ophi, y$region.df, by="grids", all.x=T, all.y=F)
head(PD_ophi)
library(RColorBrewer)

phylo_nmds<-y$NMDS
phylo_nmds<-y$NMDS$points
class(phylo_nmds)
phylo_nmds<-as.data.frame(phylo_nmds)
phylo_nmds$phyloreg<-rownames(phylo_nmds)

tiff("nmds_phyloregions_ophi.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ophi_nmds<-ggplot()+geom_point(data=phylo_nmds, aes(x=MDS1, y=MDS2, fill=phyloreg), pch=21, size=6)+
  scale_fill_brewer(palette = "Dark2") +   # customized color palette
  labs(fill = "")+
  geom_text(data=phylo_nmds, aes(x=MDS1, y=MDS2, label=phyloreg), vjust=0.3, size=4)+
  theme_minimal()+
  theme(legend.position="none")
dev.off()

plot_NMDS(y, cex=4)
text_NMDS(y)


###########
## PLOTS ##
###########
library(ggplot2)
library(ggeffects)
library(ggpubr) #for multiple lots together
library(CCAMLRGIS)
library(ggplot2)
library(viridis)
tiff("ophi_SES_PD_phyloreg2_1000_3.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(ophiura_ses_pd,!is.na(mean.ses.pd)), aes(fill=mean.ses.pd))+
  scale_fill_viridis(option="A", direction = 1)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, fill=NA, color="black", size=0.7)
dev.off()

#######################
ggplot()+
  geom_sf(data=subset(ophiura_ses_pd,pd_obs_p<0.05), fill=NA, color="red", size=0.8)

tiff("ophi_richness2.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(ophiura_ses_pd,!is.na(richness)), aes(fill=richness))+
  #scale_fill_viridis(option="A", direction = 1)+
  scale_fill_viridis(option="D", direction = 1)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", size=0.7)+
  theme_minimal()
dev.off()


