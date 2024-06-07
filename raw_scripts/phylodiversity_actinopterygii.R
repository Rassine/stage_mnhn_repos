library(rgeos) 
library(rgdal)
library(sp)
library(sf)
#library(quantarcticR)
library(ggplot2)
library(dplyr)
library(viridis)
library(picante)  
library(phyloregion)
library(tidyr) 
library(CoordinateCleaner)
library(rfishbase)
library(pez)
library(phyloregion)
library(picante)
library(vegan)
#library(SOmap)

#####  Actinopterygii data  ###### ----
#file after deleting points on land
actin<-file.choose()
actinopter<-read_sf("inside_actinopterygii.shp",
                    layer="inside_actinopterygii")
actinopter <- cbind(actinopter,st_coordinates(actinopter))#x and Y only
actinopter_df<-st_drop_geometry(actinopter)
# inside actinopterygii comes directly from 'on land' script
setdiff(actinopter$species, inside_actinopterygii$species)
setdiff(inside_actinopterygii$species,actinopter_worms$species)
actinopter$species<-gsub(" ", "_", actinopter$species)
#setdiff(inside_actinopterygii$species, animaliaBAQ_filt$species)
length(unique(inside_actinopterygii$species))
length(unique(actinopter$species))

## too many differnces 
## but I tol up to 40 degrees not 45 like last time!! so less species in actinos_worms

actinopter_worms<-file.choose()
actinopter_worms<-read.csv(actinopter_worms, dec=".", sep=";", header=T)
actinopter_worms<-actinopter_worms[,c(1,6)]
actinopter_worms<-plyr::rename(actinopter_worms, c("scientificname"="species", "scientificname.1"="newscientificname"))
actinopter_worms$newscientificname<-gsub(" ", "_", actinopter_worms$newscientificname)
actinopter_worms$species<-gsub(" ", "_", actinopter_worms$species)
#test to merge species new names (after synonymy)
#merge new names with data without points on land
actinopter<-merge(actinopter,actinopter_worms,by="species",all.x=TRUE, all.y=FALSE)
table(actinopter$source)
length(unique(actinopter$newscientificname))
actinopter<-actinopter%>%
  filter(!is.na(newscientificname))
actinopter <- actinopter %>% mutate_all(na_if,"")

#inside_actinos_baq<-merge(inside_actinos_baq,
  #                        actinopter_worms,by="species",all.x=TRUE, all.y=FALSE)
#inside_actinos_baq<-inside_actinos_baq%>%
 # filter(!is.na(species))%>%
  #select(!c(newscientificname.x, newscientificname.y ))

setdiff(inside_actinos_baq$newscientificname, phylo_fish_maker$tip.label )
setdiff(inside_actinos_baq$newscientificname, phylo_fish_maker2$tip.label )# le mieux

#phylo_fishes_baq<-congeneric.merge(phylo_fishes_baq, spp_inset, split="_")

### CHECK FOR DOUPLICATA AGAIN in the occurrence data frame in new species names!!!!####
fish_names_phy<-unique(actinopter_worms$newscientificname)

library(FishPhyloMaker)
fish_names_phy2<-FishTaxaMaker(fish_names_phy, allow.manual.insert = FALSE)
fish_s_f_o<-fish_names_phy2[[2]]#duplicates deleted according to fishbase but there are 15 spp not found in the fishbase
fish_all<-fish_names_phy2[[1]]#with duplicates 
write.csv2(fish_all, file="fish_all_fishbase.csv")
write.csv2(fish_s_f_o, file="fish_s_f_o.csv")
fish_sfo_anna<-file.choose()
fish_sfo_anna<-read.csv(fish_sfo_anna, dec=".", sep=";", header=T)
## delete duplicates!
fish_sfo_anna$s[duplicated(fish_sfo_anna$s)]
fish_sfo_anna<-fish_sfo_anna%>%
  distinct(s,f, o,  .keep_all = TRUE)

fish_not_found<-as.data.frame(fish_names_phy2[[3]])#not found in fishbase cause not upgraded lately
dup<-fish_all$valid_names[duplicated(fish_all$valid_names)]
dup#4 + 12NAs
#these names are in disaccordance between worms & fishbase 
#(eschmeyer database https://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp)

library(worms)
fish_names_phy<-gsub("_", " ", fish_names_phy)
#prepare species, fam, orred froom fishbase
#fish_s_f_o2<-unique(fish_all[,c("valid_names","Family","Order")])
#fish_s_f_o3<-subset(fish_all, !is.na(Family))
#g<-setdiff(fish_s_f_o3$valid_names, fish_s_f_o$s)
#g2<-setdiff(fish_s_f_o3$user_spp, fish_s_f_o$s)
#h<-subset(fish_all, valid_names=="Lagiacrusichthys macropinna")
#prepare species family and order by yourself
actinopter_phy<-actinopter%>%select(newscientificname, family, order)
actinopter_phy<-st_drop_geometry(actinopter_phy)
actinopter_phy$newscientificname<-gsub( " ", "_", actinopter_phy$newscientificname)
actinopter_phy<-unique(actinopter_phy[,c("newscientificname", "family", "order")])

# phylo from fish tree of life https://fishtreeoflife.org/
fish_tree<-file.choose()
fish_tree<-read.tree(fish_tree)

which_spp_add<-whichFishAdd(fish_sfo_anna)

#make a fish synthetic phylogeny for all my species ONLY with FishPhyloMaker
## will add species minus 15 not found in FishBase
res_phylo <- FishPhyloMaker(data = fish_s_f_o,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE)

res_phylo2 <- FishPhyloMaker(data = fish_sfo_anna,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE)


insertions_data<-res_phylo[[2]]
insertions_data2<-res_phylo2[[2]]
table(insertions_data2$insertions)
plot.phylo(res_phylo2[[1]], show.tip.label=T, cex = 0.2)
is.ultrametric(res_phylo2[[1]])
hist(res_phylo2[[1]]$edge.length, breaks = 200)
#write a newly made tree
#phylo_fish_maker<-res_phylo$Phylogeny
phylo_fish_maker2_nov2022<-res_phylo2$Phylogeny# phylo to use ----
write.tree(phylo_fish_maker2_nov2022, file="phylo_fish_maker_nov2022")
phylo_fish_maker2_nov2022$tip.label[duplicated(phylo_fish_maker2_nov2022$tip.label)]
# "Lampanyctus_ater"    "Lampanyctus_achirus"

fish_phylo_names<-phylo_fish_maker2$tip.label
write.csv2(fish_phylo_names, file="fish_phylo_names.csv")
#write.tree(new_phylo_fish_maker, file="final_fish_phylo.tre")

## for species not on a phylogeny ----
library(pez)
fish_species_insert<-setdiff(actinopter$newscientificname, phylo_fish_maker2$tip.label )
phylo_fish_maker2<-congeneric.merge(phylo_fish_maker2, fish_species_insert, split="_")
setdiff(fish_species_insert,phylo_fish_maker2$tip.label)
# some species changed classif in jan 2023:
actinopter$newscientificname<-gsub("Nannobrachium_atrum", "Lampanyctus_ater",actinopter$newscientificname)
actinopter$newscientificname<-gsub("Nannobrachium_phyllisae", "Lampanyctus_phyllisae",actinopter$newscientificname)
actinopter$newscientificname<-gsub("Nannobrachium_achirus", "Lampanyctus_achirus",actinopter$newscientificname)
actinopter$newscientificname<-gsub("Mesobius_antipodum", "Mesovagus_antipodum",actinopter$newscientificname)

actinos<-subset(actinopter, newscientificname %in% phylo_fish_maker2$tip.label)
setdiff(actinopter$newscientificname,phylo_fish_maker2$tip.label)

setdiff(fish_s_f_o$s, actinos$newscientificname)
setdiff(actinos$newscientificname,fish_s_f_o$s)
setdiff(fish_s_f_o$s,phylo_fish_maker2$tip.label)

# Leptoscopus_macropygus - ancienne Uranoscopus macropygus

###########################################################################
## occurrence data actinos - on phylogeny only ----
phylo_fish_maker2$tip.label<-gsub(" ", "_", phylo_fish_maker2$tip.label)
#actinos2022<-subset(inside_actinos_baq, newscientificname %in% phylo_fishes_baq$tip.label)

#1173 species and not 1203
st_crs(actinos)#3031
plot(st_geometry(actinos),  key.pos = NULL, axes = TRUE)

length(unique(actinos$newscientificname))
#check for duplicates
#actinos_nodup<-actinos %>%
# mutate(as.data.frame(st_coordinates(actinos)))%>%
#distinct(X,Y,newscientificname, .keep_all = TRUE) #Works

#actinos <- cbind(actinos,st_coordinates(actinos))#x and Y only
actinos_df<-st_drop_geometry(actinos)
actinos_df<-actinos_df%>%
  distinct(X,Y,newscientificname,.keep_all = TRUE)

actinos3031 = st_as_sf(actinos_df, coords = c("X", "Y"), crs = 3031)
setdiff(actinos_df$newscientificname, phylo_fish_maker2$tip.label)

# Molleweide projection for equal areas surfaces preservation
prj_moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"
prj_laea <- "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 

actinos_moll = st_transform(actinos3031, crs = prj_moll)
actinos_laea<-st_transform(actinos3031, crs=prj_laea)

length(unique(actinos_laea$newscientificname))
plot(st_geometry(actinos3031))
plot(st_geometry(actinos_laea))
plot(st_geometry(grid_aq_3_PF), add=T)

grid_aq_3<-st_transform(grid_aq_3, crs=prj_laea)
grid_aq_3_PF<-st_transform(grid_aq_3_PF, crs=prj_laea)
# 

actinos_gridPF<-st_intersection(actinos_laea, grid_aq_3_PF)#
length(unique(actinos_gridPF$newscientificname))
length(unique(actinos_gridPF$grids))
## attention duplication de certains points d'occ qui tombent entre deux grid cells !
# les enlever:
actinos_gridPF <- cbind(actinos_gridPF,st_coordinates(actinos_gridPF))#x and Y only
actinos_gridPF<-actinos_gridPF%>%
  distinct(X,Y,newscientificname,.keep_all = TRUE)
length(unique(actinos_gridPF$newscientificname))
#actinos_gridPF$year<-as.numeric(as.character(actinos_gridPF$year))
### check for NAs ----
summary(actinos_gridPF)
100 * sum(is.na(actinos_gridPF$year), na.rm = TRUE) / length(actinos_gridPF$year)
actinos_gridPF %>% summarise_all(~ sum(is.na(.)))
table(actinos_grid$newscientificname, useNA="always")
# do no use !!!!!!!
#actinos_grid2<-st_join(actinos3031, grid_aq_3, join=st_nearest_feature, left = T)
# DO NOT USE !!

## count species total abundance ----
actinos_abund<-actinos_gridPF%>%
  dplyr::count(newscientificname)
actinos_abund<-st_drop_geometry(actinos_abund)
actinos_gridPF<-merge(actinos_gridPF, actinos_abund, by="newscientificname", all.x=T, all.y=F)
rm(actinos_abund)

## delete species with 1 abundance -species at the limit of their range areas (could be)
table(actinos_gridPF$n)# 223 species are occ = 1
occ1_actinos<-actinos_abund%>%filter(n == 1)
write.csv2(occ1_actinos, file="occ1_actinos.csv")

table(actinos_abund$n)

actinos_gridPF_1<-actinos_gridPF%>%
  filter(!n <3)%>%
  filter(!newscientificname %in% c("Hisonotus_hungy", "Sprattus_fuegensis", 
                                   "Salmo_trutta","Parabrotula_plagiophthalma","Gonorynchus_gonorynchus",
                                  "Halichoeres_dispilus", "Notolychnus_valdiviae" ))
length(unique(actinos_gridPF_1$grids))
length(unique(actinos_gridPF_1$newscientificname))#
table(actinos_gridPF_1$n)

## most high ED species were errors, not AQ species ! Even after taking only species limited by PF
## Galaxias_maculatus - fresh water mostly n=1
## Sprattus_fuegensis, tierre de fuego n= 3
## Gonorynchus_gonorynchus -   n=1
## (Woodsia_meyerwaardeni)   n=7 South Atlantic: probably reaches southern Africa, nearest is 37°08'S, 5°23'E. Southwest Pacific: New Zealand
## Macroparalepis_macrogeneion n=2 , MNHN !! South Atlantic: discontinuous from continental slope off South America to South Africa. Indo-West Pacific: off Australia (Ref. 7300) and New Zealand (Ref. 5755). Southeast Pacific: Chile (Ref. 9068).
## Lestidiops_similis n=1 atlantic fish, until 40, 35 degrees
## Macroparalepis_affinis n=2 Atlantic Ocean: antitropical in distribution. Known only from stomachs of Alepisaurus in the South Pacific.
## Vinciguerria_lucetia n=3 Western Central Pacific: Papua New Guinea. Eastern Pacific: throughout the California Current region, usually south of Point Conception and seaward of shelf (Ref. 35839); also in Chile (Ref. 9068).
## (Vinciguerria_attenuata) n=6
## Odontognathus_mucronatus tropical and freshwater n=1
## (Phosichthys_argenteus) n=10
## (Ichthyococcus_australis) n=1 Circumglobal in the sub-tropical convergence region of the southern hemisphere with records between 30°S and 40°S in the southeast Atlantic.
## Caristius_macropus (peut etre garder) n=1 North Pacific: from subtropical waters to the Bering Sea and the coast of Alaska to where it is apparently drifted by warm currents.
## Parabrotula_plagiophthalma n=1
## Metelectrona_ventralis n=77
## Cyclothone_obscura
## Salmo_trutta n =1, european spp !!!
## Rhamdia_quelen = freshwater!
## 
table(actinos_gridPF_1$family)
length(unique(actinos_gridPF_1$family))# 51, and in Atlas there are 47

actinos_families<- actinos_gridPF_1 %>% # Applying group_by & summarise
  group_by(family,genus, .add=T) %>%
  summarise(count = dplyr::n_distinct(newscientificname))

actinos_familiesSR<- actinos_gridPF_1 %>% 
  st_drop_geometry()%>%
  group_by(grids, .add=T) %>%
  summarise(sr_fam = dplyr::n_distinct(family))

PD_actinos<-merge(PD_actinos, actinos_familiesSR, by="grids", all.x=T, all.y=F)
names(PD_actinos)
# the most frequent actino families
Nototheniidae<-subset(actinos_gridPF_1, family=="Nototheniidae")
length(unique(Nototheniidae$newscientificname))#42 species
plot(Nototheniidae["newscientificname"])
plot(Nototheniidae["genus"])

Myctophidae<-subset(actinos_gridPF_1, family=="Myctophidae")
length(unique(Myctophidae$newscientificname))#40 species 
unique(Myctophidae$newscientificname)
plot(Myctophidae)

Artedidraconidae<-subset(actinos_gridPF_1, family=="Artedidraconidae")
length(unique(Artedidraconidae$newscientificname))#23 species 
unique(Artedidraconidae$newscientificname)
plot(st_geometry(Artedidraconidae))

Liparidae<-subset(actinos_gridPF_1, family=="Liparidae")
length(unique(Liparidae$newscientificname))#33
plot(st_geometry(Liparidae))

Zoarcidae<-subset(actinos_gridPF_1, family=="Zoarcidae")
length(unique(Zoarcidae$newscientificname))#24
plot(st_geometry(Zoarcidae))

plot(actinos_gridPF_1["family"])# spatial pattern in families 
ggplot()+
geom_sf(data=actinos_gridPF_1, aes(colour=family))
## maybe delete preserved specimens ?
table(actinos_gridPF_1$bssOfRc)

# check species by its range area in fish base - some are not AQ at all
actinos_list<-unique(actinos_gridPF_1$newscientificname)
actinos_list_all<-unique(actinos_grid$newscientificname)

library(rfishbase)
actinos_distrib<-rfishbase::diet(species_list=actinos_list)
# rfishbase does not work
write.csv2(actinos_list, file="actinos_list.csv")

actinos_list_checkup<-read.csv("actinos_list.csv", dec=".", sep=";", header=T, na.strings = "")
actinos_list_checkup$newscientificname<-gsub(" ", "_", actinos_list_checkup$newscientificname)
actinos_gridPF_1<-merge(actinos_gridPF_1, actinos_list_checkup, by="newscientificname", all.x=T, all.y=F)

# keep:
#	Merluccius hubbsi -  20 occurrences
#	Harpagifer bispinis - 24 occurrences, observed in subantarctic islands - DO NOT DELETE*
fd<-actinos_gridPF_1%>%
  filter(newscientificname %in% c("Merluccius_hubbsi"))

ggplot()+geom_sf(data=fd,shape = 16, size=2, alpha = 0.4, color = "purple4" )+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red2", size=6,linetype = "dashed")

spp_to_delete<-subset(actinos_gridPF_1, !is.na(delete))
length(unique(spp_to_delete$newscientificname))
spp_to_delete<-spp_to_delete%>%
  filter(!newscientificname %in% c("Merluccius_hubbsi","Harpagifer_bispinis"))
spp_to_delete<-unique(spp_to_delete$newscientificname)
# delete species cheked in fish base as not antarctic (subtropical, temperate, NOT occsional)
actinos_gridPF_1<-actinos_gridPF_1%>%
  filter(!newscientificname %in% spp_to_delete)
length(unique(actinos_gridPF_1$newscientificname))
plot(st_geometry(actinos_gridPF_1))

table(actinos_gridPF_1$occasional)
actinos_gridPF_1$occasional<-gsub("no ", "no",  actinos_gridPF_1$occasional)
plot(actinos_gridPF_1["occasional"])# spatial pattern in families 

## species richness through years - also check for sampling effort
SR_year_actinos<- actinos_gridPF_1 %>% # Applying group_by & summarise
  group_by(year) %>%
  dplyr::summarise(count = dplyr::n_distinct(newscientificname))

#tiff("actinos336_SR_years.tiff",width=6,height=4,res=400, units="in", compression ="lzw")
#SR_year_actinos%>%
#  ggplot(aes(x=year, y=count))+geom_jitter()+geom_line()
#dev.off()

100 * sum(is.na(actinos_gridPF_1$year), na.rm = TRUE) / length(actinos_gridPF_1$year)

#st_write(actinos_grid, "actinos_grid.shp", layer="actinos_grid", driver="ESRI Shapefile")
## create presence - abscence data, not abundance as before THEN delete grid cells with one species (and not 1 abundance)
actinos_grid_pres_abs<-actinos_gridPF%>%
  dplyr::select (grids,newscientificname)%>%
  filter(!is.na(grids))%>%
  distinct(grids,newscientificname, .keep_all=T)#delete duplicated sp in unique cell==pres/abs
length(unique(actinos_grid_pres_abs$grids))
length(unique(actinos_grid_pres_abs$newscientificname))#

# delete cells with SR = 1
SR_actinos<-actinos_gridPF_1%>%
  group_by(grids, .add=T) %>%
  st_drop_geometry()%>%
  summarise(count = dplyr::n_distinct(newscientificname))
actinos_gridPF_1<-merge(actinos_gridPF_1, SR_actinos, by = "grids", all.x=T, all.y=F)
actinos_gridPF_1<-actinos_gridPF_1%>%
  filter(!count == 1)
length(unique(actinos_gridPF_1$newscientificname))
length(unique(actinos_gridPF_1$grids))
table(actinos_gridPF_1$count)

write.csv2(actinos_gridPF_1, file="actinos_gridPF_1_10823.csv")
actinos_gridPF_1<-read.csv2("actinos_gridPF_1_10823.csv")
length(unique(actinos_gridPF_1$newscientificname))
length(unique(actinos_gridPF_1$grids))

#create a community matrix ----
sparse_grid3_actinos <- long2sparse(actinos_gridPF_1, grid="grids", species="newscientificname")
dim(sparse_grid3_actinos)
rownames(sparse_grid3_actinos)
## ses_PD with phyloregion package
library(phyloregion)
set.seed(34)# for SES calculations and zscore
#delete species that  re not on the phylogeny !
com.data.actinos<-match_phylo_comm(phylo_fish_maker2, sparse_grid3_actinos, delete_empty_rows = F)
dim(com.data.actinos$comm)
comm_actinos<-com.data.actinos$comm
comm_actinos<-as.data.frame(as.matrix(comm_actinos))
setdiff(colnames(sparse_grid3_actinos),colnames(comm_actinos))
com.data.actinos$phy

ses_PD_actinos<-PD_ses(com.data.actinos$comm, com.data.actinos$phy, model="tipshuffle", reps=999)
#ses_PD_actinos2<-ses.mpd.anna(comm_actinos, cophenetic(phylo_fish_maker2), null.model="taxa.labels")
PD_actinos<-merge(grid_aq_3_PF, ses_PD_actinos, by="grids", all.x =T, all.y=F)
PD_actinos<-PD_actinos%>%
  filter(!is.na(PD_obs))#no NA allowed for following models
## richness is not richness but total abundance in each grid cell !!! Use the one from pez results SR
#library(pez)
#is.rooted(phylo_fishes_baq)
#c.data.actino<-comparative.comm(phylo_fishes_baq, as.matrix(comm_actinos))
#dim(c.data.actino$comm)
PE<-phylo_endemism(com.data.actinos$comm, com.data.actinos$phy)
WE<-weighted_endemism(com.data.actinos$comm)
PD_actinos <- merge(PD_actinos, data.frame(grids=names(PE), PE=PE), by="grids")
PD_actinos <- merge(PD_actinos, data.frame(grids=names(WE), WE=WE), by="grids")

PD_actinos$std_PE <- as.vector(scale(resid(loess(PE ~ WE, data=PD_actinos))))

#is.rooted(com.data.actinos$phy)
#class(com.data.actinos$phy)
#actino_phy=multi2di(com.data.actinos$phy)

pd_actinos_pez<-as.data.frame(picante::pd(comm_actinos, com.data.actinos$phy, include.root = F))
pd_actinos_pez$grids<-rownames(pd_actinos_pez)
PD_actinos<-merge(PD_actinos, pd_actinos_pez, by="grids", all.x =T, all.y=F)
# delete celles with one species !
PD_actinos<-PD_actinos%>%
  filter(!is.na(PD))
#also delete species that where there only
comm_actinos<-comm_actinos[(rownames(comm_actinos) %in% PD_actinos$grids),]
comm_actinos<-comm_actinos[,!(colSums(comm_actinos) == 0)]

dim(comm_actinos)

### quadratic entropy ---- 
library(adiv)
library(adephylo)
phylo_dist_actinos<-distTips(phylo_fishes_baq)
class(phylo_dist_actinos)
rao_QE<-QE(comm_actinos, phylo_dist_actinos, formula="QE")
rao_QE$grids<-rownames(rao_QE)
names(rao_QE)<-c("rao_QE","grids")
PD_actinos<-merge(PD_actinos, rao_QE, by="grids", all.x =T, all.y=F)
# with pres abs
rao_QE<-QE(comm_actinos, phylo_fishes_baq, formula="QE")
rao_QE$grids<-rownames(rao_QE)
names(rao_QE)<-c("rao_QE_PA","grids")
PD_actinos<-merge(PD_actinos, rao_QE, by="grids", all.x =T, all.y=F)
# inverse of QE
r_entropy<-Rentropy(comm_actinos, phylo_fishes_baq, scale=T)
r_entropy$grids<-rownames(r_entropy)
names(r_entropy)<-c("r_entropy","grids")
PD_actinos<-merge(PD_actinos, r_entropy, by="grids", all.x =T, all.y=F)

# stndardized PD ----
PD_actinos$std_resid_loess <- as.vector(scale(resid(loess(PD_obs ~ SR, data=PD_actinos))))
library(mgcv)
PD_actinos$std_resid_gam <- as.vector(scale(resid(gam(PD_obs ~ s(SR), data=PD_actinos))))
PD_actinos$std_zscore<-as.vector(scale(PD_actinos$zscore))

#PD_actinos$std_rao_QE<-as.vector(scale(resid(loess(rao_QE ~ SR, data=PD_actinos))))
#PD_actinos$std_rao_QE_PA<-as.vector(scale(resid(loess(rao_QE_PA ~ SR, data=PD_actinos))))
#PD_actinos$rentropy_std<-as.vector(resid(loess(r_entropy ~ SR, data=PD_actinos)))

## check for spacial autocorrelation ! ----
loess_model_actinos_pd<-loess(PD_obs ~ SR, data=PD_actinos)
plot(loess_model_actinos_pd)
centroids_mpa_pd<- cbind(centroids_mpa_pd,st_coordinates(centroids_mpa_pd))#x and Y only
plot(st_geometry(centroids_mpa_pd))

# correlation spatiale check
library(ncf)
ncf_pd_loess_actinos <- correlog(centroids_mpa_pd$X, centroids_mpa_pd$Y, 
                                 resid(loess_model_actinos_pd),
                                  increment=250, resamp=99, latlon=T)#OK !
lisa_pd_loess_actinos <- lisa(centroids_mpa_pd$X, centroids_mpa_pd$Y, 
                                 resid(loess_model_actinos_pd),
                                 neigh = 500,resamp=99, latlon=T)#OK !
plot(ncf_pd_loess_actinos)
plot(lisa_pd_loess_actinos)

bubble(spdata, "var", main = "SES PD", xlab = "X-coordinates", 
       ylab = "Y-coordinates", maxsize = 2, col = c("#d01c8b", "#4dac26"), fill = FALSE)

spline_correl<-spline.correlog(x = centroids_mpa_pd$X, y = centroids_mpa_pd$Y,
                               resid(loess_model_actinos_pd),resamp = 999)
plot.correlog(spline_correl)

moranI<-ncf_pd_loess_actinos$correlation[1:15]
Distance<-ncf_pd_loess_actinos$mean.of.class[1:15]
plot(Distance, moranI,ylim = c(-0.2, 0.2), xlab = "Distance (m)", ylab = "Moran's I")
points(Distance, moranI, col=c("red","blue","blue","blue", "blue","blue","blue","red","blue","blue","blue" ), pch=18, cex=1.4)
abline(h=0)
#abline(v=300)
ncf_pd_loess_actinos$p[1:15]
# fin de ncf ----
#################################
library(ggplot2)
PD_actinos%>%ggplot(aes(x=SR, y=richness))+geom_point()+geom_jitter()+
coord_cartesian(xlim=c(0, 110))#tSR increases with sampling effort = number of species occurrences
st_write(PD_actinos, "PD_actinos.shp", layer="PD_actinos", driver="ESRI Shapefile")
#mpd_actinos_pez<-.ses.mpd(c.data.actino, null.model="taxa.labels", permute=999)

## MPD and MNTD ----
# use of comm_pycnos lready without 1 species cells
phy.dist_actinos<-cophenetic(com.data.actinos$phy)
actino_mpd<-as.data.frame(picante::mpd(comm_actinos, phy.dist_actinos))
names(actino_mpd)<-"mpd"
actino_mpd$grids<-PD_actinos$grids
head(actino_mpd)

ses_mpd_actino<-ses.mpd(comm_actinos, phy.dist_actinos, null.model="taxa.labels", runs=999)
head(ses_mpd_actino)
ses_mpd_actino$grids<-row.names(ses_mpd_actino)
ses_actinos<-merge(grid_aq_3_PF, ses_mpd_actino, by="grids", all.x =T, all.y=F)

actino_mntd<-as.data.frame(picante::mntd(comm_actinos, phy.dist_actinos))
names(actino_mntd)<-"mntd"
actino_mntd$grids<-PD_actinos$grids
head(actino_mntd)

PD_actinos<-merge(PD_actinos, actino_mpd, by="grids", all.x=T, all.y=F)
PD_actinos<-merge(PD_actinos, actino_mntd, by="grids", all.x=T, all.y=F)

ses_mntd_actino<-ses.mntd(comm_actinos, phy.dist_actinos, null.model="taxa.labels", runs=999)
ses_mntd_actino$grids<-row.names(ses_mntd_actino)
ses_actinos<-merge(ses_actinos, ses_mntd_actino, by="grids", all.x =T, all.y=F)
ses_actinos<-ses_actinos%>%
  filter(!is.na(ses_actinos$ntaxa.x))

PD_actinos$loess_mpd <- as.vector(scale(resid(loess(mpd ~ SR, data=PD_actinos))))
PD_actinos$loess_mntd <- as.vector(scale(resid(loess(mntd ~ SR, data=PD_actinos))))

write.csv(PD_actinos, file="PD_actinos_measures.csv")

PD_actinos<-merge(PD_actinos, st_drop_geometry(ses_actinos), by="grids", all.x=T, all.y=F)

PD_actinos_noXY<-PD_actinos%>%
  st_drop_geometry(PD_actinos)%>%
  select(grids, PD, std_resid_loess, SR, mpd, mntd, 
         mpd.obs.z,mntd.obs.z,loess_mpd, loess_mntd,
         PE, WE, std_PE,)

str(PD_actinos_noXY)

pairs(PD_actinos_noXY[,c(2:13)], lower.panel = panel.smooth,  upper.panel = panel.cor, 
      cex = 1.5, pch= 1, cex.labels = 1, gap=1/4)

## will use std_resid_loess

#####################
##### species ED ----
####################
## for data without & occ species
## top 5 % ED ----
ED_actinos_PF_1<-picante::evol.distinct(com.data.actinos$phy)
names(ED_actinos_PF_1)<-c("newscientificname", "ED_tot")
actinos_gridPF_1<-merge(actinos_gridPF_1, ED_actinos_PF_1, by="newscientificname", all.x=T, all.y=F)
cor.test(actinos_gridPF_1$n, actinos_gridPF_1$ED_tot)## 
ggplot(actinos_gridPF_1, aes(ED_tot, n))+geom_point()+geom_smooth(method='lm')

data <- actinos_gridPF_1 %>% arrange(desc(actinos_gridPF_1$ED_tot))%>%select(newscientificname, ED_tot)
obs <- nrow(data) 
topED5_PF<-data %>% filter(row_number() < obs * 0.025)
names(topED5_PF)<-c("newscientificname", "ED_top2.5")
actinos_gridPF_1<-merge(actinos_gridPF_1, topED5_PF, by="newscientificname", all.x=T, all.y=F)
cor.test(actinos_gridPF_1$n, actinos_gridPF_1$ED_top5, na.omit=T)
length(unique(topED5_PF$newscientificname))

# with phyloregion
#PD_actinos_inexPD$hot_PDr_loess_phyloregion<-phyloregion::hotspots(PD_actinos_inexPD$loess_PDrare_15, prob = 10)

# PD top10% captured in hotspots
sum(PD_actinos_inexPD$PDrare_15[PD_actinos_inexPD$hot_PDr_loess_phyloregion==1], na.rm=T)
sum(PD_actinos_inexPD$PDrare_15, na.rm=T)

sum(com.data.actinos$phy$edge.length)
sum(PD_actinos$pd_added, na.rm=T)
sum(PD_actinos$PD_obs[!is.na(PD_actinos$hotspot10_resid_loess)], na.rm=T)

PD_actinos_inexPD$PDr_percentage = 
  PD_actinos_inexPD$PDrare_15/
  sum(com.data.actinos$phy$edge.length)
sum(PD_actinos_inexPD$PDr_percentage[PD_actinos_inexPD$hot_PDr_loess_phyloregion==1], na.rm=T)
#####

data <- centroids_mpa_pd %>% 
  select(SR, grids)%>%
  arrange(SR)
obs <- nrow(data) 
coldspot10_SR<-data %>% filter(row_number() < obs * 0.1)
names(coldspot10_SR)<-c("coldspot10_SR", "grids", "geometry")
plot(coldspot10_SR)

centroids_mpa_pd<-merge(centroids_mpa_pd, st_drop_geometry(hotspot10_SR), 
                        by="grids", all.x =T, all.y=F )
centroids_mpa_pd<-merge(centroids_mpa_pd, st_drop_geometry(coldspot10_SR), 
                        by="grids", all.x =T, all.y=F )

hotspots_actino<-centroids_mpa_pd%>%select(grids, hotspot10_SR,coldspot10_SR)
hotspots_actino<-centroids_mpa_pd%>%select(grids, hotspot10_SR,coldspot10_SR, hotspot10_resid_loess, coldspot10_resid_loess)
PD_actinos<-merge(PD_actinos, st_drop_geometry(hotspots_actino), by="grids", all.x=T, all.y=F)

library(dplyr)
PD_actinos<-PD_actinos %>% 
  group_by(std_resid_loess)%>%
  mutate(percent_std_loess = proportions(std_resid_loess))

## another way - 5 groups of PD importance ----
#top distinctiveness”, “high distinctiveness”,
#“moderate distinctiveness”, “low distinctiveness”, “very low distinctiveness”
#(ranked in the 5%, 25%, 50%, 75% and > 75% of the distribution, respectively).

# warning, quantiles ARE NOT percentage - they are equally separated data
hotspot_quantiles_loess <- cut(
  PD_actinos$std_resid_loess,
  breaks = quantile(PD_actinos$std_resid_loess, c(0, 0.05, 0.25, 0.5, 0.75, 1)),
  labels = c("lowest_100", "low_75", "moderate_50", "high_25", "top_5"),
  right  = FALSE,
  include.lowest = TRUE
)

PD_actinos<-cbind(PD_actinos, hotspot_quantiles_loess)

hotspot_quantiles_SR <- cut(
  PD_actinos$SR,
  breaks = quantile(unique(PD_actinos$SR), c(0, 0.05, 0.25, 0.5, 0.75, 1), names=T),
  labels = c("lowest_100", "low_75", "moderate_50", "high_25", "top_5"),
  right  = FALSE,
  include.lowest = TRUE
)
PD_actinos<-cbind(PD_actinos, hotspot_quantiles_SR)
PD_actinos<-PD_actinos%>%mutate(pd_sr_hots=if_else(hotspot_quantiles_loess == hotspot_quantiles_SR,
                                                   "congruent", "not_congruent"))
PD_actinos$congruence<- paste(PD_actinos$pd_sr_hots, PD_actinos$hotspot_quantiles_loess, sep="_")

## PCA between std_zscore pd and SR
for.pca.actinos <- PD_actinos%>%
  select(SR, std_resid_loess, PD_obs, std_mpd, std_mntd, mpd, mntd, grids )%>%
  st_drop_geometry()

rownames(for.pca.actinos)<-for.pca.actinos$grids
actinos_PCA<-prcomp(for.pca.actinos[,c(1:7)], center = TRUE, scale. = TRUE)
summary(actinos_PCA)
actinos_PCA$rotation
actinos_PCA
str(actinos_PCA)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(actinos_PCA, choices = 1:2, labels=rownames(for.pca.actinos))#labels=rownames(for.pca.actinos)
ggbiplot(actinos_PCA, choices = 1:2)#labels=rownames(for.pca.actinos)
ggbiplot(actinos_PCA, choices = 2:3)#labels=rownames(for.pca.actinos)

## DO THE SAME FOR ALL SPECIES GROUPS ----
# ith grid cells where all 4 groups are present only

## calculate percentage that falls into MPAs and proposed mpas 
# check the nmr of species spatial congruence

################################
## create phyloregions########
###########################
phylo_beta_grid3 <- phylobeta(sparse_grid3, phylo_fishes_baq,  index.family = "sorensen")
setdiff(new_phylo_fish_maker2$tip.label,colnames(sparse_grid3))
# needs a perfect match btw tips and spp in sparse data
#new_phylo_fish_maker_P <- keep.tip(new_phylo_fish_maker2, intersect(new_phylo_fish_maker2$tip.label, 
#                                                                    colnames(sparse_grid3)))

com.data.actinos<-match_phylo_comm(new_phylo_fish_maker2, sparse_grid3_actinos, delete_empty_rows = F)

### use from precedent calculations of ses PD 
## use community data without 1 or 2 species cells !!!
dim(comm_actinos_t_no1or2)
comm_actinos_no1or2_sparse<-t(dense2sparse(comm_actinos_t_no1or2))

phylo_beta_actinos <- phylobeta(com.data.actinos$comm, com.data.actinos$phy, index.family = "sorensen")
#phylo_beta_actinos <- phylobeta(comm_actinos_no1or2_sparse, com.data.actinos$phy)# 336 tips 

# simpson
class(phylo_beta_actinos[[1]])
names(phylo_beta_actinos[[1]])

library(ade4)
library(vegan)
# representer une matrice de distances avec un analyse en coordonnées principales PCoA
pcoa_sor_actino<-dudi.pco(phylo_beta_actinos[[3]], scannf=F, nf=3)
scatter(pcoa_sor_actino, posieig = "bottomright", pch=2)

#select the less distorting clustering method
select_linkage(phylo_beta_actinos[[1]])
#select optimal number of clusters with selected method
optim<-optimal_phyloregion(phylo_beta_actinos[[3]], k=30)
plot(optim$df$k, optim$df$ev, pch=20)# k - nbr of clusters VS explained variance given k

grid_aq_3_sp<-as_Spatial(grid_aq_3_PF)

#grid_aq_3_sp<-as_Spatial(grid_aq_3_PFwg)
#grid_aq_3_PFwg<-st_transform(grid_aq_3_PF, crs=4326)

proj4string(grid_aq_3_sp)

##estimate evolutionary distinctiveness of each phyloregion by computing the mean value of 
# phylogenetic beta diversity between a focal phyloregion and all other phyloregions in the study area.
y <- phyloregion(phylo_beta_actinos[[3]], shp=grid_aq_3_sp, k=8, method="average")
summary(y)
phylo_nmds<-y$NMDS

#make an sf file for
phyloreg_sf<-y$shp
st_crs(phyloreg_sf)
plot(phyloreg_sf)
actinos_phyloreg_sf<-st_as_sf(phyloreg_sf)
plot(actinos_phyloreg_sf["ED"])

names(y$region.df)
table(y$region.df$COLOURS)
#ed_palette<-c(unique(y$region.df$COLOURS))
PD_actinos<-PD_actinos%>%select(!c(cluster, ED, COLOURS))

PD_actinos<-merge(PD_actinos, y$region.df, by="grids", all.x=T, all.y=F)
names(PD_actinos)
write.csv2(PD_actinos, file="PD_actinos_all.csv")

PD_actinos %>% # Applying group_by & summarise
  group_by(cluster, .add=T) %>%
  summarise(mean_stdPD = mean(std_resid_loess), mean_stdPE= mean(std_PE), mean_MPD=mean(mpd),
            mean_ED=mean(ED), meanSR_fam=mean(sr_fam))
ggboxplot(PD_actinos, x="cluster", y="std_resid_loess")
ggboxplot(PD_actinos, x="cluster", y="std_PE")
ggboxplot(PD_actinos, x="cluster", y="mpd")
#is there any significat differnce btw clusters ?
library(rstatix)
res.kruskal <- PD_actinos %>% kruskal_test(std_resid_loess ~ cluster)
# or
kruskal.test(PD_actinos$std_resid_loess, g=PD_actinos$cluster)
#PD_actinos %>% kruskal_effsize(std_resid_loess ~ cluster)
# which groups are significantly different? -> Dunn test
pairwise.wilcox.test(PD_actinos$std_resid_loess, g = PD_actinos$cluster, p.adjust.method = "bonferroni")

phyloreg_actinos<-y$membership
phyloreg_actinos<-merge(phyloreg_actinos, actinos_phyloreg_sf, by="cluster", all.x=T, all.y=F)
head(phyloreg_actinos)

# colours for phylogerions
#"#9AB093" "#A5A6E2" "#B6A5B6" "#C0A58D" "#C79EC0" "#CA9AD4" "#D4A17A" "#D69BA5" 

tiff("ED_phylobeta_K8_actinos_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+geom_sf(data=actinos_phyloreg_sf, aes(fill=ED))+
  scale_fill_viridis(option="H", direction = 1, name="ED phyloregion", begin=0.4)+
  #scale_fill_manual(values=c("blue", "blue4",  "violet", "purple", "purple4", "yellow4", "yellow2", "yellow"))+ 
  theme_minimal()+
  #geom_sf(data=mpa_proposed, aes(colour=name), linewidth=0.5, fill=NA)+
  geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  #geom_sf(data=subset(actinos_gridPF, !is.na(ED_top5)),shape = 1, size=2, colour="red")+ 
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()  

phylo_nmds<-y$NMDS
phylo_nmds<-y$NMDS$points
class(phylo_nmds)
phylo_nmds<-as.data.frame(phylo_nmds)
phylo_nmds$phyloreg<-rownames(phylo_nmds)

tiff("nmds_phyloregions_actinos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
actino_nmds<-ggplot()+geom_point(data=phylo_nmds, aes(x=MDS1, y=MDS2, fill=phyloreg), pch=21, size=6)+
  scale_fill_brewer(palette = "Set3") +   # customized color palette
  labs(fill = "")+
  geom_text(data=phylo_nmds, aes(x=MDS1, y=MDS2, label=phyloreg), vjust=0.3, size=4)+
  theme_minimal()+
  theme(legend.position="none")
dev.off()

## create phylogen trees per cell (N>=3) and calculate ED inside -> meanED
phy_comm_1166 <-list()
for(i in rownames(sparse_grid3_no1_2)){
  print (i)
  flush.console()
  plot <- (sparse_grid3_no1_2[i, sparse_grid3_no1_2[i,] == 0])
  X <- names(plot)
  phyloX <- drop.tip(new_phylo_fish_maker, X, trim.internal = T, subtree=F, root.edge = 0) # sub-tree des esp du plot en question
  phy_comm_1166[[i]]<-phyloX
}
names(phy_comm_1166) <- rownames(sparse_grid3_no1_2)

ED_cells<-lapply(phy_comm_1166, evol.distinct, type="fair.proportion")

f <- function(x, i) {
  mean(x[i][, 2])
}
ED_mean_cells<-lapply(ED_cells, f)
ED_mean_cells<-do.call("rbind",ED_mean_cells)
ED_mean_cells<-as.data.frame(ED_mean_cells)
names(ED_mean_cells)<-"meanED"
ED_mean_cells$grids<-row.names(ED_mean_cells)

grid_aq3_sf<-merge(grid_aq3_sf, ED_mean_cells, by="grids", all.x=T, all.y=F)

ses_PD_phyloregions<-merge(ses_PD_phyloregions, ED_mean_cells,by="grids", all.x=T, all.y=F)
## correlation btw ses pd score and mean ed as alpha diversity
cor(ses_PD_phyloregions$meanED, ses_PD_phyloregions$zscore, use="complete.obs")#R2=0.8

## recalc ses_pd on cell >=3 species ----
new_phylo_fish_maker_no1_2 <- keep.tip(new_phylo_fish_maker2, intersect(new_phylo_fish_maker2$tip.label) )

                                                                        
###################                                                                                                                                             colnames(sparse_grid3_no1_2)))
## FAUX !!!!###### pas de la SR mais des abondances

## phylogenetic beta pd ses correction----
ses_Phylobeta_actinos<-phylobeta_ses(sparse_grid3_actinos, new_phylo_fish_maker_P, model="tipshuffle", reps=999,
                             index.family = "simpson")
## phylogerions with corrected ses pd (no need as Simpson index is independent of SR)
y <- phyloregion(ses_Phylobeta_actinos[[4]], shp=grid_aq_3_sp, k=15, method="average")
plot_NMDS(y, cex=3)
text_NMDS(y)
par(mar=rep(0,4))
plot(y, palette="viridis")  

