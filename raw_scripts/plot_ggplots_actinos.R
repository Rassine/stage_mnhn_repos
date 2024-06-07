###########
## PLOTS ##
###########
library(ggplot2)
library(ggeffects)
library(ggpubr) #for multiple lots together
#library(CCAMLRGIS)
library(viridis)
library(ggsn)
library(ggpattern)
library(measoshapes)
data(measo_regions05_coastline)
measo <- measo_regions05 %>% group_by(name) %>% summarize() %>% 
  inner_join(measo_names)
#> Joining, by = "name"
#coast <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% dplyr::filter(sovereignt == "Antarctica")
#coast <- st_transform(coast, st_crs(measo))
tiff("grid_cell_actinos.tiff",width=6,height=6,res=400, units="in", compression ="lzw")
ggplot()+
  #geom_sf(data=measo, aes(fill = fill)) + scale_fill_identity() +
  geom_sf(data=grid_aq_3_PF, fill = NA)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red3", size=7,linetype = "dashed")+
  theme_minimal()+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()  

Coastline=load_Coastline()
Coastline_sf<-st_as_sf(Coastline)
st_crs(Coastline_sf)=st_crs(grid_aq_3)

mpa_proposed<-st_transform(mpa_proposed, crs=prj_laea)
plot(st_geometry(mpa_proposed))

Coastline_sf<-st_transform(Coastline_sf, crs=prj_laea)

long_3d<-read_sf("3dg_longitude.shp",
                    layer="3dg_longitude")

# plot of occurrences used for inext
actinos_gridPF_comm<-actinos_gridPF_1%>%
  filter(newscientificname %in% rownames(comm_actinos_t))%>%
  filter(grids %in% colnames(comm_actinos_t))
length(unique(actinos_gridPF_comm$newscientificname))
length(unique(actinos_gridPF_comm$grids))

tiff("actinos_occ_281sp.tiff",width=6,height=6,res=400, units="in", compression ="lzw")
actinos_occ_281sp<-ggplot()+
  geom_sf(data=actinos_gridPF_1, shape = 16, size=2, alpha = 0.4, color = "purple4" )+
  geom_jitter()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="red2", size=6,linetype = "dashed")+
  theme_minimal()+
  coord_sf(datum = NA)
dev.off()
#####################################################################

## species richness
ggplot()+
  geom_sf(data=subset(PD_actinos,loess_PD>-6), aes(fill=loess_PD))+
  scale_fill_viridis(option="A", direction = 1, begin = 0.3)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, alpha=0.2)
#####################################################################

ggplot()+
  geom_sf(data=domains_sf, aes(fill=pd.obs.z.y))+
  scale_fill_viridis(option="E", direction = 1)+
  geom_sf(data=Coastline_sf, fill="white", alpha=0.5)
geom_sf(data=mpa_eez, fill="blue", alpha=0.2)

#####################################################################
## PD corrected plots ----
table(PD_actinos$SR)
PD_actinos_no1or2<-subset(PD_actinos, SR>2 )

tiff("PD_actinos_stdPD.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
stdPD_actinos<-ggplot()+  
  geom_sf(data=PD_actinos, aes(fill=std_resid_loess), color=NA)+
  #scale_fill_continuous(low="blue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="stdPD")+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme_minimal()+
  theme(legend.text = element_text(size=6),legend.title=element_text(size=6),
        legend.key.width=unit(0.2, "cm"))+
  coord_sf(datum = NA)
dev.off()

ggplot()+
  geom_sf(data=PD_actinos, aes(fill=PE), color=NA)+
  #scale_fill_continuous(low="blue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = -1, name="PE")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
    theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  theme(legend.text = element_text(size=8),legend.title=element_text(size=8))+
  coord_sf(datum = NA)

tiff("std_PE_actinos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
std_PE_actinos<-ggplot()+
  geom_sf(data=PD_actinos, aes(fill=std_PE), color=NA)+
  #scale_fill_continuous(low="blue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="stdPE")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  theme(legend.text = element_text(size=8),legend.title=element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("MPD_actinos_281.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
MPD_actinos<-ggplot()+
  geom_sf(data=PD_actinos, aes(fill=mpd), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="MPD")+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme_minimal()+
  theme(legend.text = element_text(size=8),legend.title=element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("PD_actinos_SR_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
PD_actinos_SR<-ggplot()+
  geom_sf(data=PD_actinos, aes(fill=SR), color=NA)+
  scale_fill_viridis(option="H", direction = 1, name="SR")+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  theme_minimal()+
  theme(legend.text = element_text(size=8),legend.title=element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

ggplot()+
  geom_sf(data=PD_actinos, aes(fill=sr_fam), color=NA)+
  scale_fill_viridis(option="H", direction = 1, name="SR_fam")+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  theme_minimal()+
  theme(legend.text = element_text(size=8),legend.title=element_text(size=8))+
  coord_sf(datum = NA)


tiff("PD_actinos_std_loess_hotspot_coldspot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_stdPDloess)), fill="#660000", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_stdPDloess)), fill="#000066", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot5_stdPDloess)), fill="#993333", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot5_stdPDloess)), fill="#000099", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot2.5_stdPDloess)), fill="#CC6666", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot2.5_stdPDloess)), fill="#0000FF", color="NA", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="#FFFFFF")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

actinos_stdPD_hot_cold<-ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_stdPDloess)), fill="red3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_stdPDloess)), fill="dodgerblue3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos, is.na(hotspot10_stdPDloess)& is.na(coldspot10_stdPDloess)), 
          fill="grey77", color="NA", linewidth=0.5)+theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)

tiff("MPD_actinos_hotspot_coldspot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_mpd)), fill="#660000", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_mpd)), fill="#000066", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot5_mpd)), fill="#993333", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot5_mpd)), fill="#000099", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot2.5_mpd)), fill="#CC6666", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot2.5_mpd)), fill="#0000FF", color="NA", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="#FFFFFF")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

actinos_MPD_hot_coldsplot<-
ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_mpd)), fill="red3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_mpd)), fill="dodgerblue3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos, is.na(hotspot10_mpd)& is.na(coldspot10_mpd)), 
          fill="grey77", color="NA", linewidth=0.5)+theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)

tiff("MNTD_actinos_hotspot_coldspot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+

  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_mntd)), fill="#660000", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_mntd)), fill="#000066", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot5_mntd)), fill="#993333", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot5_mntd)), fill="#000099", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot2.5_mntd)), fill="#CC6666", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot2.5_mntd)), fill="#0000FF", color="NA", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.4)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.4)+
  geom_sf(data=Coastline_sf, fill="#FFFFFF")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("SR_actinos_hotspot_coldspot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_sr)), fill="#660000", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_sr)), fill="#000066", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot5_sr)), fill="#993333", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot5_sr)), fill="#000099", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot2.5_sr)), fill="#CC6666", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot2.5_sr)), fill="#0000FF", color="NA", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="white")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("stdPE_actinos_hotspot_coldspot.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_stdPE)), fill="#660000", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_stdPE)), fill="#000066", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot5_stdPE)), fill="#993333", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot5_stdPE)), fill="#000099", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot2.5_stdPE)), fill="#CC6666", color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot2.5_stdPE)), fill="#0000FF", color="NA", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="#FFFFFF")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

actinos_stdPE_hot_coldsplot<-ggplot()+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_stdPE)), fill="red3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos,!is.na(coldspot10_stdPE)), fill="dodgerblue3",color="NA", linewidth=0.5)+
  geom_sf(data=subset(PD_actinos, is.na(hotspot10_stdPE)& is.na(coldspot10_stdPE)), 
          fill="grey77", color="NA", linewidth=0.5)+theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)

library(ggeffects)
library(ggpubr)#
tiff("actinos_hotcold.tiff",width=10,height=4,res=600, units="in", compression ="lzw")
ggarrange(actinos_stdPD_hot_cold,actinos_stdPE_hot_coldsplot, actinos_MPD_hot_coldsplot,
          nrow = 1, ncol = 3,
          labels=c("(A)", "(B)", "(C)"), 
          vjust=7, hjust=-1,
          common.legend = FALSE, legend="right",
          font.label=list(size=10, face="bold"))
dev.off()


tiff("PD_complementarity_actinos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
PD_complementarity_actino<-ggplot()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(PD_actinos,!is.na(pd_added)), aes(fill=pd_added), color=NA)+
  scale_fill_viridis(option="H", direction = 1, begin = 0.2, name = "PD ")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("SR_complementarity_actinos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
SR_complementarity_actino<-ggplot()+
  #scale_fill_continuous(low="lightblue", high="red", name="PD")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=subset(PD_actinos,!is.na(sr_added)), aes(fill=sr_added), color=NA)+
  scale_fill_viridis(option="H", direction = 1, begin = 0.2, name = "SR ")+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()


## with iNEXT SR estimate
tiff("PD_actinos_SRrare_15_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
PD_actinos_SRrare_15_281sp<-ggplot()+
  geom_sf(data=PD_actinos_inexTD, aes(fill=qD_15), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="SRr")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="white")+
  theme_minimal()+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

tiff("PD_actinos_stdPDrare_15_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=PD_actinos_inexPD, aes(fill=loess_PDrare_15), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="std PDr")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="deeppink", linewidth=0.6)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  #geom_sf(data=subset(actinos_gridPF_1, !is.na(ED_tot_PDr)),shape = 1, size=1, color="violet")+ 
  geom_sf(data=mpa_proposed, fill=NA, color="orange", linewidth=0.6)+
  geom_sf(data=Coastline_sf, fill="white")+
  theme_minimal()+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()


## plots with nbr of rarefied and extrap
tiff("PD_actinos_RareExtrap_15_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=PD_actinos_inexPD, aes(fill=method_15), color=NA)+
  scale_fill_manual(values=c("lightblue","orange","deeppink"), name="")+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  #scale_fill_viridis(option="D", direction = -1, name="std PDr")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  #geom_sf(data=subset(PD_actinos,!is.na(hotspot10_SR)), fill=NA, color="red", linewidth=1)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="white")+
  theme_minimal()+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()

# get a global comined plot ----
library(ggeffects)
library(ggpubr)#
tiff("actinos_SR_SRr_PD_PDr_hotcold_compl.tiff",width=8,height=10,res=600, units="in", compression ="lzw")
ggarrange(PD_actinos_SR_281sp,PD_actinos_SRrare_15_281sp,
          PD_actinos_std_loess_hotspot_281sp, PD_actinos_PDrare_15_10hotspot_281sp,
          PD_complementarity_actino, PDrare_actinos_std_loess_hotspot_coldspot,
          nrow = 3, ncol = 2,
          labels=c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"), 
          vjust=3, hjust=-2,
          common.legend = FALSE, legend="right",
          font.label=list(size=10, face="bold"))
dev.off()

tiff("actinos_SR_SRr_PD_PDr_hotcold.tiff",width=8,height=10,res=600, units="in", compression ="lzw")
ggarrange(PD_actinos_SR_281sp,PD_actinos_std_loess_hotspot_281sp,
          PD_actinos_SRrare_15_281sp,
          PD_actinos_PDrare_15_10hotspot_281sp,
          SRrare_actinos_std_loess_hotspot_coldspot, PDrare_actinos_std_loess_hotspot_coldspot,
          nrow = 3, ncol = 2,
          labels=c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"), 
          vjust=3, hjust=-2,
          common.legend = FALSE, legend="right",
          font.label=list(size=10, face="bold"))
dev.off()

## plot with PD bias corrected by samplig effort in grid cells
cor.test(PD_actinos$richness, PD_actinos$PD_obs, method="pearson")

tiff("PD_SR_congruence.tiff",width=6,height=5,res=600,
     units="in", compression ="lzw")
ggplot()+
  #geom_sf(data=PD_actinos, aes(fill=congruence), color=NA)+
  geom_sf(data=grid_aq_3_PF, fill=NA)+
  #scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.8)+
  geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkblue", size=1,linetype = "dashed", linewidth=0.7)+
  #geom_sf(data=subset(actinos_gridPF, !is.na(ED_top5)),shape = 1, size=2, colour="red")+ 
  geom_sf(data=mpa_proposed, fill=NA, color="blue3", linewidth=0.8)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot_quantiles_loess=="top_5")), fill="red3")+
  geom_sf(data=subset(PD_actinos, hotspot_quantiles_SR=="top_5"), fill="purple")+
  geom_sf(data=subset(PD_actinos,congruence=="congruent_top_5"), fill="pink")+
  theme_bw()+
  geom_sf(data=Coastline_sf, fill="white")+
  coord_sf(datum = NA)
dev.off()

tiff("actino_measures_hot_congruence.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=hots_congr_actinos, aes(fill=hot_congr), color=NA)+
  scale_fill_manual(values=c("white","lightblue","orange","deeppink"), name="")+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  #scale_fill_viridis(option="D", direction = -1, name="std PDr")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.5)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  #geom_sf(data=subset(PD_actinos,!is.na(hotspot10_SR)), fill=NA, color="red", linewidth=1)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.5)+
  geom_sf(data=Coastline_sf, fill="white")+
  theme_minimal()+
  theme(legend.text = element_text(size=8))+
  coord_sf(datum = NA)
dev.off()


tiff("PD_actinos_PEW.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=PD_actinos, aes(fill=PEW), color=NA)+
  #scale_fill_continuous(low="blue", high="red", name="Std PD")+
  scale_fill_viridis(option="H", direction = 1, name="PE")+
  geom_sf(data=Coastline_sf, fill="white")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="#FF61CC", linewidth=0.6)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  #geom_sf(data=subset(actinos_gridPF, !is.na(ED_top2.5)),shape = 19, size=2)+ 
  geom_sf(data=mpa_proposed, fill=NA, color="orange", linewidth=0.6)+
  #geom_sf(data=subset(PD_actinos,!is.na(hotspot5_std_zscore)), fill=NA, color="red", linewidth=1)+
  theme_minimal()+
  coord_sf(datum = NA)
dev.off()

ggplot()+
  #geom_sf(data=PD_actinos, aes(fill=congruence), color=NA)+
  geom_sf(data=grid_aq_3_PF, linewidth=0.4)+
  #scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.8)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkblue", size=1,linetype = "dashed", linewidth=0.7)+
  geom_sf(data=mpa_proposed, fill=NA, color="blue3", linewidth=0.8)+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_resid_loess)), fill="red3")+
  geom_sf(data=subset(PD_actinos,!is.na(hotspot10_SR)), colour="yellow", fill=NA,linewidth=0.8)+
  #geom_sf(data=subset(PD_actinos, hotspot_quantiles_SR=="top_5"), fill="purple")+
  #geom_sf(data=subset(PD_actinos,congruence=="congruent_top_5"), fill="pink")+
  #geom_sf(data=subset(actinos_gridPF, !is.na(ED_top2.5)),shape = 19, size=2)+ 
  geom_sf(data=subset(actinos_gridPF, !is.na(ED_top5)),shape = 19, size=2)+ 
  
  theme_bw()+
  #geom_sf(data=Coastline_sf, fill="white")+
  coord_sf(datum = NA)

## plots for rarefied PD ----
PD_actinos_100
SR_hill_100<-ggplot()+
  geom_sf(data=PD_actinos_100, aes(fill=est_asy), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=Coastline_sf, fill="white")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="gray", size=0.5)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  theme_minimal()

SR_norich<-ggplot()+
  geom_sf(data=PD_actinos_100, aes(fill=SR), color=NA)+
  #scale_fill_continuous(low="lightblue", high="red", name="Std PD")+
  scale_fill_viridis(option="D", direction = 1, name="Std PD")+
  geom_sf(data=Coastline_sf, fill="white")+
  #geom_sf(data=domains_ses_pd, fill=NA,size=0.7, color="black")+
  geom_sf(data=mpa_eez, fill=NA, color="gray", size=0.5)+
  #geom_sf(data=subset(aq_fronts, NAME == "Polar Front (PF)"), color="darkgrey", size=0.7,linetype = "dashed")+
  #geom_line(aes(st_coordinates(subset(aq_fronts, NAME == "Polar Front (PF)"))))+
  theme_minimal()


dev.off()

#####################################################################
## add cells with significant ses pd values + hotspots
ggplot()+
  geom_sf(data=subset(grid_aq_3,pd.obs.p>0.95), fill=NA, color="red", size=0.8)

#####################################################################
tiff("SES_PD_phyloreg.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ggplot()+
  geom_sf(data=subset(grid_aq3_sf,!is.na(richness)), aes(fill=zscore))+
  scale_fill_viridis(option="A", direction = 1)+
  geom_sf(data=Coastline_sf, fill="white")+
  geom_sf(data=mpa_eez, fill=NA, color="black", size=0.7)
dev.off()

