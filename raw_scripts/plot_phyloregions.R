
# ophi

#breaks <- c(0.7,0.85,1)

tiff("ED_K7_phylobeta_ophi.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ED_K7_phylobeta_ophi<-ggplot()+
  geom_sf(data=ophi_phyloreg_sf,aes(fill=ED))+
  scale_fill_viridis(option="H", direction = 1, name="ED", begin=0.4)+
  theme_minimal()+
  geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  theme(legend.text = element_text(size=6),legend.title=element_text(size=6),
        legend.key.width=unit(0.2, "cm"))+  coord_sf(datum = NA)
dev.off()  

tiff("phyloregions_ophi.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
phyloregions_ophi<-ggplot()+geom_sf(data=ophi_phyloreg_sf, aes(fill=as.factor(as.character(cluster))))+
  scale_fill_brewer(palette = "Dark2") +   # customized color palette
  theme_minimal()+
  labs(fill = "")+
  #geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  theme(legend.position = "none")+
  coord_sf(datum = NA)
dev.off()

# pycnos
tiff("ED_K7_phylobeta_pycnos_209sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ED_K7_phylobeta_pycnos<-ggplot()+geom_sf(data=pycno_phyloreg_sf, aes(fill=ED))+
  scale_fill_viridis(option="H", direction = 1, name="ED", begin=0.4)+
  #scale_fill_continuous(low="yellow1", high="red4", name="ED")+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  theme(legend.text = element_text(size=6),legend.title=element_text(size=6),
        legend.key.width=unit(0.2, "cm"))+  coord_sf(datum = NA)
dev.off()  

tiff("phyloregions_pycnos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
phyloregions_pycnos<-ggplot()+geom_sf(data=pycno_phyloreg_sf, aes(fill=as.factor(as.character(cluster))))+
  scale_fill_brewer(palette = "Dark2") +   # customized color palette
  theme_minimal()+
  labs(fill = "")+
  #geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  theme(legend.position = "none")+
  coord_sf(datum = NA)
dev.off()

# actinos

tiff("ED_phylobeta_K8_actinos_281sp.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
ED_phylobeta_K8_actinos<-ggplot()+geom_sf(data=actinos_phyloreg_sf, aes(fill=ED))+
  scale_fill_viridis(option="H", direction = 1, name="ED", begin=0.4)+
  theme_minimal()+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  theme(legend.text = element_text(size=6),legend.title=element_text(size=6),
        legend.key.width=unit(0.2, "cm"))+
  coord_sf(datum = NA)
dev.off()  

tiff("phyloregions_actinos.tiff",width=6,height=5,res=600, units="in", compression ="lzw")
phyloregions_actinos<-ggplot()+geom_sf(data=actinos_phyloreg_sf, aes(fill=as.factor(as.character(cluster))))+
  scale_fill_brewer(palette = "Set3") +   # customized color palette
  theme_minimal()+
  labs(fill = "")+
  #geom_sf(data=Coastline_sf, fill="white", linewidth=0.1)+
  geom_sf(data=mpa_eez, fill=NA, color="black", linewidth=0.2)+
  geom_sf(data=mpa_proposed, fill=NA, color="black", linewidth=0.2)+
  #theme(legend.position = "none")+
  coord_sf(datum = NA)
dev.off()

legend_1 <- get_legend(ED_phylobeta_K8_actinos)
legend_2 <- get_legend(ED_K7_phylobeta_pycnos)
legend_3 <- get_legend(iris_3)
legends <- ggarrange(legend_1, legend_2, legend_3, nrow=3)

library(ggpubr)#
tiff("ed_phyloreg_act_pycno_ophi.tiff",width=10,height=3,res=600, units="in", compression ="lzw")
ggarrange(ED_phylobeta_K8_actinos,ED_K7_phylobeta_pycnos,ED_K7_phylobeta_ophi,
          nrow = 1, ncol = 3,
          labels=c("(A)", "(B)", "(C)"), 
          vjust=3, hjust=-2,
          common.legend = F, legend="right",
          font.label=list(size=8, face="bold"))
dev.off()

tiff("phyloreg_nmds_act_pycno_ophi.tiff",width=10,height=7,res=600, units="in", compression ="lzw")
ggarrange(phyloregions_actinos, phyloregions_pycnos, phyloregions_ophi,
          actino_nmds, pycno_nmds, ophi_nmds,
          nrow = 2, ncol = 3,
          labels=c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"), 
          vjust=3, hjust=-2,
          common.legend = F,
          font.label=list(size=8, face="bold"))
dev.off()
