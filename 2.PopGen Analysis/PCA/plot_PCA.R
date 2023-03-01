library(tidyverse)
library(plotly)
library(ggrepel)
library(ggforce)

#load the data
d <- read.table('1240v2_HO_NP_final_MDH_laza2022_eurasie_focus_caucase.evec', as.is=TRUE)
colnames(d)=c('ID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Group_ID')
ind_Anc <- read.table('metadata/ancientind_clean.txt', header=T, as.is=TRUE)
ind_eurasie <- read.table('metadata/eurasie_noEastAsia.txt', header=T, as.is=TRUE)
group_complet <- read.table('metadata/group_region_time.txt', header = TRUE)

#Subset the data
Anc = merge(ind_Anc,d)
eurasie = merge(ind_eurasie,d)
MT = subset(d, Group_ID == 'Mentesh')
Anc1 = subset(Anc, (Region != "SouthAsia" & Region != "Baikal" & Region != "EasternStepppe")  & (Time == "Chalcolithic" | Time == "Mesolithic" |  Time == "Paleolithic"))
Anc2 = subset(Anc, (Region != "SouthAsia" & Region != "Baikal" & Region != "EasternStepppe") & 
                (Time == "Neolithic" | Time == "EMBA" | Time == "MLBA" | Time == "IronAge" | Time == "Historical"))


#calculate each PC %
axe = read.table("1240v2_HO_NP_final_MDH_laza2022_eurasie_focus_caucase.eval")
axe$PC <- seq.int(nrow(axe))
axe$PCnb <- paste("PC",axe$PC)
sum = sum(axe$V1)
PC = 100*axe[1]/sum
PC1 = PC[1,]
PC2 = PC[2,]
PC3 = PC[3,]
PC4 = PC[4,]
png("PC.png")
ggplot(axe, aes(PC,V1, label=PCnb))+
  geom_point()+
  coord_trans(x = "log10")+geom_line()+
  geom_text(nudge_y = 1, nudge_x = 1)+
  labs(y="eval")
dev.off()

#plot the PCA for PC1, PC2 and PC3
pdf("ACP_eurasie.pdf", width=12, height=10)
PCA12 <- ggplot()+
  geom_text(subset(eurasie), mapping = aes(PC1,PC2, label = Group_ID, text = ID), col = 'snow3', size = 2)+
  geom_point(Anc2, mapping = aes(PC1,PC2, fill = Region, col = Region, shape = Time,text = paste(Group_ID, ':',ID)), 
             size=3,alpha = 0.5)+
  geom_point(Anc1, mapping = aes(PC1,PC2, col = Region, shape = Time, text = paste(Group_ID, ':',ID)), 
             size=3, alpha = 0.5)+
  geom_point(subset(Anc, Region == 'Caucasus'), 
             mapping = aes(PC1,PC2, shape = Time, text = paste(Group_ID, ':',ID)), col = 'hotpink4', fill = 'hotpink', size = 4)+
  geom_point(MT, mapping = aes(PC1,PC2, text = paste(Group_ID, ':',ID)), col='black',fill = 'hotpink3',shape = 21, size = 5)+
  geom_label_repel(MT, mapping = aes(PC1,PC2, label = ID), col='hotpink3', size =6)+
  theme_light(base_size = 20)+
  labs(x=paste("PC1=",round(PC1, digits = 1),"%"),y=paste("PC2=",round(PC2, digits = 1),"%"))+
  scale_shape_manual(values = c(3,24,23,22,4,25,21,13))+
  scale_fill_manual(values = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                               "CentralSteppe"="gold","Europe"="darkslategray4",
                               "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                               "Turan"="limegreen","WesternSteppe"="tomato", "Mesopotamia" = "brown"))+
  scale_color_manual(values = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                                "CentralSteppe"="gold","Europe"="darkslategray4",
                                "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                                "Turan"="limegreen","WesternSteppe"="tomato", "Mesopotamia" = "brown"), guide = FALSE)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  coord_fixed(ratio = 1)

PCA32 <- ggplot()+
  geom_text(subset(eurasie), mapping = aes(PC3,PC2, label = Group_ID, text = ID), col = 'snow3',size = 2)+
  geom_point(Anc2, mapping = aes(PC3,PC2, fill = Region, col = Region, shape = Time, text = paste(Group_ID, ':',ID)), 
             size=3, alpha = 0.5)+
  geom_point(Anc1, mapping = aes(PC3,PC2, col = Region, shape = Time, text = paste(Group_ID, ':',ID)),
             size=3, alpha  = 0.5)+
  geom_point(subset(Anc, Region == 'Caucasus'), mapping = aes(PC3,PC2, shape = Time,text = paste(Group_ID, ':',ID)), 
             col = 'hotpink4', fill = 'hotpink', size = 5)+
  geom_point(MT, mapping = aes(PC3,PC2,text = paste(Group_ID, ':',ID)), col='black',fill = 'hotpink3',shape = 21, size = 6)+
  geom_label_repel(MT, mapping = aes(PC3,PC2, label = ID, text = ID), col='hotpink3', size =7)+
  theme_light(base_size = 20)+
  labs(x=paste("PC3=",round(PC3, digits = 1),"%"),y=paste("PC2=",round(PC2, digits = 1),"%"))+
  scale_shape_manual(values = c(3,24,23,22,4,25,21,13))+
  scale_fill_manual(values = c("Anatolia"="darkolivegreen3","Caucase"="hotpink",
                               "CentralSteppe"="gold","Europe"="darkslategray4",
                               "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                               "Turan"="limegreen","WesternSteppe"="tomato","Mesopotamia" = "brown"))+
  scale_color_manual(values = c("Anatolia"="darkolivegreen3","Caucase"="hotpink",
                                "CentralSteppe"="gold","Europe"="darkslategray4",
                                "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                                "Turan"="limegreen","WesternSteppe"="tomato","Mesopotamia" = "brown"), guide = FALSE)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  coord_fixed(ratio = 1)
PCA12
PCA32
dev.off()

#plot interactively 
ggplotly(PCA12)
ggplotly(PCA32)


p <- plot_ly(Anc, 
             x = ~PC1, y = ~PC2, z = ~PC3,
             symbol = ~Time,
             symbols = c('Chalcolithic'= 'circle','EMBA'='square-open','Historical'='diamond',
                         'IronAge'='diamond-open','Mesolithic'='x','MLBA'='square','Modern'='cross',
                         'Neolithic'='circle','Paleolithic'='x'),
             color = ~Region,
             colors = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                        "CentralSteppe"="gold","Europe"="darkslategray4",
                        "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="coral2","SouthAsia"="brown","Mesopotamia"="brown",
                        "Turan"="limegreen","WesternSteppe"="purple"),
             
             type="scatter3d",
             mode = 'markers', text = ~paste(Group_ID, ID), size = I(60)
)
p1 = add_markers(p, data = eurasie,
            x = ~PC1, y = ~PC2, z = ~PC3, text = ~Group_ID, color = I('grey80'), symbol = I('cross'), size = I(30))

add_markers(p1, data = MT,
                 x = ~PC1, y = ~PC2, z = ~PC3, text = ~ID, color = I('red'), symbol = I('circle'), size = I(100))


###zoom on Caucasus
Caucasus_classif = read.table(file = "Caucase_sites.txt", header = T)

Caucasus = merge(Caucase_classif,d)

zoom <- ggplot()+
  geom_text(subset(eurasie), mapping = aes(PC1,PC2, label = Group_ID, text = ID), col = 'snow3', size = 2)+
  geom_point(Anc2, mapping = aes(PC1,PC2, fill = Region, col = Region, shape = Time,text = paste(Group_ID, ':',ID)), 
             size=3,alpha = 0.5)+
  geom_point(Anc1, mapping = aes(PC1,PC2, col = Region, shape = Time, text = paste(Group_ID, ':',ID)), 
             size=3, alpha = 0.5)+
  geom_point(Caucasus, 
             mapping = aes(PC1,PC2, shape = Time, text = paste(Group_ID, ':',ID),
                           #col = C_Region,
                           ), fill = 'hotpink',  
             col = 'hotpink4',
             size = 4, stroke = 1)+
  geom_point(MT, mapping = aes(PC1,PC2, text = paste(Group_ID, ':',ID)), col='black',fill = 'hotpink3',shape = 21, size = 5)+
  #geom_label_repel(MT, mapping = aes(PC1,PC2, label = ID), col='hotpink3', size =6)+
  geom_label_repel(subset(Anc, Region == 'Caucasus' & Time == 'Neolithic'),
                   mapping = aes(PC1,PC2, label = ID), col = 'hotpink', size =4)+
  theme_light(base_size = 10)+
  labs(x=paste("PC1=",PC1,"%"),y=paste("PC2=",PC2,"%"))+
  scale_shape_manual(values = c(3,24,23,22,4,25,21,13))+
  scale_fill_manual(values = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                               "CentralSteppe"="gold","Europe"="darkslategray4",
                               "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                               "Turan"="limegreen","WesternSteppe"="purple", "Mesopotamia"="brown"))+
  scale_color_manual(values = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                                "CentralSteppe"="gold","Europe"="darkslategray4",
                                "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown", "Mesopotamia"="brown",
                                "Turan"="limegreen","WesternSteppe"="purple", 'NorthCaucase' = 'cornflowerblue', 'SouthCaucase' = 'chartreuse3'), guide = FALSE)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  xlim(-0.04, 0) + ylim(-0.04, 0.025)+
  coord_fixed(ratio = 1)

pdf('zoom_Caucasus.pdf') 
zoom
dev.off()
