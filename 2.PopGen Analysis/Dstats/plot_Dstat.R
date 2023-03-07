library(tidyverse)

#Supplementary Fig 3

D = read.table('out_Dstat_MT_Paleo_Neo_select', col.names = c("file", "Outgroup","Group_ID", "MT",
                                                                 "test", "Dstat", "Z","ABBA","ABAB","SNPs"))
D[1] = NULL
popgroup = read.table('metadata/group_region_time_all.txt', header = T)

D_group = merge(D, popgroup)
anti_join(D, D_group)
stderr = D_group$Dstat/D_group$Z
D_group$Dmin = D_group$Dstat - 3*stderr
D_group$Dmax = D_group$Dstat + 3*stderr
D_group$p = 2*pnorm(D_group$Z, mean = 0, sd = 1, lower.tail=FALSE) # mulitply per 2 cause Dstat is a two-tailed test
D_group$padj = p.adjust(D_group$p, method = "BH")
D_group$Zadj = qnorm(0.5*D_group$padj, lower.tail = FALSE) # divise per 2 cause Dstat is a two-tailed test

Dgroup = D_group %>% mutate_all(funs(str_replace(.,'Luxembourg_Loschbour_published.DG','Loschbour')))  %>%
  mutate_all(funs(str_replace(.,'Belgium_UP_GoyetQ-2_published','GoyetQ2'))) %>%
  mutate_all(funs(str_replace(.,'Italy_North_Villabruna_HG','Villabruna'))) %>%
  mutate_all(funs(str_replace(.,'Iberia_HG_published','Iberia_HG'))) %>%
  mutate_all(funs(str_replace(.,'Russia_Caucasus_Eneolithic','NCaucasus_EN'))) %>%
  mutate_all(funs(str_replace(.,'Russia_Caucasus_EBA_Yamnaya','NCaucasus_Yamnaya')))%>%
  mutate_all(funs(str_replace(.,'Caucasus_lowlands_LC','SCaucasus_LC')))%>%
  mutate_all(funs(str_replace(.,'Serbia_Mesolithic_IronGates','IronGates_M')))%>%
  mutate_all(funs(str_replace(.,'Greece_Peloponnese_N','Greece_N')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_Epipaleolithic','Anatolia_EP')))%>%
  mutate_all(funs(str_replace(.,'ARM_Aknashen_N','Aknashen_N')))%>%
  mutate_all(funs(str_replace(.,'ARM_Masis_Blur_N','Masis_Blur_N')))%>%
  mutate_all(funs(str_replace(.,'TUR_SE_Mardin_PPN','Mardin_PPN')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_TellKurdu_EC','TellKurdu_EC')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_TepecikCiftlik_N.SG','TepecikCiftlik_N')))
Dgroup$Z = as.numeric(Dgroup$Z)
Dgroup$Zadj = as.numeric(Dgroup$Zadj)
Dgroup$Dstat = as.numeric(Dgroup$Dstat)
Dgroup$Dmin = as.numeric(Dgroup$Dmin)
Dgroup$Dmax = as.numeric(Dgroup$Dmax)
Dgroup$SNPs = as.numeric(Dgroup$SNPs)

Dgroup$Group_ID <- factor(Dgroup$Group_ID, levels=c("EEHG",'WSHG','NCaucasus_Yamnaya',
                                                    'Iran_GanjDareh_N','Iran_LN_SehGabi','IRQ_Nemrik9_PPN','PPN','Israel_C','Latvia_HG',
                                                    'Iberia_HG','Loschbour','GoyetQ2',
                                                    'Villabruna','IronGates_M','Serbia_EN','Greece_N',
                                                    'Barcin_N','Anatolia_EP','CHG','Masis_Blur_N','Aknashen_N',
                                                    'NCaucasus_EN','Armenia_C',
                                                    'SCaucasus_LC', 'Polutepe_002','Mentesh_001','MT23','MT26','MT7'))

Dgroup$test <- factor(Dgroup$test, levels=c('CHG','Aknashen_N','Masis_Blur_N','NCaucasus_EN','Polutepe_002','SCaucasus_LC',
                                            'Armenia_C','IRQ_Nemrik9_PPN','IRQ_Shanidar_N','IRQ_Bestansur_PPN',
                                            'Mardin_PPN','Barcin_N','TellKurdu_EC','TepecikCiftlik_N','MT23','MT26','MT7'))

plot_select <- ggplot(Dgroup[Dgroup$test != 'MT23' & Dgroup$test != 'MT26' & Dgroup$test != 'MT7',],
                            aes(Group_ID, Dstat, ymin=Dmin, ymax=Dmax, size = SNPs, col = Region),)+
  geom_pointrange(shape = ifelse(abs(Dgroup[Dgroup$test != 'MT23' & Dgroup$test != 'MT26' & Dgroup$test != 'MT7',]$Z)>3,16,1))+
  scale_colour_manual(values = c("Anatolia"="darkolivegreen3","Caucase"="hotpink",
                                 "CentralSteppe"="gold","Europe"="darkslategray4",
                                 "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                                 "Turan"="limegreen","WesternSteppe"="purple"))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.position = "right")+
  scale_size_binned_area(breaks = c(10000), max_size = 2)+
  geom_hline(yintercept=0, col = "black",linetype = 3, size =1)+
  facet_grid(test~MT)+
  ggtitle('D(Mbuti, X; Y, MT)')
plot_select

pdf("Dstat_corrected_MT_samples_select_final.pdf", width=11, height=22)
plot_select
dev.off()

############################
# Figure 2.c #Dstat with mentesh as a group

D = read.table('out_Dstat_Mentesh_Paleo_Neo_select', col.names = c("file", "Outgroup","Group_ID", "MT",
                                                                          "test", "Dstat", "Z","ABBA","ABAB","SNPs"))
D[1] = NULL
popgroup = read.table('group_region_time_all.txt', header = T)

D_group = merge(D, popgroup)
anti_join(D, D_group)
stderr = D_group$Dstat/D_group$Z
D_group$Dmin = D_group$Dstat - 3*stderr
D_group$Dmax = D_group$Dstat + 3*stderr
D_group$p = 2*pnorm(abs(D_group$Z), mean = 0, sd = 1, lower.tail=FALSE) # mulitply per 2 cause Dstat is a two-tailed test
D_group$padj = p.adjust(D_group$p, method = "BH")
D_group$Zadj = ifelse(D_group$padj !=0,
                qnorm(0.5*D_group$padj, lower.tail = FALSE), # divise per 2 cause Dstat is a two-tailed test
                abs(D_group$Z))

Dgroup = D_group %>% mutate_all(funs(str_replace(.,'Luxembourg_Loschbour_published.DG','Loschbour')))  %>%
  mutate_all(funs(str_replace(.,'Belgium_UP_GoyetQ-2_published','GoyetQ2'))) %>%
  mutate_all(funs(str_replace(.,'Italy_North_Villabruna_HG','Villabruna'))) %>%
  mutate_all(funs(str_replace(.,'Iberia_HG_published','Iberia_HG'))) %>%
  mutate_all(funs(str_replace(.,'Russia_Caucasus_Eneolithic','NCaucasus_EN'))) %>%
  mutate_all(funs(str_replace(.,'Russia_Caucasus_EBA_Yamnaya','NCaucasus_Yamnaya')))%>%
  mutate_all(funs(str_replace(.,'Caucasus_lowlands_LC','SCaucasus_LC')))%>%
  mutate_all(funs(str_replace(.,'Serbia_Mesolithic_IronGates','IronGates_M')))%>%
  mutate_all(funs(str_replace(.,'Greece_Peloponnese_N','Greece_N')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_Epipaleolithic','Anatolia_EP')))%>%
  mutate_all(funs(str_replace(.,'ARM_Aknashen_N','Aknashen_N')))%>%
  mutate_all(funs(str_replace(.,'ARM_Masis_Blur_N','Masis_Blur_N')))%>%
  mutate_all(funs(str_replace(.,'TUR_SE_Mardin_PPN','Mardin_PPN')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_TellKurdu_EC','TellKurdu_EC')))%>%
  mutate_all(funs(str_replace(.,'Anatolia_TepecikCiftlik_N.SG','TepecikCiftlik_N')))
Dgroup$Z = as.numeric(Dgroup$Z)
Dgroup$Zadj = as.numeric(Dgroup$Zadj)
Dgroup$Dstat = as.numeric(Dgroup$Dstat)
Dgroup$Dmin = as.numeric(Dgroup$Dmin)
Dgroup$Dmax = as.numeric(Dgroup$Dmax)
Dgroup$SNPs = as.numeric(Dgroup$SNPs)

Dgroup$Group_ID <- factor(Dgroup$Group_ID, levels=c("EEHG",'WSHG','NCaucasus_Yamnaya',
                                                    'Iran_GanjDareh_N','Iran_LN_SehGabi','IRQ_Nemrik9_PPN','PPN','Israel_C','Latvia_HG',
                                                    'Iberia_HG','Loschbour','GoyetQ2',
                                                    'Villabruna','IronGates_M','Serbia_EN','Greece_N',
                                                    'Barcin_N','Anatolia_EP','CHG','Masis_Blur_N','Aknashen_N',
                                                    'NCaucasus_EN','Armenia_C',
                                                    'SCaucasus_LC', 'Polutepe_002'))


Dgroup$test <- factor(Dgroup$test, levels=c('CHG','Aknashen_N','Masis_Blur_N','NCaucasus_EN','Polutepe_002','SCaucasus_LC',
                                            'Armenia_C','IRQ_Nemrik9_PPN','IRQ_Shanidar_N','IRQ_Bestansur_PPN',
                                            'Mardin_PPN','Barcin_N','TellKurdu_EC','TepecikCiftlik_N'))

plot_select <- ggplot(Dgroup[Dgroup$test != 'IRQ_Bestansur_PPN',],aes(Group_ID, -Dstat, ymin=-Dmin, ymax=-Dmax, size = SNPs, col = Region),)+
  geom_pointrange(shape = ifelse(abs(Dgroup[Dgroup$test != 'IRQ_Bestansur_PPN',]$Zadj)>2,16,1))+
  scale_colour_manual(values = c("Anatolia"="darkolivegreen3","Caucasus"="hotpink",
                                 "CentralSteppe"="gold","Europe"="darkslategray4",
                                 "Iran"="forestgreen","Levant"="mistyrose2","Steppe"="red","SouthAsia"="brown",
                                 "Turan"="limegreen","WesternSteppe"="purple","Mesopotamia" = "brown"))+
  coord_flip(ylim=c(-0.1,0.1))+
  theme_bw(base_size = 14)+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 13),
        strip.background = element_rect(fill = 'white'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1)
        )+
  scale_size_binned_area(breaks = c(10000), max_size = 2)+
  geom_hline(yintercept=0, col = "black",linetype = 3, size =1)+
  facet_grid(.~test)+
  ylab("D(Mbuti, Y; Z, Mentesh)")
  #ggtitle('D(Mbuti, Y; Z, Mentesh)')+xlab('Y')
plot_select

pdf("Figure2c.pdf", width=20, height=6)
plot_select
dev.off()


###################################
## Supp Figure 6 # Mentesh group in caucase and anat BA

D = read.table('out_Dstat_mentesh_caucase_BA', col.names = c("file", "Outgroup","MT", "Group_ID",
                                                        "Caucase2", "Dstat", "Z","ABBA","ABAB","SNPs"))
D[1] = NULL
popgroup = read.table('group_region_time_all.txt', header = T)

D_group = merge(D, popgroup)
D_group$Group_ID <- gsub("_", " ", D_group$Group_ID)
D_group$Caucase2 <- gsub("_", " ", D_group$Caucase2)

stderr = D_group$Dstat/D_group$Z
D_group$Dmin = D_group$Dstat - stderr
D_group$Dmax = D_group$Dstat + stderr
D_group$p = 2*pnorm(abs(D_group$Z), mean = 0, sd = 1, lower.tail=FALSE) # mulitply per 2 cause Dstat is a two-tailed test
D_group$padj = p.adjust(D_group$p, method = "BH")
D_group$Zadj = ifelse(D_group$padj !=0,
                      qnorm(0.5*D_group$padj, lower.tail = FALSE), # divise per 2 cause Dstat is a two-tailed test
                      abs(D_group$Z))

Dsign = D_group[abs(D_group$Zadj) >3,]

Dsign$pair = paste(Dsign$Group_ID,"_",Dsign$Caucase2)
inverse = Dsign[Dsign$Group_ID == 'Russia Steppe Maikop'| Dsign$Group_ID == 'Russia Steppe Lola',]
inverse$Dstat = -inverse$Dstat
inverse$Dmin = -inverse$Dmin
inverse$Dmax = -inverse$Dmax
inverse$Group_ID = inverse$Caucase2
inverse$Caucase2 = Dsign[Dsign$Group_ID == 'Russia Steppe Maikop'| Dsign$Group_ID == 'Russia Steppe Lola',]$Group_ID
Dsign = rbind(Dsign, inverse)
Dsign = Dsign[Dsign$Group_ID != 'Russia Steppe Maikop'& Dsign$Group_ID != 'Russia Steppe Lola',]

plot_select_BA <- ggplot(Dsign,
                         aes(Group_ID, Dstat, ymin=Dmin, ymax=Dmax, size = SNPs, col = MT, label = MT),)+
  geom_pointrange(col = 'hotpink')+
  coord_flip()+
  theme_bw()+
  scale_size_binned_area(breaks = c(10000), max_size = 2)+
  theme(axis.title.y = element_blank(),
        legend.position = "left")+
  geom_hline(yintercept=0, col = "black",linetype = 3, size =1)+
  facet_grid(Caucase2~., scales = 'free', space = 'free',
             labeller=label_wrap_gen(width=6, multi_line=T))+
  ggtitle('D(Mbuti, MT; Caucase BA1, Caucase BA2)')
plot_select_BA

D_anat = read.table('out_Dstat_mentesh_anat_BA', col.names = c("file", "Outgroup","MT", "caucase",
                                                                  "anat", "Dstat", "Z","ABBA","ABAB","SNPs"))
D_anat[1] = NULL

D_anat$anat <- gsub("_", " ", D_anat$anat)
D_anat$caucase <- gsub("_", " ", D_anat$caucase)

stderr = D_anat$Dstat/D_anat$Z
D_anat$Dmin = D_anat$Dstat - stderr
D_anat$Dmax = D_anat$Dstat + stderr

Dsign_anat = D_anat[abs(D_anat$Z) >3,]
Dsign_anat$pair = paste(Dsign_anat$Group_ID,"_",Dsign_anat$Caucase2)

plot_select_anat_BA <- ggplot(Dsign_anat,
                              aes(caucase, Dstat, ymin=Dmin, ymax=Dmax, size = SNPs, col = MT, label = MT),)+
  geom_pointrange(col = 'hotpink')+
  coord_flip()+
  theme_bw()+
  scale_size_binned_area(breaks = c(10000), max_size = 2)+
  theme(axis.title.y = element_blank(),
        legend.position = 'none')+
  geom_hline(yintercept=0, col = "black",linetype = 3, size =1)+
  facet_grid(anat~., scales = 'free',
             labeller=label_wrap_gen(width=6, multi_line=T))+
  ggtitle('D(Mbuti, MT; Caucase BA, Anatolia BA)')
plot_select_anat_BA


library(gridExtra)

pdf("SuppFig6.pdf", width=14, height=11)
grid.arrange(plot_select_BA, plot_select_anat_BA ,
             layout_matrix = rbind(c(1,2)))

dev.off()
