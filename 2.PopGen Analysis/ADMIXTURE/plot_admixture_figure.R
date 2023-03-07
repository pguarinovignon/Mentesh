library(tidyverse)
library(car)
library(cowplot)
library(Hmisc)
library(reshape2)
library(gridExtra)

cv=read.table("data/sum_CV_HO.txt" , header = T)
ggplot(cv, aes(K,CV))+geom_jitter(aes(col = rep, shape = rep),stroke = 2, alpha = 0.5)+geom_boxplot(alpha = 0)+
  scale_color_manual(values = c('rep01'= '#a6cee3','rep11'= '#a6cee3',
                                'rep02'='#1f78b4','rep12'='#1f78b4',
                                'rep03'='#b2df8a','rep13'='#b2df8a',
                                'rep04'='#33a02c','rep14'='#33a02c',
                                'rep05'='#fb9a99','rep15'='#fb9a99',
                                'rep06'='#e31a1c','rep16'='#e31a1c',
                                'rep07'='#fdbf6f','rep17'='#fdbf6f',
                                'rep08'='#ff7f00','rep18'='#ff7f00',
                                'rep09'='#cab2d6','rep19'='#cab2d6',
                                'rep10'='#6a3d9a','rep20'='#6a3d9a'))+
  scale_shape_manual(values = c('rep01'= 3,'rep11'= 4,
                                'rep02'= 3,'rep12'= 4,
                                'rep03'= 3,'rep13'= 4,
                                'rep04'= 3,'rep14'= 4,
                                'rep05'= 3,'rep15'= 4,
                                'rep06'= 3,'rep16'= 4,
                                'rep07'= 3,'rep17'= 4,
                                'rep08'= 3,'rep18'= 4,
                                'rep09'= 3,'rep19'= 4,
                                'rep10'= 3,'rep20'= 4))+
  theme_minimal()+labs(colour='Replicat', shape = 'Replicat')
  

tbl=read.table("MATRIX_adm_pruned.9.rep05.Q")
fam=read.table("MATRIX_adm_pruned.fam")
ind=read.table("data/ind_group_eurasie_HO.txt", header = T)
pop=read.table("data/group_region_time_figure.txt", header = T)
col9 = c('darkolivegreen3','red','forestgreen','hotpink','gold','darkslategray4','brown','mistyrose2','purple')
names(col9) = levels(mdat$Ancestry)


tbl$Instance_ID = fam$V2
tbl1=merge(tbl,ind, by = "Instance_ID")
anti_join(tbl,tbl1)
tbl2=merge(tbl1,pop, by = "Group_ID", no.dups = T)
anti_join(tbl1,tbl2)

Steppe = tbl2[tbl2$Region == 'Steppe', ]
SCaucase = tbl2[tbl2$Region == 'SCaucasus'& tbl2$Pop != 'IronAge',]
NCaucase = tbl2[tbl2$Region == 'NCaucasus',]
Iran = tbl2[tbl2$Region == 'Iran',]
Anatolia = tbl2[tbl2$Region == 'Anatolia' & tbl2$Group_ID != 'Anatolia_Epipaleolithic',]
Anatolia_EP = tbl2[tbl2$Group_ID == 'Anatolia_Epipaleolithic',]
Levant = tbl2[tbl2$Region == 'Levant' ,]
Mesopotamia = tbl2[tbl2$Region == 'Mesopotamia',]


mdat_Steppe = melt(Steppe, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
               variable.name="Ancestry", value.name="Fraction")
mdat_SCaucase = melt(SCaucase, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
               variable.name="Ancestry", value.name="Fraction")
mdat_NCaucase = melt(NCaucase, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
               variable.name="Ancestry", value.name="Fraction")
mdat_Iran = melt(Iran, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
                   variable.name="Ancestry", value.name="Fraction")
mdat_Anatolia = melt(Anatolia, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
                     variable.name="Ancestry", value.name="Fraction")
mdat_Anatolia_EP = melt(Anatolia_EP, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
                     variable.name="Ancestry", value.name="Fraction")
mdat_Levant = melt(Levant, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
                     variable.name="Ancestry", value.name="Fraction")
mdat_Mesopotamia = melt(Mesopotamia, id.vars=c("Instance_ID", "Group_ID", "Region","Time",'Pop'), 
                   variable.name="Ancestry", value.name="Fraction")


mdat_Steppe$Ancestry = factor(mdat_Steppe$Ancestry, 
                          levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_SCaucase$Ancestry = factor(mdat_SCaucase$Ancestry, 
                          levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_NCaucase$Ancestry = factor(mdat_NCaucase$Ancestry, 
                          levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_Iran$Ancestry = factor(mdat_Iran$Ancestry, 
                              levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_Anatolia$Ancestry = factor(mdat_Anatolia$Ancestry, 
                                levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_Anatolia_EP$Ancestry = factor(mdat_Anatolia_EP$Ancestry, 
                                levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_Levant$Ancestry = factor(mdat_Levant$Ancestry, 
                                levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))
mdat_Mesopotamia$Ancestry = factor(mdat_Mesopotamia$Ancestry, 
                              levels=names(sort(colSums(tbl[, 1:9]), decreasing=TRUE)))


mdat_SCaucase$Pop_f <- factor(mdat_SCaucase$Pop, levels = c('CHG','Neo','Mentesh','LateC','Armenia_C','EBA',
                                                'KuraAraxes','MLBA','IronAge','Med'))
plot_SCaucase <- ggplot(mdat_SCaucase, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_SCaucase

mdat_NCaucase$Pop_f <- factor(mdat_NCaucase$Pop, levels = c('Cauca_EN','KuraAraxes','BA','Maikop'))
plot_NCaucase <- ggplot(mdat_NCaucase, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_NCaucase

mdat_Anatolia$Pop_f <- factor(mdat_Anatolia$Pop, levels = c('Neo','Chalco','BA'))
plot_Anatolia <- ggplot(mdat_Anatolia, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Anatolia

plot_Anatolia_EP <- ggplot(mdat_Anatolia_EP, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Anatolia_EP

mdat_Iran$Pop_f <- factor(mdat_Iran$Pop, levels = c('Meso','Neo','Chalco','BA'))
plot_Iran <- ggplot(mdat_Iran, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Iran

mdat_Levant$Pop_f <- factor(mdat_Levant$Pop, levels = c('Natufian','PPN','Chalco','BA'))
plot_Levant <- ggplot(mdat_Levant, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Levant

mdat_Mesopotamia$Pop_f <- factor(mdat_Mesopotamia$Pop, levels = c('Neo','BA'))
plot_Mesopotamia <- ggplot(mdat_Mesopotamia, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Mesopotamia


mdat_Steppe$Pop_f <- factor(mdat_Steppe$Pop, levels = c('EEHG','WSHG','Neolithic','Yamnaya'))
plot_Steppe <- ggplot(mdat_Steppe, aes(x=Instance_ID, y=Fraction, fill=Ancestry, order=Ancestry)) +
  geom_bar(stat="identity", position="fill", width=1) +
  facet_grid( ~ Pop_f, drop=TRUE, space="free", scales="free", switch = 'both') +
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=15, angle = 90, hjust = 1),
        legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=col9) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))+
  theme(legend.position = "none")
plot_Steppe

#####ALIGN PLOTS TO HAVE SAME HEIGHT
##FIRST CREATE FUNCTION THAT DO WHAT WE WANT
same.size.ggplot <- function(vector.string.graph, # a vector of strings which correspond to Robject ggplot graphs
                             reference.string.graph, # a string of a  Robject ggplot graphs where height and/or height will be taken for reference
                             width = F, # if you wanna adapat only the width
                             height = T # if you wanna adapat only the height
) {
  
  # example: same.size.ggplot(p0rep(c("a", "b"), thre), "a30")
  
  
  which(vector.string.graph %in% reference.string.graph)
  
  newref <- ggplotGrob(get(reference.string.graph))
  ref.width <- newref$widths
  ref.height <- newref$heights
  
  assign(reference.string.graph, newref, env = parent.frame(1))
  
  for(i in seq_along(vector.string.graph)) {
    if(vector.string.graph[i] != reference.string.graph) {
      new <- ggplotGrob(get(vector.string.graph[i]))
      if( width ) {
        new$widths <- ref.width
      }
      if( height ) {
        new$heights <- ref.height
      }
      assign(vector.string.graph[i], new, env = parent.frame(1))
    }
  }
}

#CREATE AN VECTOR OF STRING CALLING THE NAME OF THE GRAPH
myplots.vec=c("plot_Anatolia","plot_Anatolia_EP", "plot_Iran", "plot_Levant", "plot_NCaucase", "plot_Steppe","plot_Mesopotamia")

#use the function
same.size.ggplot(myplots.vec,"plot_SCaucase") #use function define before


plot_grid(plot_SCaucase, plot_NCaucase, plot_Steppe,plot_Iran,plot_Anatolia_EP,plot_Anatolia,plot_Levant,plot_Mesopotamia, nrow = 1)
grid.arrange(plot_SCaucase, plot_NCaucase, plot_Steppe,plot_Iran,plot_Anatolia_EP,plot_Anatolia,plot_Levant,plot_Mesopotamia,
             layout_matrix = rbind(c(01,01,01,01,01,01,01,01,01,01,01,01,01,01,01,01,01,01,02,02,02,02,02,02,02,03,03,03,03,03,04,04,04,04,04,
                                     05,06,06,06,06,07,07,07,07,07,08,08)))
