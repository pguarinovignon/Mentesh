setwd("D:/Cerveau/These/data/samplesMDH/2021/Analyses_Mentesh/ADMIXTURE/")
library(tidyverse)
library(data.table)
library(plotly)
`%notin%` <- Negate(`%in%`)

rawind<-read.table("ID_GroupID_full.txt",header=TRUE)

#let's begin to play, en resamplant les pops avec + de 20 inds
taille_pop = as.data.frame(table(rawind$Group_ID))
small_pops_ind = as.data.frame(rawind[rawind$Group_ID %in% taille_pop[taille_pop$Freq <=20,]$Var1,]$ID)
#boucle pour garder 20 inds au hasard dans les pops qui en contiennent + de 20
sub_inds = list()
j=1
for (i in taille_pop[taille_pop$Freq > 20,]$Var1){
  pop = taille_pop[taille_pop$Var1 == i,]
  inds = rawind[rawind$Group_ID == pop$Var1,]
  sub = as.data.frame(sample(inds$ID,20))
  sub_inds[[j]] <- sub
  j=j+1
}
big_pops_ind = rbindlist(sub_inds)
#la liste finale d'inds
colnames(small_pops_ind) <- colnames(big_pops_ind)
subsample_inds = rbind(small_pops_ind,big_pops_ind)
colnames(subsample_inds) <- c('ind')
#reprendre et ne garder que les col et row correspondant au subsample
rawind_subsample <- rawind[rawind$ID %in% subsample_inds$ind,]



#####faire la même chose pour les pops modernes

d <- read.table('1240v2_HO_MDHm2021_1_NP_eurasie.evec', as.is=TRUE)
colnames(d)=c('ID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','Group_ID')
ind_eurasie <- read.table('eurasie_HO.txt', header=T, as.is=TRUE)
eurasie = merge(ind_eurasie,d)
eurasie_ind = eurasie[c('ID','Group_ID')]
eurasie_taille_pop = as.data.frame(table(eurasie_ind$Group_ID))
eurasie_small_pops_ind = as.data.frame(eurasie_ind[eurasie_ind$Group_ID %in% eurasie_taille_pop[eurasie_taille_pop$Freq <=20,]$Var1,]$ID)
#boucle pour garder 20 inds au hasard dans les pops qui en contiennent + de 20
sub_inds = list()
j=1
for (i in eurasie_taille_pop[eurasie_taille_pop$Freq > 20,]$Var1){
  pop = eurasie_taille_pop[eurasie_taille_pop$Var1 == i,]
  inds = eurasie_ind[eurasie_ind$Group_ID == pop$Var1,]
  sub = as.data.frame(sample(inds$ID,20))
  sub_inds[[j]] <- sub
  j=j+1
}
eurasie_big_pops_ind = rbindlist(sub_inds)
#la liste finale d'inds
colnames(eurasie_small_pops_ind) <- colnames(eurasie_big_pops_ind)
eurasie_subsample_inds = rbind(eurasie_small_pops_ind,eurasie_big_pops_ind)
colnames(eurasie_subsample_inds) <- c('ind')
#reprendre et ne garder que les col et row correspondant au subsample
eurasie_subsample <- eurasie_ind[eurasie_ind$ID %in% eurasie_subsample_inds$ind,]

admixture_subsample = rbind(eurasie_subsample,rawind_subsample)

write.table(admixture_subsample, file = "ind_group_eurasie_HO.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(admixture_subsample$ID, file = "list_ind_eurasie_HO.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

