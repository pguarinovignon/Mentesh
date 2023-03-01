ind <- read.table("list_ind_eurasie_HO.txt")
colnames(ind) = c("ID")

pedind=read.table("1240v2_HO_NP_final_MDH_laza2022.pedind")

ped=data.frame(pedind$V1, pedind$V2)
colnames(ped)=c("nbID","ID")
mergeind = merge(ind, ped)

admixture_ind = data.frame(mergeind$nbID,mergeind$ID)
write.table(admixture_ind,file="tokeep_admixture.txt", sep = " ", dec = ".",row.names = FALSE, col.names = FALSE, quote = FALSE)
