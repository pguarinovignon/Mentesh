args<-commandArgs(TRUE)
comb = read.table(args[1])
popleft = read.table(args[2])
popright = read.table(args[3])
target = args[4]
N = args[5]

if (N==2) {
for (i in 1:nrow(comb)){
  l=t(comb[i,])
  l=rbind(target,l)
  r=as.data.frame(popleft[popleft$V1 != as.character(comb[i,1]) & popleft$V1 != as.character(comb[i,2]),])
  colnames(r) = c("V1")
  r=rbind(popright,r)
  write.table(l, file = paste0("left_",target,"_set1_2pop_", i, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(r, file = paste0("right_",target,"_set1_2pop_", i, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}
}
else {
for (i in 1:nrow(comb)){
  l=t(comb[i,])
  l=rbind(target,l)
  r=as.data.frame(popleft[popleft$V1 != as.character(comb[i,1]) & popleft$V1 != as.character(comb[i,2]) & popleft$V1 != as.character(comb[i,3]),])
  colnames(r) = c("V1")
  r=rbind(popright,r)
  write.table(l, file = paste0("left_",target,"_set1_3pop_", i, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(r, file = paste0("right_",target,"_set1_3pop_", i, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}
}
q()
