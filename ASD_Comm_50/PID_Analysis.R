# Loading data
mi = as.numeric(read.csv(file = "MI.csv", header = F))
c_mi = read.csv(file = "C_mi.csv", header = F)
c_tot = read.csv(file = "C_tot.csv", header = F)
red = as.numeric(t(read.csv(file = "RED.csv", header = F)))
syn = as.numeric(t(read.csv(file = "SYN.csv", header = F)))
load("DB_50.RData")
load("feat_50.RData")
DB = DB[, which(colnames(DB) %in% feat)]

######## 
# If genes as nodes
nsub = ncol(DB)
# If subjects as nodes
# nsub = nrow(DB)

##### MI
mi_mat = mat.or.vec(nr = nsub, nc = nsub)
for(i in 1:length(mi)){
  mi_mat[c_mi[i, 1], c_mi[i, 2]] = mi[i]
}


##### RED
# TARGET
red_mat = list()
for(i in 1:nsub){
  mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
  tmp = which(c_tot[,1] == i)
  for(j in tmp){
    mat_tmp[c_tot[j, 2], c_tot[j, 3]] = red[j]
  }
  mat_tmp_t = t(mat_tmp)
  mat_tmp_t[lower.tri(mat_tmp)] <- mat_tmp[lower.tri(mat_tmp)]
  # mat_tmp_t = mat_tmp_t[-which(colMeans(mat_tmp_t) == 0),-which(colMeans(mat_tmp_t) == 0)]
  red_mat[[i]] = mat_tmp_t
}
## A general-purpose adder:
add <- function(x) Reduce("+", x)
red_sum = add(red_mat)
red_sum = red_sum/(length(red_mat)-1)


##### SYNERGY

# TARGET
syn_mat = list()
for(i in 1:nsub){
  mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
  tmp = which(c_tot[,1] == i)
  for(j in tmp){
    mat_tmp[c_tot[j, 2], c_tot[j, 3]] = syn[j]
  }
  mat_tmp_t = t(mat_tmp)
  mat_tmp_t[lower.tri(mat_tmp)] <- mat_tmp[lower.tri(mat_tmp)]
  # mat_tmp_t = mat_tmp_t[-which(colMeans(mat_tmp_t) == 0),-which(colMeans(mat_tmp_t) == 0)]
  syn_mat[[i]] = mat_tmp_t
}
## A general-purpose adder:
add <- function(x) Reduce("+", x)
syn_sum = add(syn_mat)
syn_sum = syn_sum/(length(syn_mat)-1)


###### Cluster Analysis
library(factoextra)
library(cluster)
set.seed(123)


colnames(mi_mat) <- rownames(mi_mat) <- DB
colnames(syn_sum) <- rownames(syn_sum) <- DB
ann = read.csv(file = "GPL6883-11606.txt", header = T, sep = "\t")
# If genes as nodes
# Synergy
fviz_nbclust(syn_sum, FUNcluster = pam, method = "silhouette")
pam.syn <- pam(syn_sum, k = 2, diss = T, metric = "euclidean", stand = T)
colnames(DB)[which(pam.syn$clustering == 1)]
colnames(DB)[which(pam.syn$clustering == 2)]

library(clusterProfiler)
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.syn$clustering == 1)])]
enrichKEGG(tmp)
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.syn$clustering == 2)])]
enrichKEGG(tmp)


# MI
fviz_nbclust(mi_mat, FUNcluster = pam, method = "silhouette")
pam.mi <- pam(mi_mat, 2, diss = T,  metric = "euclidean", stand = T)
colnames(DB)[which(pam.mi$clustering == 1)]
colnames(DB)[which(pam.mi$clustering == 2)]

tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.mi$clustering == 1)])]
enrichKEGG(tmp)
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.mi$clustering == 2)])]
enrichKEGG(tmp)




# If subjects as nodes
# Synergy cluster 1
pos1 = which(pam.syn$clustering == 1)
lab_syn1 = label[pos1]
DB_syn1 = DB[pos1, ]
for(i in 1:ncol(DB)){
  t = kruskal.test(list(DB_syn1[which(lab_syn1 == 0), i], DB_syn1[which(lab_syn1 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# Synergy cluster 2
pos2 = which(pam.syn$clustering == 2)
lab_syn2 = label[pos2]
DB_syn2 = DB[pos2, ]
for(i in 1:ncol(DB)){
  t = kruskal.test(list(DB_syn2[which(lab_syn2 == 0), i], DB_syn2[which(lab_syn2 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# MI cluster 1
pos1 = which(pam.mi$clustering == 1)
lab_mi1 = label[pos1]
DB_mi1 = DB[pos1, ]
for(i in 1:ncol(DB)){
  t = kruskal.test(list(DB_mi1[which(lab_mi1 == 0), i], DB_mi1[which(lab_mi1 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# MI cluster 2
pos2 = which(pam.mi$clustering == 2)
lab_mi2 = label[pos2]
DB_mi2 = DB[pos2, ]
for(i in 1:ncol(DB)){
  t = kruskal.test(list(DB_mi2[which(lab_mi2 == 0), i], DB_mi2[which(lab_mi2 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}

