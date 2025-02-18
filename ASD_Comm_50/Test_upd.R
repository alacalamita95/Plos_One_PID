# library(rhdf5)
# mi = h5read(file = "MI.mat", name = "MI")
# c_mi = h5read(file = "C_mi.mat", name = "C_mi")
# c_tot = h5read(file = "C_tot.mat", name = "C_tot")
# red = h5read(file = "RED.mat", name = "RED")
# syn = h5read(file = "SYN.mat", name = "SYN")

# Loading data
# mi = as.numeric(read.csv(file = "MI.csv", header = F))
# c_mi = read.csv(file = "C_mi.csv", header = F)
# c_tot = read.csv(file = "C_tot.csv", header = F)
# red = as.numeric(t(read.csv(file = "RED.csv", header = F)))
# syn = as.numeric(t(read.csv(file = "SYN.csv", header = F)))
# df = read.csv(file = "DB_29.csv")
# label = c(0, as.numeric(t(read.csv(file = "label_29.csv"))))
# 
# mi_mat = mat.or.vec(nr = 51, nc = 51)
# for(i in 1:length(mi)){
#   mi_mat[c_mi[i, 1], c_mi[i, 2]] = mi[i]
# }
# save(mi_mat, file = "mi_mat_g.RData")
# # h_mi = heatmap(mi_mat, main = "Mutual Information", symm = T, keep.dendro=TRUE)
# library("gplots")
# h_mi = heatmap.2(mi_mat, sepwidth = 0, 
#           trace = "none", density.info = "none",
#           main = "Mutual Information", col = my_palette)
# save(h_mi, file = "h_mi_g.RData")
# 
# 
# ##### RED
# 
# # TARGET
# nsub = 140
# # n = length(unique(c_tot[,1]))
# red_mat = list()
# for(i in 1:nsub){
#   mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
#   tmp = which(c_tot[,1] == i)
#   for(j in tmp){
#     mat_tmp[c_tot[j, 2], c_tot[j, 3]] = red[j]
#   }
#   mat_tmp_t = t(mat_tmp)
#   mat_tmp_t[lower.tri(mat_tmp)] <- mat_tmp[lower.tri(mat_tmp)]
#   # mat_tmp_t = mat_tmp_t[-which(colMeans(mat_tmp_t) == 0),-which(colMeans(mat_tmp_t) == 0)]
#   red_mat[[i]] = mat_tmp_t
# }
# # for(i in 1:length(red_mat)){
# #   print(isSymmetric(red_mat[[i]]))
# # }
# ## A general-purpose adder:
# add <- function(x) Reduce("+", x)
# red_sum = add(red_mat)
# red_sum = red_sum/(length(red_mat)-1)
# # h_red = heatmap(mat_sum, main = "Redundancy Target", symm = T, keep.dendro=TRUE)
# library("gplots")
# h_red = heatmap.2(red_sum, sepwidth = 0, 
#                   trace = "none", density.info = "none",
#                   main = "Redundancy", col = my_palette)
# save(red_sum, file = "red_sum.RData")
# save(h_red, file = "h_red.RData")
# 
# # DRIVER
# red_mat = list()
# for(i in 1:nsub){
#   mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
#   tmp = c(which(c_tot[,2] == i))
#   for(j in tmp){
#     mat_tmp[c_tot[j, 1], c_tot[j, 3]] = red[j]
#   }
#   tmp = c(which(c_tot[,3] == i))
#   for(j in tmp){
#     mat_tmp[c_tot[j, 1], c_tot[j, 2]] = red[j]
#   }
#   # mat_tmp = mat_tmp[-which(colMeans(mat_tmp) == 0),-which(colMeans(mat_tmp) == 0)]
#   red_mat[[i]] = mat_tmp
# }
# # for(i in 1:length(red_mat)){
# #   print(isSymmetric(red_mat[[i]]))
# # }
# ## A general-purpose adder:
# add <- function(x) Reduce("+", x)
# red_sum_d = add(red_mat)
# red_sum_d = red_sum_d/(length(red_mat)-1)
# # h_red = heatmap(mat_sum, main = "Redundancy Target", symm = T, keep.dendro=TRUE)
# library("gplots")
# h_red_d = heatmap.2(red_sum_d, sepwidth = 0,
#                     trace = "none", density.info = "none",
#                     main = "Redundancy Driver")
# save(red_sum_d, file = "red_sum_d_g.RData")
# save(h_red_d, file = "h_red_d_g.RData")
# 
# 
# ############ SYNERGY
# 
# # TARGET
# syn_mat = list()
# for(i in 1:nsub){
#   mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
#   tmp = which(c_tot[,1] == i)
#   for(j in tmp){
#     mat_tmp[c_tot[j, 2], c_tot[j, 3]] = syn[j]
#   }
#   mat_tmp_t = t(mat_tmp)
#   mat_tmp_t[lower.tri(mat_tmp)] <- mat_tmp[lower.tri(mat_tmp)]
#   # mat_tmp_t = mat_tmp_t[-which(colMeans(mat_tmp_t) == 0),-which(colMeans(mat_tmp_t) == 0)]
#   syn_mat[[i]] = mat_tmp_t
# }
# ## A general-purpose adder:
# add <- function(x) Reduce("+", x)
# syn_sum = add(syn_mat)
# syn_sum = syn_sum/(length(syn_mat)-1)
# library("gplots")
# library("RColorBrewer")
# my_palette = colorRampPalette(c("darkred", "white", "darkblue"))(n = 1000)
# h_syn = heatmap.2(syn_sum, sepwidth = 0,
#                   trace = "none", density.info = "none",
#                   main = "Synergy", col = my_palette)
# save(syn_sum, file = "syn_sum_g.RData")
# save(h_syn, file = "h_syn_g.RData")
# 
# # DRIVER
# syn_mat = list()
# for(i in 1:nsub){
#   mat_tmp = mat.or.vec(nr = nsub, nc = nsub)
#   tmp = c(which(c_tot[,2] == i))
#   for(j in tmp){
#     mat_tmp[c_tot[j, 1], c_tot[j, 3]] = syn[j]
#   }
#   tmp = c(which(c_tot[,3] == i))
#   for(j in tmp){
#     mat_tmp[c_tot[j, 1], c_tot[j, 2]] = syn[j]
#   }
#   # mat_tmp = mat_tmp[-which(colMeans(mat_tmp) == 0),-which(colMeans(mat_tmp) == 0)]
#   syn_mat[[i]] = mat_tmp
# }
# ## A general-purpose adder:
# add <- function(x) Reduce("+", x)
# syn_sum_d = add(syn_mat)
# syn_sum_d = syn_sum_d/(length(syn_mat)-1)
# h_syn_d = heatmap.2(syn_sum_d, sepwidth = 0,
#                     trace = "none", density.info = "none",
#                     main = "Synergy Driver")
# save(syn_sum_d, file = "syn_sum_d_g.RData")
# save(h_syn_d, file = "h_syn_d_g.RData")

###### Loads
load("mi_mat.RData")
# load("red_sum.RData")
# load("red_sum_d.RData")
load("syn_sum.RData")
# load("syn_sum_d.RData")
# load("h_mi.RData")
# load("h_red.RData")
# load("h_red_d.RData")
# load("h_syn.RData")
# load("h_syn_d.RData")
# 
# mi_rows = h_mi$rowInd
# row.cluster_mi = as.hclust(h_mi$rowDendrogram)
# sub_grup_mi <- cutree(row.cluster_mi, k = 3)
# red_rows = h_red$rowInd
# row.cluster_red = as.hclust(h_red$rowDendrogram)
# sub_grup_red <- cutree(row.cluster_red, k = 3)
# red_rows_d = h_red_d$rowInd
# row.cluster_red_d = as.hclust(h_red_d$rowDendrogram)
# sub_grup_red_d <- cutree(row.cluster_red_d, k = 3)
# syn_rows = h_syn$rowInd
# row.cluster_syn = as.hclust(h_syn$rowDendrogram)
# sub_grup_syn <- cutree(row.cluster_syn, k = 3)
# syn_rows_d = h_syn_d$rowInd
# row.cluster_syn_d = as.hclust(h_syn_d$rowDendrogram)
# sub_grup_syn_d <- cutree(row.cluster_syn_d, k = 3)
# save(sub_grup_mi, file = "sub_grup_mi.RData")
# save(sub_grup_red, file = "sub_grup_red_g.RData")
# save(sub_grup_red_d, file = "sub_grup_red_d_g.RData")
# save(sub_grup_syn, file = "sub_grup_syn_g.RData")
# save(sub_grup_syn_d, file = "sub_grup_syn_d_g.RData")

#######################

# load("sub_grup_red.RData")
# load("sub_grup_red_d.RData")
# load("sub_grup_syn.RData")
# load("sub_grup_syn_d.RData")

load("DB_50.RData")
load("feat_50.RData")
load("label.RData")
DB = DB[, which(colnames(DB) %in% feat)]
# df = read.csv(file = "DB_29.csv")
# label = c(0, as.numeric(t(read.csv(file = "label_29.csv"))))

###################################
# cor.test(mi_rows, red_rows)
# cor.test(mi_rows, red_rows_d)
# cor.test(mi_rows, syn_rows)
# cor.test(mi_rows, syn_rows_d)
# 
# cor.test(mi_mat, red_sum)
# cor.test(mi_mat, red_sum_d)
# cor.test(mi_mat, syn_sum)
# cor.test(mi_mat, syn_sum_d)
# cor.test(data_thr, red_sum)
# cor.test(data_thr, syn_sum)

############################
library(factoextra)
library(NbClust)
library(cluster)
library(clustertend)
set.seed(314)
fviz_nbclust(syn_sum, FUNcluster = pam, method = "silhouette")
# fviz_nbclust(syn_sum, FUNcluster = kmeans, method = "silhouette")
# clusternum <- NbClust((syn_sum), distance="euclidean", method="kmeans")
pam.syn <- pam(syn_sum, 2, diss = T, metric = "euclidean", stand = T)
# km.syn <- kmeans(syn_sum, 2)
fviz_cluster(pam.syn, data = syn_sum, palette = c("#FC4E07", "#00AFBB", "#E7B800"), ellipse.type = "euclid", 
             star.plot = TRUE, 
             repel = TRUE, 
             ggtheme = theme_minimal(), main = "Synergy")
# label
table(label[which(pam.syn$clustering == 1)])
table(label[which(pam.syn$clustering == 2)])
# table(label, km.syn$cluster)
# table(label[which(pam.syn$clustering == 1)])
# table(label[which(pam.syn$clustering == 2)])
# table(label, pam.syn$clustering)


# library(polycor)
# library(OmicsMarkeR)
# polychor(label, km.syn$cluster)
# polychor(label, km.mi$cluster)
# sorensen(label, km.syn$cluster)
# sorensen(label, km.mi$cluster)
# jaccard(label, km.syn$cluster)
# jaccard(label, km.mi$cluster)

# table(label[which(km.syn$cluster == 3)])
# table(label[which(sub_grup_syn == 4)])
# table(label[which(sub_grup_syn == 5)])
fviz_nbclust(mi_mat, FUNcluster = pam, method = "silhouette")
# fviz_nbclust(mi_mat, FUNcluster = kmeans, method = "silhouette")
pam.mi <- pam(mi_mat, 2, diss = T, metric = "euclidean", stand = T)
# km.mi <- kmeans(mi_mat, 2)
fviz_cluster(pam.mi, data = mi_mat, palette = c("#FC4E07", "#00AFBB", "#E7B800"), ellipse.type = "euclid", 
             star.plot = TRUE, 
             repel = TRUE, 
             ggtheme = theme_minimal(), main = "Mutual Information")
# table(label[which(km.mi$cluster == 1)])
# table(label[which(km.mi$cluster == 2)])
# table(label[which(km.mi$cluster == 3)])
# table(label, km.mi$cluster)
table(label[which(pam.mi$clustering == 1)])
table(label[which(pam.mi$clustering == 2)])
# table(label[which(pam.mi$clustering == 3)])


# library(Hmisc)
# mat_cor = rcorr(as.matrix(t(df)))
# data_thr = mat_cor$r
# diag(data_thr) <- 0
# data_thr[which(mat_cor$P > 0.01)] <- 0
# h_cor = heatmap.2(data_thr, sepwidth = 0,
#                     trace = "none", density.info = "none",
#                     main = "Correlation")
# cor_rows = h_cor$rowInd
# row.cluster_cor = as.hclust(h_cor$rowDendrogram)
# sub_grup_cor <- cutree(row.cluster_cor, k = 3)


# fviz_nbclust(data_thr, FUNcluster = pam, method = "silhouette")
# pam.cor <- pam(data_thr, 2,  metric = "euclidean", stand = FALSE)
# km.cor <- kmeans(data_thr, 2)
# fviz_cluster(pam.cor, data = data_thr, palette = c("#FC4E07", "#00AFBB", "#E7B800"), ellipse.type = "euclid", 
#              star.plot = TRUE, 
#              repel = TRUE, 
#              ggtheme = theme_minimal() )
# table(label[which(km.cor$cluster == 1)])
# table(label[which(km.cor$cluster == 2)])
# table(label[which(sub_grup_cor == 3)])
# table(label[which(sub_grup_cor == 4)])
# table(label[which(sub_grup_cor == 5)])



# jaccard <- function(x,y){
#   # assert_is_character(x)
#   # assert_is_character(y)
#   
#   x <- as.vector(x)
#   y <- as.vector(y)
#   index <- length(intersect(x,y))/length(union(x, y))
#   return(index)
# }
# 
# sorensen <- function(x,y){
#   # assert_is_character(x)
#   # assert_is_character(y)
#   index <- 
#     2*(length(intersect(x,y)))/(2*(length(intersect(x,y)))+
#                                   length(setdiff(x,y))+
#                                   length(setdiff(y,x)))
#   return(index)
# }
# 

library(limma)

boxplot(DB[which(label == 0), 5], DB[which(label == 1), 5])
# Entire dataset
for(i in 1:44){
  # t = t.test(DB[which(label == 0), i], DB[which(label == 1), i])
  t = kruskal.test(list(DB[which(label == 0), i], DB[which(label == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# Synergy cluster 1
pos1 = which(pam.syn$clustering == 1)
lab_syn1 = label[pos1]
DB_syn1 = DB[pos1, ]
for(i in 1:44){
  # t = t.test(DB_syn1[which(lab_syn1 == 0), i], DB_syn1[which(lab_syn1 == 1), i])
  t = kruskal.test(list(DB_syn1[which(lab_syn1 == 0), i], DB_syn1[which(lab_syn1 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# Synergy cluster 2
pos2 = which(pam.syn$clustering == 2)
lab_syn2 = label[pos2]
DB_syn2 = DB[pos2, ]
for(i in 1:44){
  # t = t.test(DB_syn2[which(lab_syn2 == 0), i], DB_syn2[which(lab_syn2 == 1), i])
  t = kruskal.test(list(DB_syn2[which(lab_syn2 == 0), i], DB_syn2[which(lab_syn2 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# MI cluster 1
pos1 = which(pam.mi$clustering == 1)
lab_mi1 = label[pos1]
DB_mi1 = DB[pos1, ]
for(i in 1:44){
  # t = t.test(DB_mi1[which(lab_mi1 == 0), i], DB_mi1[which(lab_mi1 == 1), i])
  t = kruskal.test(list(DB_mi1[which(lab_mi1 == 0), i], DB_mi1[which(lab_mi1 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}
# MI cluster 2
pos2 = which(pam.mi$clustering == 2)
lab_mi2 = label[pos2]
DB_mi2 = DB[pos2, ]
for(i in 1:44){
  # t = t.test(DB_mi1[which(lab_mi1 == 0), i], DB_mi1[which(lab_mi1 == 1), i])
  t = kruskal.test(list(DB_mi2[which(lab_mi2 == 0), i], DB_mi2[which(lab_mi2 == 1), i]))
  if(t$p.value < 0.01)
    print(colnames(DB)[i])
}






# DGE
sample = as.factor(label)
design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(DB), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg <- topTable(fit3, p.value = 0.01, adjust.method = 'fdr', number = Inf)
boxplot(DB$PBX1[which(label == 0)], DB$PBX1[which(label == 1)])
t.test(DB[which(label == 0), 1], DB[which(label == 1), 1])

# Synergy
pos1 = which(km.syn$cluster == 1)
sample = as.factor(label[pos1])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(DB[pos1, ]), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg.syn1 <- topTable(fit3, p.value = 0.05,
                adjust.method = 'fdr', lfc = log2(1.5), number = Inf)
pos2 = which(km.syn$cluster == 2)
sample = as.factor(label[pos2])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(DB[pos2, ]), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg.syn2 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)

# MI
pos1 = which(km.mi$cluster == 1)
sample = as.factor(label[pos1])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(DB[pos1, ]), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg.mi1 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)
pos2 = which(km.mi$cluster == 2)
sample = as.factor(label[pos2])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(DB[pos2, ]), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg.mi2 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)

save(deg, file = "dge_tot.RData")
save(deg.syn1, file = "dge_syn1.RData")
save(deg.mi1, file = "dge_mi1.RData")

library(xlsx)
write.xlsx(deg, file = "dge_tot.xlsx")
write.xlsx(deg.syn1, file = "dge_syn1.xlsx")
write.xlsx(deg.mi1, file = "dge_mi1.xlsx")






#############################
# loading LOI
# loi_gene = as.numeric(t(read.csv(file = "loi_test_gene.csv", header = F)))
# loi_sort = sort(loi_gene)
# for(i in 1:length(loi_sort)){
#   cat("Gene: ", colnames(df)[which(loi_gene == loi_sort[i])], " LOI: ", loi_sort[i], "\n")
# }
# 
# loi_gene_norm = 0
# for(i in 1:length(loi_gene)){
#   
#   loi_gene_norm[i] = (loi_gene[i] - min(loi_gene))/(max(loi_gene) - min(loi_gene))
#   
# }
# mean(loi_gene_norm)
# mean(loi_gene_norm[which(colnames(df) %in% rownames(deg.mi1))])
# loi_sort_norm = sort(loi_gene_norm)
# for(i in 1:length(loi_sort_norm)){
#   cat("Gene: ", colnames(df)[which(loi_gene_norm == loi_sort_norm[i])], " LOI: ", loi_sort_norm[i], "\n")
# }
# 
# 
# 
# # loi_subj = as.numeric(t(read.csv(file = "loi_test_subj.csv", header = F)))
# loi_subj = as.numeric(t(read.csv(file = "loi_test_gaus.csv", header = F)))
# loi_sort_subj = sort(loi_subj)
# 
# loi_subj_norm = 0
# for(i in 1:length(loi_subj)){
#   
#   loi_subj_norm[i] = (loi_subj[i] - min(loi_subj))/(max(loi_subj) - min(loi_subj))
#   
# }
# 
# mean(loi_subj)
# sd(loi_subj)
# range(loi_subj[which(km.syn$cluster == 1)])
# range(loi_subj[which(km.syn$cluster == 2)])
# mean(loi_subj[which(km.syn$cluster == 2)])
# sd(loi_subj[which(km.syn$cluster == 1)])
# # mean(loi_subj[which(km.syn$cluster == 2)])
# 
# 
# range(loi_subj[which(km.mi$cluster == 1)])
# # range(loi_subj[which(km.mi$cluster == 2)])
# # range(loi_subj[which(km.mi$cluster == 3)])
# mean(loi_subj[which(km.mi$cluster == 3)])
# sd(loi_subj[which(km.mi$cluster == 3)])
# # mean(loi_subj[which(km.mi$cluster == 2)])
# # mean(loi_subj[which(km.mi$cluster == 3)])
# 
# mean(loi_subj_norm)
# sd(loi_subj_norm)
# mean(loi_subj_norm[which(km.syn$cluster == 1)])
# sd(loi_subj_norm[which(km.syn$cluster == 1)])
# mean(loi_subj_norm[which(km.mi$cluster == 3)])
# sd(loi_subj_norm[which(km.mi$cluster == 3)])
# 
# kruskal.test(loi_subj_norm, loi_subj_norm[which(km.syn$cluster == 1)])
# kruskal.test(loi_subj_norm[which(km.mi$cluster == 3)], loi_subj_norm[which(km.syn$cluster == 1)])
# tmp = list()
# tmp[[1]] = loi_subj_norm
# tmp[[2]] = loi_subj_norm[which(km.syn$cluster == 1)]
# tmp[[3]] = loi_subj_norm[which(km.mi$cluster == 3)]
# kruskal.test(tmp[c(1,3)])
# tmp2 = list()
# tmp2[[1]] = loi_subj_norm[which(km.syn$cluster == 1)]
# tmp2[[2]] = loi_subj_norm[which(km.mi$cluster == 3)]
# kruskal.test(tmp2)
# 
# 
# boxplot(loi_subj_norm, loi_subj_norm[which(km.syn$cluster == 1)], loi_subj_norm[which(km.mi$cluster == 3)],
#         col = c("darkblue" , "darkgreen", "darkred"), grid = T, names = c("Entire Dataset", "Synergy Cluster", "MI Cluster"),
#         main = "Local O Information"
#         )
# 
# 
# b = boxplot(loi_subj_norm)
# b$out
# tmp = loi_subj_norm[-which(loi_subj_norm %in% b$out)]
# tmp_syn = loi_subj_norm[which(km.syn$cluster == 1)]
# tmp_syn = tmp_syn[-which(tmp_syn %in% b$out)]
# tmp_mi = loi_subj_norm[which(km.mi$cluster == 3)]
# # tmp_mi = tmp_mi[-which(tmp_mi %in% b$out)]
# tmp_syn = tmp_syn/max(tmp)
# tmp_mi = tmp_mi/max(tmp)
# tmp = tmp/max(tmp)
# 
# mean(tmp)
# sd(tmp)
# mean(tmp_syn)
# sd(tmp_syn)
# mean(tmp_mi)
# sd(tmp_mi)
# 
# 
# #Boxplot accuracy
# library(ggplot2)
# library(reshape2)
# library(wesanderson)
# 
# # tmp = cbind(loi_subj_norm, loi_subj_norm[which(km.syn$cluster == 1)], loi_subj_norm[which(km.mi$cluster == 1)])
# tmp = cbind(rep(1, length(loi_subj_norm)), loi_subj_norm)
# tmp = rbind(tmp, cbind(rep(2, length(loi_subj_norm[which(km.syn$cluster == 1)])), loi_subj_norm[which(km.syn$cluster == 1)]))
# tmp = rbind(tmp, cbind(rep(3, length(loi_subj_norm[which(km.mi$cluster == 1)])), loi_subj_norm[which(km.mi$cluster == 1)]))
# colnames(tmp) <- c("Var2", "value")
# tmp = as.data.frame(tmp)
# # box_mat = melt(c(loi_subj_norm, ))
# # levels(box_mat$Var2) <- comm_text
# # box_mat$Var2 = as.vector(sapply(comm_text, function(x) rep(x, 100)))
# # acc_test = mat.or.vec(nr = 10, nc = 2)
# # acc_test[, 1] = seq(1,10)
# # acc_test[, 2] = c(0.5901, 0.6149, 0.6646, 0.7081, 0.5093, 0.6398, 0.8137, 0.5093, 0.7453, 0.7453)
# # colnames(acc_test) <- c("Var2", "value")
# 
# ggplot(tmp, aes(x = Var2, group = Var2, y = value, fill = as.factor(Var2))) +
#   geom_boxplot(notch = F) +
#   scale_fill_manual(values=wes_palette(n=8, name="Zissou1", type = "continuous")) +
#   # theme(axis.text.x = element_text(as.character(pos_imp_gen), size = 16)) +
#   # theme(axis.text.x = element_text(as.character(pos_imp_gen), vjust = 1.025, hjust=1, size = 16)) +
#   theme(axis.title.x = element_blank()) +
#   # scale_x_discrete(limits = c(1,2,3), labels=c("a", "b", "c")) +
#   theme(axis.text.x = element_text(size = 16)) +
#   ylab(label = "Local O Information") +
#   theme(axis.text.y = element_text(size = 17.5)) +
#   theme(axis.title.y = element_text(size = 22, face = "bold")) +
#   theme(legend.position="none")
#   # geom_point() +
#   # geom_point(data = as.data.frame(acc_test), col = 'blue') + 
#   # ylim(0, 1)
# # ggsave(file="ML_Accuracy.eps")



##################################
# load("sub_grup_red_g.RData")
# load("sub_grup_red_d_g.RData")
# load("sub_grup_syn_g.RData")
# load("sub_grup_syn_d_g.RData")
load("syn_sum_g.RData")
load("mi_mat_g.RData")
load("DB_50.RData")
load("feat_50.RData")
load("label.RData")
DB = DB[, which(colnames(DB) %in% feat)]

colnames(mi_mat) <- rownames(mi_mat) <- colnames(DB)
colnames(syn_sum) <- rownames(syn_sum) <- colnames(DB)
ann = read.csv(file = "GPL6883-11606.txt", header = T, sep = "\t")
library(clusterProfiler)
library(factoextra)
library(NbClust)
library(cluster)
library(clustertend)
set.seed(314)
# Entire
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB))]
enrichKEGG(tmp)
enricher(colnames(DB))



# set.seed(123)
# Synergy
fviz_nbclust(syn_sum, FUNcluster = pam, method = "silhouette")
# fviz_nbclust(syn_sum, FUNcluster = kmeans, method = "silhouette")
# clusternum <- NbClust((syn_sum), distance="euclidean", method="kmeans")
pam.syn <- pam(syn_sum, 2, diss = TRUE,  metric = "euclidean", stand = T)
# km.syn <- kmeans(syn_sum, 2)
fviz_cluster(pam.syn, data = syn_sum, palette = c("#FC4E07", "#00AFBB", "#E7B800"), ellipse.type = "euclid", 
             star.plot = TRUE, 
             repel = TRUE, 
             ggtheme = theme_minimal(), main = "Synergy")
colnames(DB)[which(pam.syn$clustering == 1)]
colnames(DB)[which(pam.syn$clustering == 2)]

t1 = colnames(DB)[which(pam.syn$clustering == 1)]
t2 = colnames(DB)[which(pam.syn$clustering == 2)]
save(t1, file = "asd_50_syn1.RData")
save(t2, file = "asd_50_syn2.RData")

# write.csv(colnames(DB)[which(pam.syn$clustering == 1)])

tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.syn$clustering == 1)])]
enrichKEGG(tmp)
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.syn$clustering == 2)])]
t = enrichKEGG(tmp)
t@result$Description
t@result$ID



fviz_nbclust(mi_mat, FUNcluster = pam, method = "silhouette")
# fviz_nbclust(syn_sum, FUNcluster = kmeans, method = "silhouette")
# clusternum <- NbClust((syn_sum), distance="euclidean", method="kmeans")
pam.mi <- pam(mi_mat, 2, diss = T,  metric = "euclidean", stand = T)
# km.mi <- kmeans(mi_mat, 2)
fviz_cluster(pam.mi, data = mi_mat, palette = c("#FC4E07", "#00AFBB", "#E7B800"), ellipse.type = "euclid", 
             star.plot = TRUE, 
             repel = TRUE, 
             ggtheme = theme_minimal(), main = "Mutual Information")
colnames(DB)[which(pam.mi$clustering == 1)]
colnames(DB)[which(pam.mi$clustering == 2)]
t1 = colnames(DB)[which(pam.mi$clustering == 1)]
t2 = colnames(DB)[which(pam.mi$clustering == 2)]
save(t1, file = "asd_50_mi1.RData")
save(t2, file = "asd_50_mi2.RData")



tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.mi$clustering == 1)])]
enrichKEGG(tmp)
tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[which(pam.mi$clustering == 2)])]
enrichKEGG(tmp)




# studying genes communities in synergy clusters
colnames(df)[which(sub_grup_syn == 1)]
colnames(df)[which(sub_grup_syn == 2)]
colnames(df)[which(sub_grup_syn == 3)]
colnames(df)[which(sub_grup_syn == 4)]

library(clusterProfiler)
ge = enrichKEGG(tmp)
tmp[i] = length(which(ge@result$p.adjust < 0.05))
ge_desc[[i]] = ge@result$Description[1:tmp[i]]
ann = read.csv(file = "GPL570-55999.txt", header = T, sep = "\t")




tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_syn == 1)])]
enrichKEGG(tmp)
tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_syn == 2)])]
enrichKEGG(tmp)
tmp2 = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_syn == 3)])]
enrichKEGG(tmp2)
tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_syn == 4)])]
enrichKEGG(tmp)

tmp = unique(ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df))])
intersect(tmp, tmp2)
enrichKEGG(tmp)


df = DB
#################
# GEPHI

### SYN

# Cluster 1
# Nodi
syn_sum_1 = syn_sum[which(pam.syn$clustering == 1), which(pam.syn$clustering == 1)]
nodi_syn1 = mat.or.vec(nr = nrow(syn_sum_1), nc = 2)
colnames(nodi_syn1) <- c("Id", "Label")
nodi_syn1[,1] <- as.numeric(1:ncol(syn_sum_1))
nodi_syn1[,2] <- colnames(df)[which(pam.syn$cluster == 1)]
# nodi_syn[,3] <- as.numeric(pam.syn$clustering)
# Archi
archi_syn1 = mat.or.vec(nr = length(which(lower.tri(syn_sum_1))), nc = 3)
colnames(archi_syn1) <- c("Source", "Target", "Weight")
k=1
for(i in 1:(nrow(syn_sum_1)-1)){
  for(j in (i+1):nrow(syn_sum_1)){
    
    archi_syn1[k,1] <- i
    archi_syn1[k,2] <- j
    archi_syn1[k,3] <- syn_sum_1[i,j]
    k=k+1
  }
}
write.csv(nodi_syn1, file = "GEPHI/nodi_syn1_gene_50.csv", row.names = F)
write.csv(archi_syn1, file = "GEPHI/archi_syn1_gene_50.csv", row.names = F)




# Cluster 2
# Nodi
syn_sum_2 = syn_sum[which(pam.syn$clustering == 2), which(pam.syn$clustering == 2)]
nodi_syn2 = mat.or.vec(nr = nrow(syn_sum_2), nc = 2)
colnames(nodi_syn2) <- c("Id", "Label")
nodi_syn2[,1] <- as.numeric(1:ncol(syn_sum_2))
nodi_syn2[,2] <- colnames(df)[which(pam.syn$cluster == 2)]
# nodi_syn[,3] <- as.numeric(pam.syn$clustering)
# Archi
archi_syn2 = mat.or.vec(nr = length(which(lower.tri(syn_sum_2))), nc = 3)
colnames(archi_syn2) <- c("Source", "Target", "Weight")
k=1
for(i in 1:(nrow(syn_sum_2)-1)){
  for(j in (i+1):nrow(syn_sum_2)){
    
    archi_syn2[k,1] <- i
    archi_syn2[k,2] <- j
    archi_syn2[k,3] <- syn_sum_2[i,j]
    k=k+1
  }
}
write.csv(nodi_syn2, file = "GEPHI/nodi_syn2_gene_50.csv", row.names = F)
write.csv(archi_syn2, file = "GEPHI/archi_syn2_gene_50.csv", row.names = F)






















load("DB_exp_29.RData")
write.csv(DB_exp, file = "DB_29_com.csv", row.names = F, col.names = F)
loi_gene_com = as.numeric(t(read.csv(file = "loi_test_gene_com.csv", header = F)))
loi_sort = sort(loi_gene_com)
for(i in 1:length(loi_sort)){
  cat("Gene: ", colnames(DB_exp)[which(loi_gene_com == loi_sort[i])], " LOI: ", loi_sort[i], "\n")
}
hist(loi_gene_com)

mat_genloi = 0
for(i in 1:ncol(df)){
  # rownames(mat_genloi)[i] = colnames(df)[i]
  if(length(loi_gene_com[which(colnames(DB_exp) == colnames(df)[i])]) == 0)
    mat_genloi[i] = 0
  else mat_genloi[i] = as.numeric(loi_gene_com[which(colnames(DB_exp) == colnames(df)[i])])
  cat("Gene: ", colnames(df)[i], " LOI: ", loi_gene_com[which(colnames(DB_exp) == colnames(df)[i])], "\n")
}
length(which(mat_genloi < 100))



########################
library(Hmisc)
mat_cor = rcorr(as.matrix(df))
data_thr = mat_cor$r
diag(data_thr) <- 0
data_thr[which(mat_cor$P > 0.01)] <- 0
library("gplots")
h_cor = heatmap.2(data_thr, sepwidth = 0,
                  trace = "none", density.info = "none",
                  main = "Correlation")
cor_rows = h_cor$rowInd
row.cluster_cor = as.hclust(h_cor$rowDendrogram)
sub_grup_cor <- cutree(row.cluster_cor, k = 4)

tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_cor == 1)])]
enrichKEGG(tmp)
tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_cor == 2)])]
enrichKEGG(tmp)
tmp2 = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_cor == 3)])]
enrichKEGG(tmp2)
tmp = ann$ENTREZ_GENE_ID[which(ann$Gene.Symbol %in% colnames(df)[which(sub_grup_cor == 4)])]
enrichKEGG(tmp)



##### Bootstrapping
t = list()
for(i in 1:10){
  
  set.seed(1)
  btstr = sample(1:ncol(DB), 9)
  tmp = ann$Entrez_Gene_ID[which(ann$Symbol %in% colnames(DB)[btstr])]
  t[[i]] = enrichKEGG(tmp)
  
}
t[[10]]


