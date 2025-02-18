# Loading data
mi = as.numeric(read.csv(file = "MI.csv", header = F))
c_mi = read.csv(file = "C_mi.csv", header = F)
c_tot = read.csv(file = "C_tot.csv", header = F)
red = as.numeric(t(read.csv(file = "RED.csv", header = F)))
syn = as.numeric(t(read.csv(file = "SYN.csv", header = F)))
df = read.csv(file = "HCC_Comm_29.csv")
label = c(0, as.numeric(t(read.csv(file = "label_29.csv"))))


######## 
# If genes as nodes
nsub = ncol(df)
# If subjects as nodes
# nsub = nrow(df)

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


# If genes as nodes
# Synergy
fviz_nbclust(syn_sum, FUNcluster = pam, method = "silhouette")
pam.syn <- pam(syn_sum, k = 2, diss = T, metric = "euclidean", stand = T)
pam.syn$data = syn_sum
colnames(df)[which(pam.syn$clustering == 1)]
colnames(df)[which(pam.syn$clustering == 2)]
library(xlsx)
# Export genes sub communities in order to analyse on GSEA
# write.xlsx(colnames(df)[which(pam.syn$clustering == 1)], file = "syn_genes_1.xlsx")
# write.xlsx(colnames(df)[which(pam.syn$clustering == 2)], file = "syn_genes_2.xlsx")
# MI
fviz_nbclust(mi_mat, FUNcluster = pam, method = "silhouette")
pam.mi <- pam(mi_mat, k = 2, diss = TRUE, metric = "euclidean", stand = T)
# write.xlsx(colnames(df)[which(pam.mi$clustering == 1)], file = "mi_genes_1.xlsx")
# write.xlsx(colnames(df)[which(pam.mi$clustering == 2)], file = "mi_genes_2.xlsx")





# If subjects as nodes
library(limma)
# DGE
sample = as.factor(label)
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(df), design.mat)
# fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit)
deg <- topTable(fit3, p.value = 0.05,
                   adjust.method = 'fdr', lfc = log2(1.5), number = Inf)

# Synergy
# Cluster 1
pos1 = which(pam.syn$clustering == 1)
sample = as.factor(label[pos1])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(df[pos1, ]), design.mat)
fit3 <- eBayes(fit)
deg.syn1 <- topTable(fit3, p.value = 0.05,
                adjust.method = 'fdr', lfc = log2(1.5), number = Inf)
# Cluster 2
pos2 = which(pam.syn$clustering == 2)
sample = as.factor(label[pos2])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(df[pos2, ]), design.mat)
fit3 <- eBayes(fit)
deg.syn2 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)


# MI
# Cluster 1
pos1 = which(pam.mi$clustering == 1)
sample = as.factor(label[pos1])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(df[pos1, ]), design.mat)
fit3 <- eBayes(fit)
deg.mi1 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)
# Cluster 2
pos2 = which(pam.mi$clustering == 2)
sample = as.factor(label[pos2])
design.mat <- model.matrix(~ 0 + sample)
colnames(design.mat) <- levels(sample)
design.mat
# Fitting Model
fit <- lmFit(t(df[pos2, ]), design.mat)
fit3 <- eBayes(fit)
deg.mi2 <- topTable(fit3, p.value = 0.05,
                     adjust.method = 'fdr', lfc = log2(1.5), number = Inf)

