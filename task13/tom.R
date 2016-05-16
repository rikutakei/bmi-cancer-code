## Add color bar matrices to both sides:

source("heatmap-mik.2.R") should be equivalent to

devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")

## Manually reorder heatmap, trees, and color bars for column:

source("heatmap-mik-mod.2.R") should be equivalent to 
devtools::install_github("TomKellyGenetics/heatmap.2x", ref="supr")

# This example takes a subtype "#EF559E" subset and splits the heatmap be mutation status "CDH!_Mut"

bb3_Stat<-bb3[,bb3[8,]=="#EF559E" & cut == 1 & is.na(CDH1_Mt3)==F] ##bb3 is my column colour matrix of clinical and mutation data
dim(dataset)
## tree for genes
tree_exprSL_voom_corr_dist<-as.dendrogram(hclust(as.dist(1-cor(t(dataset)))))
## split columns by mutation status
data_low<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==1 & is.na(CDH1_Mt3)==F & cut == 1])
data_high<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==0 & is.na(CDH1_Mt3)==F & cut == 1])
## trees for each split dataset (correlation distance)
dist_low<-dist(as.dist(1-cor(data_low)))
hc_low<-hclust(dist_low,method='complete')
dist_high<-dist(as.dist(1-cor(data_high)))
hc_high<-hclust(dist_high,method='complete')
## join trees together
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high), height=max(hc_low$height, hc_high$height)+1)
#rr<-RowCols[match(rownames(dataset), rownames(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR)),]
# plot clusters of genes
Cluster<-cutree(as.hclust(tree_exprSL_voom_corr_dist),4)
ColCluster<-c("red", "blue", "green", "orange")[Cluster]
#rr<-cbind(rr, ColCluster)  ## rr is my row colour data for genes
heatmap.mik.mod2(as.matrix(dataset[,match(labels(hc), colnames(dataset))]), scale='none', trace='none', col=bluered(50), ColSideColors=bb3_Stat[,match(labels(hc), colnames(dataset))], Colv=hc, Rowv=tree_exprSL_voom_corr_dist, margin=c(12, 12), dendrogram='both', main = "TCGA Breast Gene Expression Gatza 2011", xlab = "Sample", ylab = "Pathway", cexCol=1.15, cexRow=1.15, keysize=2.25)
## add legend to heatmap
legend("topleft",legend=c(rep("", 11), "Normal", "Tumour", "Metastasis", "", "Ductal", "Lobular", "", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "", "Positive","Negative", "", "Basal","Her2","LumA","LumB","Normal","","Somatic Mutation","Somatic Mutation (Slient)","Wildtype","Somatic Mutation (Gene 1)","Somatic Mutation (Gene 2)", "Somatic Mutation (Gene 3)", "Somatic Mutation (Gene 1 Silent)","Somatic Mutation (Gene 2 Silent)", "Somatic Mutation (Gene 3 Silent)", "Somatic Mutation (Gene 1 and 2)", "", "CDH1 Low",  "CDH1 High"),
       fill=c(rep(NA, 11), "green", "red", "black", "white", "orange", "blue", "white", "green", "yellow", "orange", "red", "white", "magenta","cyan", "white", "#CE2427", "#EF559E", "#423996", "#8FBCE5", "#50A547","white","black","grey50","grey85","red","blue", "green", "palevioletred", "lightblue", "lightgreen", "purple", "white", "green", "red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.9, xjust=0)

tree_exprSL_voom_eu_dist<-as.dendrogram(hclust(dist(dataset)))
data_low<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==1 & is.na(CDH1_Mt3)==F & cut == 1])
data_high<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==0 & is.na(CDH1_Mt3)==F & cut == 1])
dist_low<-dist(t(data_low), method = "euclidean")
hc_low<-hclust(dist_low,method='complete')
dist_high<-dist(t(data_high), method = "euclidean")
hc_high<-hclust(dist_high,method='complete')
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high))
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high), height=max(hc_low$height, hc_high$height)+1)#rr<-RowCols[match(rownames(dataset), rownames(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR)),]
Cluster<-cutree(as.hclust(tree_exprSL_voom_eu_dist),4)
ColCluster<-c("red", "blue", "green", "orange")[Cluster]
#rr[, ncol(rr)]<-ColCluster
heatmap.mik.mod2(as.matrix(dataset[,match(labels(hc), colnames(dataset))]), scale='none', trace='none', col=bluered(50), ColSideColors=bb3_Stat[,match(labels(hc), colnames(dataset))], Colv=hc, Rowv=tree_exprSL_voom_eu_dist, margin=c(12, 12), dendrogram='both', main = "TCGA Breast Gene Expression Gatza 2011", xlab = "Sample", ylab = "Pathway", cexCol=1.15, cexRow=1.15, keysize=2.25)
legend("topleft",legend=c(rep("", 11), "Normal", "Tumour", "Metastasis", "", "Ductal", "Lobular", "", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "", "Positive","Negative", "", "Basal","Her2","LumA","LumB","Normal","","Somatic Mutation","Somatic Mutation (Slient)","Wildtype","Somatic Mutation (Gene 1)","Somatic Mutation (Gene 2)", "Somatic Mutation (Gene 3)", "Somatic Mutation (Gene 1 Silent)","Somatic Mutation (Gene 2 Silent)", "Somatic Mutation (Gene 3 Silent)", "Somatic Mutation (Gene 1 and 2)", "", "CDH1 Low",  "CDH1 High"),
       fill=c(rep(NA, 11), "green", "red", "black", "white", "orange", "blue", "white", "green", "yellow", "orange", "red", "white", "magenta","cyan", "white", "#CE2427", "#EF559E", "#423996", "#8FBCE5", "#50A547","white","black","grey50","grey85","red","blue", "green", "palevioletred", "lightblue", "lightgreen", "purple", "white", "green", "red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.9, xjust=0)
