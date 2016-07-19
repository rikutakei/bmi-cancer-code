###############################################################################

                        ##R code for task 14

###############################################################################

## metagene direction stuff

###############################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task14/')

#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task13/functions.R')

# install packages
source('http://bioconductor.org/biocLite.R')
install.packages("data.table")
install.packages("gplots")
install.packages("RColorBrewer")
biocLite()
biocLite("WGCNA")
biocLite("hgu133a.db")
biocLite("hgu133b.db")
biocLite("hgu133plus2.db")
biocLite("tidyr")
biocLite("dplyr")
biocLite("limma")
biocLite("edgeR")
biocLite("DESeq")
biocLite("affy")
biocLite("org.Hs.eg.db")
biocLite("KEGG.db")
biocLite("GO.db")
biocLite("reactome.db")
biocLite("GMD")
biocLite("sva")

# load libraries:

library(data.table)
library(gplots)
library(lattice)
library(GMD)
library(RColorBrewer)
library(mclust)
library(colorRamps)

#heatmap stuff
library("devtools")
devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")
library("heatmap.2x")

#libraries for data wrangling
library(WGCNA)
library(hgu133a.db)
library(hgu133b.db)
library(hgu133plus2.db)
library(tidyr)
library(dplyr)
library(sva)

#libraries for DEG analysis
library(limma)
library(edgeR)
library(DESeq)
library(affy)

#libraries for pathway analysis
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(reactome.db)

##for venn diagram
source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R')

###############################################################################
## Use data from Gatza paper, and see if I get the same results as in their
# paper, using svd

files = readLines('./raw/gatza/hgu133a/files.txt')
files = paste('./raw/gatza/hgu133a/', files, sep='')
gatzaraw = ReadAffy(filenames = files)
gatzarawrma = rma(gatzaraw) ## RMA normalise the data
gatzarawmas5 = mas5(gatzaraw) ## MAS5 normalise the data
gatzarawrma = exprs(gatzarawrma) ## change the format into matrix
gatzarawmas5 = exprs(gatzarawmas5) ## change the format into matrix
gatzarawmas5 = log2(gatzarawmas5) ## need to log2 the data, as in their paper
dim(gatzarawrma) # 22283 by 1060
dim(gatzarawmas5) # 22283 by 1060

# Import the sample names included in different GSE experiments:
files = readLines('./raw/gatza/sample_names/files.txt')
files = paste('./raw/gatza/sample_names/', files, sep='')

sample_names = c()
batch = c()
for (i in 1:length(files)) {
	tmp = readLines(files[i])
	txt = gsub('.*/', '', files[i])
	txt = gsub('_.*', '', txt)
	sample_names = c(sample_names, tmp)
	tmpbatch = rep(txt, length(tmp))
	batch = c(batch, tmpbatch)
}
batch = cbind(sample_names, batch)
colnames(batch) = c("samples", "batch")
rownames(batch) = gsub('.gz', '', batch[,1])

# Only interested in the samples that are from HGU133A microarray chip:
batch = batch[which(rownames(batch) %in% colnames(gatzarawmas5)),]

# Correct for the batch effect:
# First make a model matrix for the different batches:
mod = model.matrix(~1, data=as.data.frame(batch))

# Adjust for batch effect:
gatzabatchrma = ComBat(dat = gatzarawrma, batch=batch[,2], mod=mod)
gatzabatchmas5 = ComBat(dat = gatzarawmas5, batch=batch[,2], mod=mod)

# make a matrix of gene symbols for checking the direction of specific genes
gatzasymmas5 = gatzabatchmas5
tmpgenes = mapIds(hgu133a.db, keys = rownames(gatzasymmas5), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(gatzasymmas5) = tmpgenes
gatzasymmas5 = gatzasymmas5[-(which(is.na(rownames(gatzasymmas5)))),]
gatzasymmas5 = collapseRows(gatzasymmas5, unique(rownames(gatzasymmas5)), unique(rownames(gatzasymmas5)))
gatzasymmas5 = gatzasymmas5$datETcollapsed
dim(gatzasymmas5) #13031 genes

gatzasymrma = gatzabatchrma
tmpgenes = mapIds(hgu133a.db, keys = rownames(gatzasymrma), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(gatzasymrma) = tmpgenes
gatzasymrma = gatzasymrma[-(which(is.na(rownames(gatzasymrma)))),]
gatzasymrma = collapseRows(gatzasymrma, unique(rownames(gatzasymrma)), unique(rownames(gatzasymrma)))
gatzasymrma = gatzasymrma$datETcollapsed
dim(gatzasymrma) #13031 genes

## import pathway gene list from the Gatza paper
files = readLines('./gatzagenelist/pathlist.txt')
files = paste('./gatzagenelist/', files, sep='')
for (i in 1:length(files)) {
	txt = gsub('./gatzagenelist/','',files[i])
	txt = gsub('.txt','',txt)
	genes = read.csv(files[i])
	genes = as.vector(genes[,1])
	genes = genes[which(genes %in% rownames(gatzaraw))]
	assign(txt, genes)
}
head(akt_probes)

# make a variable with list of the name of the pathways
paths = gsub('./gatzagenelist/', '', files)
paths = gsub('.txt', '', paths)
paths

# list of genes related to/representing the pathway:
checkgene = c('AKT1', 'CTNNB1', 'E2F1', 'EGFR', 'ESR1', 'ERBB2', 'IFNA1', 'IFNG', 'MYC', 'TP53', 'TP63', 'PIK3CA', 'PGR', 'HRAS', 'SRC', 'STAT3', 'TGFB1', 'TNF')

# make data matrix for the heatmap:
# This matrix is only used to visualise the direction of the metagene, and is not used in the creation of metagene or the transformation matrix
#matheat = gatzasymrma[checkgene,]
matheat = gatzasymmas5[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# Show the direction of metagenes:
#pdf('pdf/gatzarmametadir.pdf')
pdf('pdf/gatzamas5metadir.pdf')
mgfliprma = c(
			   #'akt_probes',
			   'bcat_probes',
			   'e2f1_probes',
			   'egfr_probes',
			   'er_probes',
			   'her2_probes',
			   'ifna_probes',
			   'ifng_probes',
			   'myc_probes',
			   'p53_probes',
			   #'p63_probes',
			   'pi3k_probes',
			   'pr_probes',
			   'ras_probes',
			   #'src_probes',
			   #'stat3_probes',
			   'tgfb_probes'#,
			   #'tnfa_probes'
			   )
mgflipmas5 = c(
			   #'akt_probes',
			   'bcat_probes',
			   'e2f1_probes',
			   'egfr_probes',
			   'er_probes',
			   'her2_probes',
			   'ifna_probes',
			   'ifng_probes',
			   'myc_probes',
			   'p53_probes',
			   #'p63_probes',
			   'pi3k_probes',
			   'pr_probes',
			   'ras_probes',
			   #'src_probes',
			   #'stat3_probes',
			   'tgfb_probes'#,
			   #'tnfa_probes'
			   )

gatzametalist = list()
gatzacor = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	#mat = gatzabatchrma[gene,]
	mat = gatzabatchmas5[gene,]
	#mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	tmpsvd = svd(mat)
	tmpmeta = tmpsvd$v[,1]
	#if (paths[i] %in% mgfliprma) {
	if (paths[i] %in% mgflipmas5) {
		tmpmeta = 1-tmpmeta
	}
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	ord = order(tmpmeta)
	col = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(matheat[checkgene[i],])]
	col = rbind(col2, Metagene=col)
	rownames(col)[1] = checkgene[i]
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = col[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	gatzametalist[[i]] = tmpmeta
	cor = cor(tmpmeta, matheat[i,])
	gatzacor[[i]] = cor
}
gatzametalist = as.data.frame(gatzametalist)
colnames(gatzametalist) = gsub('_probes', '', paths)
gatzametalist = t(as.matrix(gatzametalist))

names(gatzacor) = paths
gatzacor = as.data.frame(gatzacor)

tmpord = c('er', 'pr', 'p53', 'bcat', 'e2f1', 'pi3k', 'myc', 'ras', 'ifna', 'ifng', 'akt', 'p63', 'src', 'her2', 'egfr', 'tgfb', 'stat3', 'tnfa')

heatmap.2(gatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1, main='Gatza metagenes')
heatmap.2(gatzametalist[tmpord,], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F, main='Gatza metagenes ordered as in their paper')

cor = cor(t(gatzametalist), method='pearson')
heatmap.2(cor, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(cor[tmpord, tmpord], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F, Colv=F)
dev.off()













origgatzacor = cor(t(gatzametalist), method='pearson')
origgatzacor2 = cor(t(gatzametalist), method='spearman')

pdf('pdf/gatzacheckoriginal.pdf', height=7, width=14)
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2(gatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(gatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F)
ord = hclust(dist(gatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(gatzametalist))
#ord = rev(ord)
dev.off()
pdf('pdf/gatzacororiginal.pdf')
heatmap.2(origgatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(origgatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()
