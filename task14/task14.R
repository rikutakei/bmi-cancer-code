###############################################################################

                        ##R code for task 14

###############################################################################

## metagene direction stuff

###############################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task15/')

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

## First load all of the necessary data into R:

###############################################################################
## Creighton et al data:

files = readLines('./raw/creighton/files.txt')
files = paste('./raw/creighton/', files, sep='')
cr_raw = ReadAffy(filenames = files)
cr_raw = rma(cr_raw) ## RMA normalise the data
cr_raw = exprs(cr_raw) ## change the format into matrix

cr_obsgene = read.csv('./obsgenes/crobsgenes.txt', header=F)
cr_obsgene = as.vector(cr_obsgene[,1])

crclin = read.csv('./clindata/crclin.csv', sep=',', header=T)

###############################################################################
## Fuentes-Mattei et al data:

files = readLines('./raw/fuentes-mattei/files.txt')
files = paste('./raw/fuentes-mattei/', files, sep='')
fm_raw = ReadAffy(filenames = files)
fm_raw = rma(fm_raw) ## RMA normalise the data
fm_raw = exprs(fm_raw) ## change the format into matrix

fm_obsgene = read.csv('./obsgenes/fmobsgenes.txt', header=T)

fmclin = read.csv('./clindata/fmclin.csv', sep=',', header=T)

###############################################################################
## Creighton gene expression analysis stuff

## Make group:
group = crclin$bmiStatus
group = ifelse(group=='obese', 'obese', 'normal/overweight')
design = model.matrix(~group)
fit = lmFit(cr_raw, design)
fit = eBayes(fit)
tt = topTable(fit, coef = 2, adjust = 'BH', n=nrow(cr_raw))
sum(tt$P.Value < 0.05)   # 5278 genes
sum(tt$P.Value < 0.01)   # 1781 genes
sum(tt$adj.P.Val < 0.05) # 9 genes
sum(tt$adj.P.Val < 0.01) # 0 genes
length(which(2^tt$logFC > 1.2)) #only 132 genes

## How many Creighton genes are p < 0.01
rawobsgenes = tt[tt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(rawobsgenes) )) ## 390
## What about in top 799 genes?
rawobsgenes = rownames(rawobsgenes[1:799,])
crolgenes = rawobsgenes[which(rawobsgenes %in% cr_obsgene)]
length(crolgenes) ## 239

## Do the same, but in residual data:

## Fit linear model on the data using all the other clinical variables:
residuals = cr_raw
residuals = t(apply(residuals, 1, function(x) lm(x ~ crclin$age + crclin$race + crclin$menopause + crclin$grade + crclin$LNstatus + crclin$ERstatus + crclin$PRstatus + crclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
resfit = lmFit(residuals, design)
resfit = eBayes(resfit)
restt = topTable(resfit, coef=2, n = nrow(residuals))
sum(restt$P.Value < 0.01)   # 1104 genes
sum(restt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
resobsgenes = restt[restt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(resobsgenes ) )) ## 216
## What about in top 799 genes?
resobsgenes = rownames(resobsgenes[1:799,])
rescrolgenes = resobsgenes[which(resobsgenes %in% cr_obsgene)]
length(rescrolgenes) ## 168

pdf('pdf/venn1.pdf')
## make venn diagram with res, raw, and cr obesity genes:
setlist = list(Residual = resobsgenes, My_genes = rawobsgenes, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)
dev.off()

## Do the same, but in Caucasian-only data:

camat =cr_raw
ind = which(crclin$race == 'w')
camat = camat[,ind] ## 77 caucasian samples in total
caclin = crclin[ind,]

group = caclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'normal/overweight')
design = model.matrix(~group)
cafit = lmFit(camat, design)
cafit = eBayes(cafit)
catt = topTable(cafit, coef = 2, n = nrow(camat))
sum(catt$P.Value < 0.01)   # 2129 genes
sum(catt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
caobsgenes = catt[catt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caobsgenes ) )) ## 271
## What about in top 799 genes?
caobsgenes = rownames(caobsgenes[1:799,])
cacrolgenes = caobsgenes[which(caobsgenes %in% cr_obsgene)]
length(cacrolgenes) ## 148

## Fit linear model on the Caucasian-only data  after removing all the other clinical variables:
caresiduals = camat
caresiduals = t(apply(caresiduals, 1, function(x) lm(x ~ caclin$age + caclin$menopause + caclin$grade + caclin$LNstatus + caclin$ERstatus + caclin$PRstatus + caclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
caresfit = lmFit(caresiduals, design)
caresfit = eBayes(caresfit)
carestt = topTable(caresfit, coef=2, n = nrow(caresiduals))
sum(carestt$P.Value < 0.01)   # 1558 genes
sum(carestt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
caresobsgenes = carestt[carestt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caresobsgenes ) )) ## 166
## What about in top 799 genes?
caresobsgenes = rownames(caresobsgenes[1:799,])
carescrolgenes = caresobsgenes[which(caresobsgenes %in% cr_obsgene)]
length(carescrolgenes) ## 92

pdf('pdf/venn2.pdf')
## make a Venn diagram with ca, cares, and cr obesity genes:
setlist = list(Caucasian_residual = caresobsgenes , Caucasian= caobsgenes, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)
dev.off()

## for each metagenes, identify the genes that are in the ICGC data and use these genes for validataion
cr_symmat = cr_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(cr_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_symmat) = tmpgenes
cr_symmat = cr_symmat[-(which(is.na(rownames(cr_symmat)))),]
cr_symmat = collapseRows(cr_symmat, unique(rownames(cr_symmat)), unique(rownames(cr_symmat)))
cr_symmat = cr_symmat$datETcollapsed
dim(cr_symmat) #13031 genes
icgcgenes = colnames(BLCA)

allobsname = c("rawobsgenes","crolgenes","resobsgenes","rescrolgenes", "caobsgenes","cacrolgenes","caresobsgenes","carescrolgenes")
metalength = matrix(1,8)
rownames(metalength) = allobsname
for (i in 1:length(allobsname)) {
	tmp = get(allobsname[i])
	tmp = mapIds(hgu133a.db, keys = tmp, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	tmp = unique(tmp)
	if (length(which(is.na(tmp))) > 0){
		tmp = tmp[-which(is.na(tmp))]
	}
	tmp = tmp[which(tmp %in% rownames(cr_symmat))]
	tmp = tmp[which(tmp %in% icgcgenes)]
	assign(allobsname[i], tmp)
	metalength[i,1] = length(tmp)
}

metalength
##rawobsgenes     678
##crolgenes       199
##resobsgenes     655
##rescrolgenes    147
##caobsgenes      657
##cacrolgenes     128
##caresobsgenes   651
##carescrolgenes  86

# need to check if the metagenes are going in the same direction

# find the genes that are common across all of the metagenes:
namelist = c("Raw metagene","CR overlap metagene (raw)","Residual metagene","CR overlap metagene (residual)","Raw metagene (Caucasian/raw)","CR overlap metagene (Caucasian/raw)","Residual metagene (Caucasian/residual)","CR overlap metagene (Caucasian/residual)")
checkgenes= c(rawobsgenes,crolgenes,resobsgenes,rescrolgenes, caobsgenes,cacrolgenes,caresobsgenes,carescrolgenes)
checkgenes = table(checkgenes)[which(table(checkgenes) == 8)]
checkgenes = names(checkgenes)

# check for the direction of each of the metagene, using the common genes
# identified above
pdf('pdf/metadirection.pdf')
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	tmpmeta = 1-tmpmeta
	mat = mat[checkgenes,]
	ord = order(tmpmeta)
	heatmap.2(mat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=namelist[i])
}
dev.off()

# the seventh metagene (caucasian+residual) must be flipped so all the metagenes are going in the same direction

## validate these metagenes in creighton data:
metalist = list()
pdf('pdf/degmetacr.pdf')
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (i == 7) {
		tmpmeta = 1-tmpmeta
	}
	metaplot3(mat, tmpmeta, crclin, name = namelist[i])
	metalist[[i]] = tmpmeta

	txt = gsub("genes", "transmat", allobsname[i])
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	assign(txt, trmat)
}
dev.off()

pdf('pdf/cor.pdf')
## make correlation matrix
x = matrix(1,103)
for(i in 1:length(metalist)) {
	x = cbind(x, metalist[[i]])
}
x = x[,-1]
colnames(x) = gsub('genes', '', allobsname)
cormat = cor(x, method = 'pearson')
cormatadj = x[,-which(colnames(x) == "caresobs")]
cormatadj = cor(cormatadj, method = 'pearson')
heatmap.2(cormat, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
heatmap.2(cormatadj, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
dev.off()

## check if this works on ICGC data:
alltransname = gsub('genes', 'transmat', allobsname)
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
bmimetalist = list()

for (i in 1:length(allobsname)) {
	pdfname = gsub('pdf/genes', 'meta.pdf', allobsname[i])
	pdf(pdfname)
	gene = get(allobsname[i])
	transmat = get(alltransname[i])
	tmpmeta = vector()
	for (j in 1:length(cancertypes)) {
		cancer = t(get(cancertypes[j]))
		mat = cancer[gene,]
		mat = standardise_data(mat)
		bmi = get(paste(cancertypes[j], 'bmi', sep=''))
		tmpsvd = t(transmat %*% mat)
		#tmpsvd = rank(tmpsvd[,1])/ncol(mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if (i == 7) {
			tmpsvd = 1-tmpsvd
		}
		main = paste(namelist[i], '(')
		main = paste(main, cancertypes[j], sep = '')
		main = paste(main, ')', sep = '')
		metaplot3(mat, tmpsvd, bmi, name = main)
		tmpmeta = c(tmpmeta, tmpsvd)
	}
	bmimetalist[[i]] = tmpmeta
	dev.off()
}
bmimetalist = as.data.frame(bmimetalist)
colnames(bmimetalist) = gsub('genes', '', allobsname)

###############################################################################
## Gatza pathway stuff

## import pathway gene list from the Gatza paper
files = readLines('./gatzagenelist/pathlist.txt')
files = paste('./gatzagenelist/', files, sep='')
for (i in 1:length(files)) {
	txt = gsub('./gatzagenelist/','',files[i])
	txt = gsub('.txt','',txt)
	genes = read.csv(files[i])
	genes = as.vector(genes[,1])
	assign(txt, genes)
}
paths = gsub('./gatzagenelist/','',files)
paths = gsub('.txt','',paths)

## convert the gene probes into gene symbols:
pathlength = matrix(1,length(paths))
rownames(pathlength) = paths
for (i in 1:length(paths)) {
	pathgenes = get(paths[i])
	pathgenes = mapIds(hgu133a.db, keys = pathgenes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	pathgenes = unique(pathgenes)
	if (length(which(is.na(pathgenes))) > 0){
		pathgenes = pathgenes[-which(is.na(pathgenes))]
	}
	pathgenes = pathgenes[which(pathgenes %in% colnames(BLCA))]
	assign(paths[i], pathgenes)
	pathlength[i,1] = length(pathgenes)
}
pathlength
#			no. of probes	no. common probes
#akt_probes    206				191
#bcat_probes    77  			 75
#e2f1_probes   128  			120
#egfr_probes   419  			399
#er_probes     102  			 97
#her2_probes   212  			199
#ifna_probes    82  			 81
#ifng_probes    88  			 84
#myc_probes    425  			394
#p53_probes    218  			206
#p63_probes    277  			254
#pi3k_probes   220  			209
#pr_probes     212  			202
#ras_probes    300  			281
#src_probes     81  			 77
#stat3_probes  107  			 98
#tgfb_probes   101  			 93
#tnfa_probes    95  			 90

# Check the correlation between the  metagenes created from raw and standardised
pdf('pdf/gatzarawvsstdmeta.pdf')
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	stdmat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
#	stdmat[stdmat < -3] = -3
#	stdmat[stdmat > 3] = 3
	rawsvd = svd(mat)
	rawmeta = rank(rawsvd$v[,1])/ncol(mat)
	stdsvd = svd(stdmat)
	stdmeta = rank(stdsvd$v[,1])/ncol(stdmat)
	plot(stdmeta, rawmeta, main=paths[i])
}
dev.off()




# need to check if the metagenes are going in the same direction

# find the genes that are common across all of the metagenes:
checkgzgenes= c(akt_probes, bcat_probes, e2f1_probes, egfr_probes, er_probes, her2_probes, ifna_probes, ifng_probes, myc_probes, p53_probes, p63_probes, pi3k_probes, pr_probes, ras_probes, src_probes, stat3_probes, tgfb_probes, tnfa_probes)
checkgzgenes = table(checkgzgenes)[which(table(checkgzgenes) == 8)]
checkgzgenes = names(checkgzgenes)

# check for the direction of each of the gatza metagene

# list of genes related to/representing the pathway:
checkgene = c('AKT1', 'CTNNB1', 'E2F1', 'EGFR', 'ESR1', 'ERBB2', 'IFNA1', 'IFNG', 'MYC', 'TP53', 'TP63', 'PIK3CA', 'PGR', 'HRAS', 'SRC', 'STAT3', 'TGFB1', 'TNF')

# make data matrix for the heatmap:
matheat = cr_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# make row colours for heatmap
col = make_col(as.vector(crclin$ERstatus), continuous=F)
col = rbind(ER = col, PR = make_col(as.vector(crclin$PRstatus), continuous=F), HER2 = make_col(as.vector(crclin$HER2status), continuous=F), LN = make_col(as.vector(crclin$LNstatus), continuous=F))

pdf('pdf/gatzametadirection.pdf')
crgatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(cr_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	crgatzametalist[[i]] = tmpmeta
}
dev.off()
crgatzametalist = as.data.frame(crgatzametalist)
colnames(crgatzametalist) = gsub('_probes', '', paths)
crgatzametalist = t(as.matrix(crgatzametalist))

# calculate the correlation of the metagenes with their representing genes:
gatzacorval = c()
for (i in 1:nrow(crgatzametalist)) {
	tmpcor = cor(crgatzametalist[i,], cr_symmat[checkgene[i],])
	gatzacorval = c(gatzacorval, tmpcor)
}
names(gatzacorval) = rownames(crgatzametalist)

# BMI metagenes were all going in the same direction as the gatza metagenes,
# and these metagenes from Gatza paper needs to be flipped:
mgflip = c('akt_probes', 'e2f1_probes', 'egfr_probes', 'er_probes', 'ifna_probes', 'myc_probes', 'p53_probes', 'ras_probes', 'src_probes', 'stat3_probes', 'tnfa_probes')
rawmgflip = c('akt_probes', 'e2f1_probes', 'egfr_probes', 'er_probes', 'ifna_probes', 'myc_probes', 'p53_probes', 'ras_probes', 'src_probes', 'stat3_probes', 'tnfa_probes')

# Check if it works:
pdf('pdf/gatzametadirectioncheck.pdf')
crgatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (paths[i] %in% mgflip) {
		tmpmeta = 1-tmpmeta
	}
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(cr_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	crgatzametalist[[i]] = tmpmeta
}
crgatzametalist = as.data.frame(crgatzametalist)
colnames(crgatzametalist) = gsub('_probes', '', paths)
crgatzametalist = t(as.matrix(crgatzametalist))
heatmap.2x(crgatzametalist, trace='none',scale='none', col='bluered', ColSideColors = col, main='All Gatza Pathway Metagenes', cexRow=1.0)
dev.off()

gatzacor = cor(t(crgatzametalist), method='pearson')
gatzacor2 = cor(t(crgatzametalist), method='spearman')

gatzaord = c('er', 'pr', 'p53', 'bcat', 'e2f1', 'pi3k', 'myc', 'ras', 'ifna', 'ifng', 'akt', 'p63', 'src', 'her2', 'egfr', 'tgfb', 'stat3', 'tnfa')

# check if I get similar clustering as Gatza paper:
pdf('pdf/gatzacheck.pdf')
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2(crgatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(crgatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F)
ord = hclust(dist(crgatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(crgatzametalist))
#ord = rev(ord)
heatmap.2(gatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(gatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
#heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()

###############################################################################
## Use data from Gatza paper, and see if I get the same results as in their
# paper, using svd

files = c('raw/gatza/hgu133a/files.txt', 'raw/gatza/hgu133b/files.txt')

for (i in 1:length(files)) {
	varname = gsub('raw/gatza/', '', files[i])
	varname = gsub('/files.txt', '', varname)
	assign(varname, readLines(files[i]))

	var = get(varname)
	varname2 = gsub('files.txt', '', files[i])
	var = paste(varname2, var, sep='')
	raw = ReadAffy(filenames = var)
	raw = rma(raw)
	raw = exprs(raw)
	assign(paste(varname, '_raw', sep=''), raw)
}

# just use the HGU133A samples for now:
dim(hgu133a_raw) # 22283 by 1060

## import pathway gene list from the Gatza paper
files = readLines('./gatzagenelist/pathlist.txt')
files = paste('./gatzagenelist/', files, sep='')
for (i in 1:length(files)) {
	txt = gsub('./gatzagenelist/','',files[i])
	txt = gsub('.txt','',txt)
	genes = read.csv(files[i])
	genes = as.vector(genes[,1])
	genes = genes[which(genes %in% rownames(hgu133a_raw))]
	assign(txt, genes)
}

# make data matrix for the heatmap:
matheat = hgu133a_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gatzametadirectionoriginal.pdf')
hgu133agatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = hgu133a_raw[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
#	mat[mat < -3] = -3
#	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	ord = order(tmpmeta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	hgu133agatzametalist[[i]] = tmpmeta
}
dev.off()

hgu133agatzametalist = as.data.frame(hgu133agatzametalist)
colnames(hgu133agatzametalist) = gsub('_probes', '', paths)
hgu133agatzametalist = t(as.matrix(hgu133agatzametalist))

origgatzacor = cor(t(hgu133agatzametalist), method='pearson')
origgatzacor2 = cor(t(hgu133agatzametalist), method='spearman')

pdf('pdf/gatzacheckoriginal.pdf', height=7, width=14)
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2(hgu133agatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(hgu133agatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F)
ord = hclust(dist(hgu133agatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(hgu133agatzametalist))
#ord = rev(ord)
dev.off()
pdf('pdf/gatzacororiginal.pdf')
heatmap.2(origgatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(origgatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()




























hgu133a_symmat = hgu133a_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(hgu133a_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(hgu133a_symmat) = tmpgenes
hgu133a_symmat = hgu133a_symmat[-(which(is.na(rownames(hgu133a_symmat)))),]
hgu133a_symmat = collapseRows(hgu133a_symmat, unique(rownames(hgu133a_symmat)), unique(rownames(hgu133a_symmat)))
hgu133a_symmat = hgu133a_symmat$datETcollapsed
dim(hgu133a_symmat) #13031 genes

# make data matrix for the heatmap:
matheat = hgu133a_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gatzametadirectionoriginal.pdf')
hgu133agatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = hgu133a_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	ord = order(tmpmeta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	hgu133agatzametalist[[i]] = tmpmeta
}
dev.off()

hgu133agatzametalist = as.data.frame(hgu133agatzametalist)
colnames(hgu133agatzametalist) = gsub('_probes', '', paths)
hgu133agatzametalist = t(as.matrix(hgu133agatzametalist))

origgatzacor = cor(t(hgu133agatzametalist), method='pearson')
origgatzacor2 = cor(t(hgu133agatzametalist), method='spearman')

pdf('pdf/gatzacheckoriginal.pdf')
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2(hgu133agatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(hgu133agatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F)
ord = hclust(dist(hgu133agatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(hgu133agatzametalist))
#ord = rev(ord)
heatmap.2(origgatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(origgatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()

###############################################################################
## Gatza and BMI metagenes in ICGC

## make metagene with Gatza pathways in ICGC data samples:
gatzametalist = list()
for (i in 1:length(paths)) {
	pdfname = gsub('_probes', 'meta.pdf', paths[i])
	pdfname = paste('pdf/', pdfname, sep='')
	pdf(pdfname)
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	tmpmeta = vector()
	for (j in 1:length(cancertypes)) {
		cancer = t(get(cancertypes[j]))
		mat = cancer[gene,]
		mat = standardise_data(mat)
		bmi = get(paste(cancertypes[j], 'bmi', sep=''))
		tmpsvd = t(trmat %*% mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if (paths[i] %in% mgflip) {
			tmpsvd = 1-tmpsvd
		}
		main = gsub('_probes', '', paths[i])
		main = toupper(main)
		main = paste(main, 'metagene (')
		main = paste(main, cancertypes[j], sep = '')
		main = paste(main, ')', sep = '')
		metaplot3(mat, tmpsvd, bmi, name = main)
		tmpmeta = c(tmpmeta, tmpsvd)
	}
	gatzametalist[[i]] = tmpmeta
	dev.off()
}
gatzametalist = as.data.frame(gatzametalist)
colnames(gatzametalist) = paths

## stick the BMI metagenes and Gatza metagenes together:
allmetagenes = cbind(bmimetalist, gatzametalist)
allmetagenes = t(data.matrix(allmetagenes))

## need to make a combined clinical data for all the cancer types
x = c('bcr_patient_barcode', 'tumor_tissue_site', 'weight', 'height')
allclin = BLCAclin[,x]
for(i in 2:length(cancertypes)) {
	txt = paste(cancertypes[i], 'clin', sep='')
	clin = get(txt)
	clin = clin[,x]
	allclin = rbind(allclin, clin)
}
rownames(allclin) = allclin$bcr_patient_barcode
allclin = allclin[colnames(allmetagenes),-1]
allclin$weight = as.numeric(as.character(allclin$weight))
allclin$height = as.numeric(as.character(allclin$height))
allclin[,4] = allclin$weight/((allclin$height/100)^2)
colnames(allclin)[4] = 'bmi'
x = allclin[,4]
for (i in 1:length(x)) {
	if(x[i] >= 30){
		x[i] = "obese"
	} else if(x[i] < 25){
		x[i] = "normal"
	} else {
		x[i] = "overweight"
	}
}
allclin[,5] = x
colnames(allclin)[5] = 'bmiStatus'

## prepare colours for cancer types:
brewercol = brewer.pal(8, "Set1")
cancertypecol = c()
type = c()
for (i in 1:length(cancertypes)) {
	type = c(type, rep(cancertypes[i], nrow(get(cancertypes[i]))))
	cancertypecol = c(cancertypecol, rep(brewercol[i], nrow(get(cancertypes[i]))))
}
allclin[,1] = type

## prepare colours for BMI status:
brewercol = brewer.pal(3, "Set1")
bmicol = allclin[,4]
for (i in 1:length(bmicol)) {
	if(bmicol[i] >= 30){
		bmicol[i] = brewercol[1]
	} else if(bmicol[i] < 25){
		bmicol[i] = brewercol[2]
	} else {
		bmicol[i] = brewercol[3]
	}
}

#bmicol = ifelse(allclin[,5] == "obese", brewercol[1], brewercol[2])

## make a vector of the order of the samples (by BMI values), per cancer type.
bmiord = vector()
bmivalcol = vector()
count = 0
for (i in 1:length(cancertypes)) {
	type = cancertypes[i]
	cl = allclin[which(allclin$tumor_tissue_site == type),]
	cl = cl$bmi
	ord = order(rank(cl)) + count
	col = bluered(length(ord))[rank(cl)]

	bmiord = c(bmiord, ord)
	bmivalcol = c(bmivalcol, col)
	count = length(bmiord)
}

allcol = rbind(BMI=bmivalcol, BMIStatus=bmicol, CancerType=cancertypecol)

pdf(file='pdf/gatzabmimeta.pdf',width=14, height=7)
heatmap.2x(allmetagenes, scale='none', trace='none', col=bluered(2000), ColSideColors=allcol)
heatmap.2x(allmetagenes[, 1:1872], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol, Colv=F)
heatmap.2x(allmetagenes[, bmiord], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol[,bmiord], Colv=F)
dev.off()

pdf('pdf/allmetacor.pdf')
x = cor(t(allmetagenes), method='pearson')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main="Pearson correlation")
x = cor(t(allmetagenes), method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
x = cor(gatzametalist, method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
dev.off()

# Plot heatmaps of all the metagenes for each cancer type
pdf(file='pdf/cancersepallmeta.pdf', width=7, height=7)
for (i in 1:length(cancertypes)) {
	type = cancertypes[i]
	ind = which(allclin$tumor_tissue_site == type)
	samples = rownames(allclin)[ind]
	mat = allmetagenes[,samples]
	heatmap.2(x=mat, scale='none', trace='none', col='bluered', main=cancertypes[i])
	heatmap.2x(allmetagenes[, bmiord[ind]], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol[c("BMI", "BMIStatus"),bmiord[ind]], Colv=F)
}
dev.off()


