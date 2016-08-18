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

# Import the sample names included in different GSE experiments:
files = readLines('./raw/gatza/sample_names/files.txt')
files = paste('./raw/gatza/sample_names/', files, sep='')

sample_list = list()
for (i in 1:length(files)) {
	tmp = readLines(files[i])
	txt = gsub('.*/', '', files[i])
	txt = gsub('_.*', '', txt)
	sample_list[[i]] = tmp
	names(sample_list)[i] = txt
}
sample_list = lapply(sample_list, function(x) gsub('.gz','', x))

# Import the hgu133a samples:
files = readLines('./raw/gatza/hgu133a/files.txt')
sample_list = lapply(sample_list, function(x) x[which(x %in% files)])

# remove GSE dataset with no HGU133A dataset:
sample_list = sample_list[lapply(sample_list, length) > 0]
sample_list = lapply(sample_list, function(x) paste('./raw/gatza/hgu133a/', x, sep=''))

# import each dataset separately
rawlist = sample_list
rawlist = lapply(rawlist, function(x) ReadAffy(filenames = x))

# RMA normalise:
rmalist = lapply(rawlist, function(x) rma(x))
rmalist = lapply(rmalist, function(x) exprs(x))

# MAS5 normalise:
mas5list = lapply(rawlist, function(x) mas5(x))
mas5list = lapply(mas5list, function(x) exprs(x))
mas5list = lapply(mas5list, function(x) log2(x))

# combine the independently normalised data into a single matrix:
gatzarma = do.call(cbind, rmalist)
gatzamas5 = do.call(cbind, mas5list)

dim(gatzarma) # 22283 by 1060
dim(gatzamas5) # 22283 by 1060

# generate batch info:
batch = sample_list
batch = lapply(batch, function(x) gsub('./raw/gatza/hgu133a/', '' , x))
batch = as.vector(unlist(batch))
groups = c()
for (i in 1:length(sample_list)) {
	tmp = rep(names(sample_list)[[i]], length(sample_list[[i]]))
	groups = c(groups, tmp)
}
batch = cbind(batch, groups)

# Correct for the batch effect:
# First make a model matrix for the different batches:
mod = model.matrix(~1, data=as.data.frame(batch))

# Adjust for batch effect:
gatzabatchrma = ComBat(dat = gatzarma, batch=batch[,2], mod=mod)
gatzabatchmas5 = ComBat(dat = gatzamas5, batch=batch[,2], mod=mod)

# standardise the data:
gtrmastd = t(apply(gatzabatchrma, 1, function(x) (x-mean(x))/sd(x)))
gtmasstd = t(apply(gatzabatchmas5, 1, function(x) (x-mean(x))/sd(x)))

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
	genes = readLines(files[i])
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
# template = c(
# 			   'akt_probes',
# 			   'bcat_probes',
# 			   'e2f1_probes',
# 			   'egfr_probes',
# 			   'er_probes',
# 			   'her2_probes',
# 			   'ifna_probes',
# 			   'ifng_probes',
# 			   'myc_probes',
# 			   'p53_probes',
# 			   'p63_probes',
# 			   'pi3k_probes',
# 			   'pr_probes',
# 			   'ras_probes',
# 			   'src_probes',
# 			   'stat3_probes',
# 			   'tgfb_probes',
# 			   'tnfa_probes'
# 			   )


# Show the direction of metagenes:
# make data matrix for the heatmap:
# This matrix is only used to visualise the direction of the metagene, and is not used in the creation of metagene or the transformation matrix
matheat = gatzasymrma
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat = gatzasymrma[checkgene,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# RMA/Standardised/rank-based
pdf('pdf/gtrmastdrank.pdf')
mgfliprma = c(
			   # 'akt_probes',
			   'bcat_probes',
			   'e2f1_probes',
			   'egfr_probes',
			   # 'er_probes',
			   'her2_probes',
			   'ifna_probes',
			   'ifng_probes',
			   # 'myc_probes',
			   'p53_probes',
			   'p63_probes',
			   'pi3k_probes',
			   'pr_probes',
			   # 'ras_probes',
			   # 'src_probes',
			   'stat3_probes',
			   'tgfb_probes',
			   'tnfa_probes'
			  )
gtrmametacor1 = gatzaPath(mat = gtrmastd, matheat = matheat, pathlist = paths, flip = mgfliprma, metalist = "gatzametarma", corlist = "gatzacorrma",  rank = T, checkgene = checkgene)
dev.off()

# RMA/Standardised/probit
pdf('pdf/gtrmastdprobit.pdf')
mgfliprma = c( 'bcat_probes', 'er_probes', 'ifna_probes', 'ifng_probes', 'myc_probes', 'pi3k_probes', 'pr_probes', 'ras_probes', 'src_probes', 'tnfa_probes')
gtrmametacor2 = gatzaPath(mat = gtrmastd, matheat = matheat, pathlist = paths, flip = mgfliprma, metalist = "gatzametarma", corlist = "gatzacorrma",  rank = F, checkgene = checkgene)
dev.off()

# RMA/non-standardised/rank-based
pdf('pdf/gtrmarank.pdf')
mgfliprma = c( 'bcat_probes', 'e2f1_probes', 'er_probes', 'ifna_probes', 'ifng_probes', 'myc_probes', 'p53_probes', 'pi3k_probes', 'pr_probes', 'ras_probes')
gtrmametacor3 = gatzaPath(mat = gatzabatchrma, matheat = matheat, pathlist = paths, flip = mgfliprma, metalist = "gatzametarma", corlist = "gatzacorrma",  rank = T, checkgene = checkgene)
dev.off()

# RMA/non-standardised/probit
pdf('pdf/gtrmaprobit.pdf')
mgfliprma = c( 'bcat_probes', 'e2f1_probes', 'er_probes', 'ifna_probes', 'ifng_probes', 'myc_probes', 'p53_probes', 'pi3k_probes', 'pr_probes', 'ras_probes', 'tgfb_probes')
gtrmametacor4 = gatzaPath(mat = gatzabatchrma, matheat = matheat, pathlist = paths, flip = mgfliprma, metalist = "gatzametarma", corlist = "gatzacorrma",  rank = F, checkgene = checkgene)
dev.off()

pdf(file='pdf/rmametaplot.pdf', width=7, height=7)
	x1 = gtrmametacor1$metagene
	x2 = gtrmametacor2$metagene
	x3 = gtrmametacor3$metagene
	x4 = gtrmametacor4$metagene
	for (i in 1:nrow(x1)) {
		main = rownames(x1)[i]
		plot(x1[i,], x2[i,], main=main, xlab="Std/rank", ylab="Std/probit")
		plot(x3[i,], x4[i,], main=main, xlab="Non-std/rank", ylab="Non-std/probit")
		plot(x1[i,], x3[i,], main=main, xlab="Std/rank", ylab="Non-std/rank")
		plot(x2[i,], x4[i,], main=main, xlab="Std/probit", ylab="Non-std/probit")
	}
dev.off()


# make data matrix for the heatmap:
# This matrix is only used to visualise the direction of the metagene, and is not used in the creation of metagene or the transformation matrix
matheat = gatzasymmas5
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat = matheat[checkgene,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# MAS5/Standardised/rank-based
pdf('pdf/gtmasstdrank.pdf')
mgflipmas5 = c(
 			   'akt_probes',
 			   'bcat_probes',
 			   # 'e2f1_probes',
			   # 'egfr_probes',
 			   # 'er_probes',
 			   # 'her2_probes',
			   'ifna_probes',
 			   'ifng_probes',
 			   # 'myc_probes',
 			   'p53_probes',
 			   'p63_probes',
 			   'pi3k_probes',
 			   'pr_probes',
 			   # 'ras_probes',
 			   'src_probes',
			   'stat3_probes',
 			   'tgfb_probes',
			   'tnfa_probes'
 			   )
gtmas5metacor1 = gatzaPath(mat = gtmasstd, matheat = matheat, pathlist = paths, flip =mgflipmas5, metalist = "gatzametamas5", corlist = "gatzacormas5", rank = T, checkgene = checkgene)
dev.off()

# MAS5/Standardised/probit
pdf('pdf/gtmasstdprobit.pdf')
mgflipmas5 = c( 'bcat_probes', 'egfr_probes', 'er_probes', 'ifna_probes', 'ifng_probes', 'myc_probes', 'p63_probes', 'pr_probes', 'ras_probes', 'stat3_probes', 'tgfb_probes', 'tnfa_probes')
mgflipmas5 = c(
 			   'akt_probes',
 			   'bcat_probes',
 			   # 'e2f1_probes',
			   # 'egfr_probes',
 			   # 'er_probes',
 			   # 'her2_probes',
			   'ifna_probes',
 			   'ifng_probes',
 			   # 'myc_probes',
 			   'p53_probes',
 			   'p63_probes',
 			   'pi3k_probes',
 			   'pr_probes',
 			   # 'ras_probes',
 			   'src_probes',
			   'stat3_probes',
 			   'tgfb_probes',
			   'tnfa_probes'
 			   )
gtmas5metacor2 = gatzaPath(mat = gtmasstd, matheat = matheat, pathlist = paths, flip =mgflipmas5, metalist = "gatzametamas5", corlist = "gatzacormas5",  rank = F, checkgene = checkgene)
dev.off()

# MAS5/non-standardised/rank-based
pdf('pdf/gtmasrank.pdf')
mgflipmas5 = c( 'bcat_probes', 'e2f1_probes', 'er_probes', 'ifna_probes', 'ifng_probes', 'myc_probes', 'p53_probes', 'pi3k_probes', 'pr_probes', 'ras_probes')
mgflipmas5 = c(
 			   # 'akt_probes',
 			   'bcat_probes',
 			   'e2f1_probes',
			   # 'egfr_probes',
 			   'er_probes',
 			   # 'her2_probes',
			   # 'ifna_probes',
 			   # 'ifng_probes',
 			   'myc_probes',
 			   'p53_probes',
 			   # 'p63_probes',
 			   'pi3k_probes',
 			   'pr_probes',
 			   'ras_probes'#,
 			   'src_probes',
			   # 'stat3_probes',
 			   # 'tgfb_probes',
			   # 'tnfa_probes'
 			   )
gtmas5metacor3 = gatzaPath(mat = gatzabatchmas5, matheat = matheat, pathlist = paths, flip =mgflipmas5, metalist = "gatzametamas5", corlist = "gatzacormas5",  rank = T, checkgene = checkgene)
dev.off()

# MAS5/non-standardised/probit
pdf('pdf/gtmasprobit.pdf')
mgflipmas5 = c( 'bcat_probes', 'e2f1_probes', 'er_probes', 'myc_probes', 'p53_probes', 'pi3k_probes', 'pr_probes', 'ras_probes')
mgflipmas5 = c(
 			   # 'akt_probes',
 			   'bcat_probes',
 			   'e2f1_probes',
			   # 'egfr_probes',
			   'er_probes',
 			   # 'her2_probes',
			   # 'ifna_probes',
 			   # 'ifng_probes',
 			   'myc_probes',
 			   'p53_probes',
 			   # 'p63_probes',
 			   'pi3k_probes',
 			   'pr_probes',
 			   'ras_probes'#,
 			   # 'src_probes',
			   # 'stat3_probes',
 			   # 'tgfb_probes',
			   # 'tnfa_probes'
 			   )
gtmas5metacor4 = gatzaPath(mat = gatzabatchmas5, matheat = matheat, pathlist = paths, flip =mgflipmas5, metalist = "gatzametamas5", corlist = "gatzacormas5",  rank = F, checkgene = checkgene)
dev.off()

pdf(file='pdf/mas5metaplot.pdf', width=7, height=7)
	x1 = gtmas5metacor1$metagene
	x2 = gtmas5metacor2$metagene
	x3 = gtmas5metacor3$metagene
	x4 = gtmas5metacor4$metagene
	for (i in 1:nrow(x1)) {
		main = rownames(x1)[i]
		plot(x1[i,], x2[i,], main=main, xlab="Std/rank", ylab="Std/probit")
		plot(x3[i,], x4[i,], main=main, xlab="Non-std/rank", ylab="Non-std/probit")
		plot(x1[i,], x3[i,], main=main, xlab="Std/rank", ylab="Non-std/rank")
		plot(x2[i,], x4[i,], main=main, xlab="Std/probit", ylab="Non-std/probit")
	}
dev.off()

pdf(file='pdf/combmetaplot.pdf', width=7, height=7)
	x1 = gtrmametacor1$metagene
	x2 = gtrmametacor2$metagene
	x3 = gtrmametacor3$metagene
	x4 = gtrmametacor4$metagene
	y1 = gtmas5metacor1$metagene
	y2 = gtmas5metacor2$metagene
	y3 = gtmas5metacor3$metagene
	y4 = gtmas5metacor4$metagene
	for (i in 1:nrow(x1)) {
		txt = rownames(x1)[i]
		main = paste(txt, ' (Std/rank)', sep='')
		plot(x1[i,], y1[i,], main=main, xlab="rma", ylab="mas5")
		main = paste(txt, ' (Std/probit)', sep='')
		plot(x2[i,], y2[i,], main=main, xlab="rma", ylab="mas5")
		main = paste(txt, ' (Non-std/rank)', sep='')
		plot(x3[i,], y3[i,], main=main, xlab="rma", ylab="mas5")
		main = paste(txt, ' (Non-std/probit)', sep='')
		plot(x4[i,], y4[i,], main=main, xlab="rma/rank/NS", ylab="mas5")
	}
dev.off()

###############################################################################
# Make transformation matrix and apply it to other datasets

# import all the data:

## Creighton et al data:

files = readLines('./raw/creighton/files.txt')
files = paste('./raw/creighton/', files, sep='')
crraw = ReadAffy(filenames = files)

crrma = rma(crraw) ## RMA normalise the data
crrma = exprs(crrma) ## change the format into matrix

crmas = mas5(crraw) ## Mas5 normalise the data
crmas = exprs(crmas) ## change the format into matrix
crmas = log2(crmas) ## log2 the data

crstdrma = t(apply(crrma, 1, function(x) (x-mean(x))/sd(x)))
crstdmas = t(apply(crmas, 1, function(x) (x-mean(x))/sd(x)))

crobsgenes = read.csv('./obsgenes/crobsgenes.txt', header=F)
crobsgenes = as.vector(crobsgenes[,1])

crclin = read.csv('./clindata/crclin.csv', sep=',', header=T)

crsymrma = crrma
tmpgenes = mapIds(hgu133a.db, keys = rownames(crsymrma), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(crsymrma) = tmpgenes
crsymrma = crsymrma[-(which(is.na(rownames(crsymrma)))),]
crsymrma = collapseRows(crsymrma, unique(rownames(crsymrma)), unique(rownames(crsymrma)))
crsymrma = crsymrma$datETcollapsed
dim(crsymrma) #13031 genes

crsymmas = crmas
tmpgenes = mapIds(hgu133a.db, keys = rownames(crsymmas), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(crsymmas) = tmpgenes
crsymmas = crsymmas[-(which(is.na(rownames(crsymmas)))),]
crsymmas = collapseRows(crsymmas, unique(rownames(crsymmas)), unique(rownames(crsymmas)))
crsymmas = crsymmas$datETcollapsed
dim(crsymmas) #13031 genes

###############################################################################
## Fuentes-Mattei et al data:

files = readLines('./raw/fuentes-mattei/files.txt')
files = paste('./raw/fuentes-mattei/', files, sep='')
fmraw = ReadAffy(filenames = files)

fmrma = rma(fmraw) ## RMA normalise the data
fmrma = exprs(fmrma) ## change the format into matrix

fmmas = mas5(fmraw) ## Mas5 normalise the data
fmmas = exprs(fmmas) ## change the format into matrix
fmmas = log2(fmmas) ## log2 the data

fmstdrma = t(apply(fmrma, 1, function(x) (x-mean(x))/sd(x)))
fmstdmas = t(apply(fmmas, 1, function(x) (x-mean(x))/sd(x)))

fmobsgene = read.csv('./obsgenes/fmobsgenes.txt', header=T)

fmclin = read.csv('./clindata/fmclin.csv', sep=',', header=T)

fmsymrma = fmrma
tmpgenes = mapIds(hgu133a.db, keys = rownames(fmsymrma), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(fmsymrma) = tmpgenes
fmsymrma = fmsymrma[-(which(is.na(rownames(fmsymrma)))),]
fmsymrma = collapseRows(fmsymrma, unique(rownames(fmsymrma)), unique(rownames(fmsymrma)))
fmsymrma = fmsymrma$datETcollapsed
dim(fmsymrma) #13031 genes

fmsymmas = fmmas
tmpgenes = mapIds(hgu133a.db, keys = rownames(fmsymmas), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(fmsymmas) = tmpgenes
fmsymmas = fmsymmas[-(which(is.na(rownames(fmsymmas)))),]
fmsymmas = collapseRows(fmsymmas, unique(rownames(fmsymmas)), unique(rownames(fmsymmas)))
fmsymmas = fmsymmas$datETcollapsed
dim(fmsymmas) #13031 genes

###############################################################################
## Cris Print's Breast cancer data:

files = readLines('./raw/cris/files.txt')
files = paste('./raw/cris/', files, sep='')
crisraw = ReadAffy(filenames = files)

crisrma = rma(crisraw) ## RMA normalise the data
crisrma = exprs(crisrma) ## change the format into matrix

crismas = mas5(crisraw) ## Mas5 normalise the data
crismas = exprs(crismas) ## change the format into matrix
crismas = log2(crismas) ## log2 the data

crisstdrma = t(apply(crisrma, 1, function(x) (x-mean(x))/sd(x)))
crisstdmas = t(apply(crismas, 1, function(x) (x-mean(x))/sd(x)))

crisclin = read.csv('./clindata/crisclin2.csv', sep=',', header=T)

crissymrma = crisrma
tmpgenes = mapIds(hgu133a.db, keys = rownames(crissymrma), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(crissymrma) = tmpgenes
crissymrma = crissymrma[-(which(is.na(rownames(crissymrma)))),]
crissymrma = collapseRows(crissymrma, unique(rownames(crissymrma)), unique(rownames(crissymrma)))
crissymrma = crissymrma$datETcollapsed
dim(crissymrma) #13031 genes

crissymmas = crismas
tmpgenes = mapIds(hgu133a.db, keys = rownames(crissymmas), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(crissymmas) = tmpgenes
crissymmas = crissymmas[-(which(is.na(rownames(crissymmas)))),]
crissymmas = collapseRows(crissymmas, unique(rownames(crissymmas)), unique(rownames(crissymmas)))
crissymmas = crissymmas$datETcollapsed
dim(crissymmas) #13031 genes

###############################################################################
## Import my obesity genes from task 13:

allobsname = c("rawobsgenes","crolgenes","resobsgenes","rescrolgenes", "caobsgenes","cacrolgenes","caresobsgenes","carescrolgenes")
for (i in 1:length(allobsname)) {
	txt = paste('./obsgenes/', allobsname[i], sep='')
	txt = paste(txt, '.txt', sep='')
	assign(allobsname[i], dget(txt))
}

###############################################################################
## Make transformation matrix:

gtrmatransmat = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = gtrmastd[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	gtrmatransmat[[i]] = trans
	names(gtrmatransmat)[i] = paths[i]
}

gtmastransmat = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = gtmasstd[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	gtmastransmat[[i]] = trans
	names(gtmastransmat)[i] = paths[i]
}

mgfliprma = c( 'bcat_probes', 'e2f1_probes', 'egfr_probes', 'her2_probes', 'ifna_probes', 'ifng_probes', 'p53_probes', 'p63_probes', 'pi3k_probes', 'pr_probes', 'stat3_probes', 'tgfb_probes', 'tnfa_probes')

mgflipmas5 = c( 'akt_probes', 'bcat_probes', 'ifna_probes', 'ifng_probes', 'p53_probes', 'p63_probes', 'pi3k_probes', 'pr_probes', 'src_probes', 'stat3_probes', 'tgfb_probes', 'tnfa_probes')

###############################################################################
## Apply transformation matrix in Creighton's data

pdf(file='pdf/test.pdf', width=7, height=7)
test1 = gttransfun(crstdrma, paths, gtrmatransmat, mgfliprma, main = "Creighton data (RMA)")
test2 = gttransfun(crstdmas, paths, gtmastransmat, mgflipmas5, main = "Creighton data (MAS5)")
test3 = gttransfun(fmstdrma, paths, gtrmatransmat, mgfliprma, main = "FM data (RMA)")
test4 = gttransfun(fmstdmas, paths, gtmastransmat, mgflipmas5, main = "FM data (MAS5)")
test5 = gttransfun(crisstdrma, paths, gtrmatransmat, mgfliprma, main = "Cris data (RMA)")
test6 = gttransfun(crisstdmas, paths, gtmastransmat, mgflipmas5, main = "Cris data (MAS5)")
test7 = gttransfun(gtrmastd, paths, gtrmatransmat, mgfliprma, main = "Cris data (RMA)")
test8 = gttransfun(gtmasstd, paths, gtmastransmat, mgflipmas5, main = "Cris data (MAS5)")
dev.off()

# make all transformation matrices in Gatza data (for now)
allmeta = c(paths, allobsname)

allrmatransmat = list()
for (i in 1:length(allmeta)) {
	gene = get(allmeta[i])
	mat = gtrmastd[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	allrmatransmat[[i]] = trans
	names(allrmatransmat)[i] = allmeta[i]
}

allmastransmat = list()
for (i in 1:length(allmeta)) {
	gene = get(allmeta[i])
	mat = gtmasstd[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	allmastransmat[[i]] = trans
	names(allmastransmat)[i] = allmeta[i]
}

mgfliprmaall = c( 'bcat_probes', 'e2f1_probes', 'egfr_probes', 'her2_probes', 'ifna_probes', 'ifng_probes', 'p53_probes', 'p63_probes', 'pi3k_probes', 'pr_probes', 'stat3_probes', 'tgfb_probes', 'tnfa_probes', 'resobsgenes')

mgflipmas5all = c( 'akt_probes', 'bcat_probes', 'ifna_probes', 'ifng_probes', 'p53_probes', 'p63_probes', 'pi3k_probes', 'pr_probes', 'src_probes', 'stat3_probes', 'tgfb_probes', 'tnfa_probes', 'crolgenes', 'rescrolgenes', 'cacrolgenes', 'carescrolgenes')

pdf(file='pdf/test2.pdf', width=7, height=7)
test1 = gttransfun(crstdrma, allmeta, allrmatransmat, mgfliprmaall, main = "Creighton data (RMA)")
test2 = gttransfun(crstdmas, allmeta, allmastransmat, mgflipmas5all, main = "Creighton data (MAS5)")
test3 = gttransfun(fmstdrma, allmeta, allrmatransmat, mgfliprmaall, main = "FM data (RMA)")
test4 = gttransfun(fmstdmas, allmeta, allmastransmat, mgflipmas5all, main = "FM data (MAS5)")
test5 = gttransfun(crisstdrma, allmeta, allrmatransmat, mgfliprmaall, main = "Cris data (RMA)")
test6 = gttransfun(crisstdmas, allmeta, allmastransmat, mgflipmas5all, main = "Cris data (MAS5)")
test7 = gttransfun(gtrmastd, allmeta, allrmatransmat, mgfliprmaall, main = "Cris data (RMA)")
test8 = gttransfun(gtmasstd, allmeta, allmastransmat, mgflipmas5all, main = "Cris data (MAS5)")
dev.off()

###############################################################################
## Quick metagene direction check for the obesity associated genes in Gatza data:

commongenes = table(c(rawobsgenes,crolgenes,resobsgenes,rescrolgenes, caobsgenes,cacrolgenes,caresobsgenes,carescrolgenes))
commongenes = commongenes[commongenes==8]
commongenes = names(commongenes)

matheat = gatzabatchrma
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat = matheat[commongenes,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gtobrma.pdf')
for (i in 1:length(allobsname)) {
	genes = get(allobsname[i])
	tmp = gtrmastd[genes,]
	tmpsvd = svd(tmp)
	tmpmeta = tmpsvd$v[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	ord = order(tmpmeta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=allobsname[i], cexRow=1.0)
}
dev.off()

matheat = gatzabatchmas
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat = matheat[commongenes,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gtobmas.pdf')
for (i in 1:length(allobsname)) {
	genes = get(allobsname[i])
	tmp = gtmasstd[genes,]
	tmpsvd = svd(tmp)
	tmpmeta = tmpsvd$v[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	ord = order(tmpmeta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=allobsname[i], cexRow=1.0)
}
dev.off()

###############################################################################
## Look at the correlation between the metagene produced from transformation matrix and SVD

# check if function works:
test1 = getmeta(gtrmastd, paths, gtrmatransmat)
test2 = getmeta(gtmasstd, paths, gtmastransmat)
test3 = getmeta(gtmasstd, paths, gtrmatransmat)
test4 = getmeta(gtrmastd, paths, gtmastransmat)
# correlation is 1 for both test1 and test2 - as expected since the transformation matrices were derived in this dataset
# correlation is slightly off for test3 and test4 - probably due to the difference in the normalisation method

# Look at the correlation of the TM vs SVD metagenes in other data:
test5  = getmeta(fmstdrma, paths, gtrmatransmat)
test6  = getmeta(fmstdmas, paths, gtmastransmat)
test7  = getmeta(crstdrma, paths, gtrmatransmat)
test8  = getmeta(crstdmas, paths, gtmastransmat)
test9  = getmeta(crisstdrma, paths, gtrmatransmat)
test10 = getmeta(crisstdmas, paths, gtmastransmat)

txt = c('FM', 'Creighton', 'Cris')

gtrmatransres = cbind(test5$correlation, test7$correlation, test9$correlation)
colnames(gtrmatransres) = txt
gtmastransres = cbind(test6$correlation, test8$correlation, test10$correlation)
colnames(gtmastransres) = txt

###############################################################################
# Repeat the same thing, but with the BMI metagenes:

allobsname = c(allobsname, 'crobsgenes')

# make transformation matrix in Creighton's data:
obsrmatransmat = list()
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = crstdrma[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	obsrmatransmat [[i]] = trans
	names(obsrmatransmat)[i] = allobsname[i]
}

obsmastransmat = list()
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = crstdmas[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	obsmastransmat [[i]] = trans
	names(obsmastransmat)[i] = allobsname[i]
}

test11 = getmeta(crstdrma, allobsname, obsrmatransmat)
test12 = getmeta(crstdmas, allobsname, obsmastransmat)

# Look at the correlation of the TM vs SVD metagenes in other data:
test13 = getmeta(gtrmastd, allobsname, obsrmatransmat)
test14 = getmeta(gtmasstd, allobsname, obsmastransmat)
test15 = getmeta(fmstdrma, allobsname, obsrmatransmat)
test16 = getmeta(fmstdmas, allobsname, obsmastransmat)
test17 = getmeta(crisstdrma, allobsname, obsrmatransmat)
test18 = getmeta(crisstdmas, allobsname, obsmastransmat)

txt = c('Gatza', 'FM', 'Cris')

obsrmatransres = cbind(test13$correlation, test15$correlation, test17$correlation)
colnames(obsrmatransres) = txt
obsmastransres = cbind(test14$correlation, test16$correlation, test18$correlation)
colnames(obsmastransres) = txt

###############################################################################
# Have a quick look at the effect of using RMA on MAS5 data and vice versa

test3 = getmeta(gtmasstd, paths, gtrmatransmat)
test4 = getmeta(gtrmastd, paths, gtmastransmat)
# correlation is 1 for both test1 and test2 - as expected since the transformation matrices were derived in this dataset
# correlation is slightly off for test3 and test4 - probably due to the difference in the normalisation method
# test3 and test4 would be the 'reference' point for the expected correlation

test19  = getmeta(fmstdrma, paths, gtmastransmat)
test20  = getmeta(fmstdmas, paths, gtrmatransmat)
test21  = getmeta(crstdrma, paths, gtmastransmat)
test22  = getmeta(crstdmas, paths, gtrmatransmat)
test23  = getmeta(crisstdrma, paths, gtmastransmat)
test24 = getmeta(crisstdmas, paths, gtrmatransmat)

txt = c('FM', 'Creighton', 'Cris')

revmasonrmagt = cbind(test19$correlation, test21$correlation, test23$correlation)
colnames(revmasonrmagt) = txt
revrmaonmasgt = cbind(test20$correlation, test22$correlation, test24$correlation)
colnames(revrmaonmasgt) = txt

# try it with obesity data

# 'reference' correlation
test25 = getmeta(crstdrma, paths, gtmastransmat)
test26 = getmeta(crstdmas, paths, gtrmatransmat)

test27  = getmeta(gtrmastd,   allobsname, obsmastransmat)
test28  = getmeta(gtmasstd,   allobsname, obsrmatransmat)
test29  = getmeta(fmstdrma,   allobsname, obsmastransmat)
test30  = getmeta(fmstdmas,   allobsname, obsrmatransmat)
test31  = getmeta(crisstdrma, allobsname, obsmastransmat)
test32 = getmeta(crisstdmas,  allobsname, obsrmatransmat)

txt = c('FM', 'Creighton', 'Cris')

revmasonrmaobs = cbind(test27$correlation, test29$correlation, test31$correlation)
colnames(revmasonrmaobs) = txt
revrmaonmasobs = cbind(test28$correlation, test30$correlation, test32$correlation)
colnames(revrmaonmasobs) = txt











