###############################################################################

                        ##R code for task 13

###############################################################################

## Combine tasks 1 to 12 into this document

###############################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task13/')

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
## Load all the data using loadData.R:

source('loadData.R')

###############################################################################
## Creighton metagene stuff

source('creightonGE.R')

###############################################################################
## Fuentes-Mattei metagene stuff

source('fmmeta.R')

###############################################################################
## Gatza pathway stuff

source('gatzameta.R')

###############################################################################
## Common genes across ICGC data stuff

## use codes from task11?

###############################################################################
## Pathway enrichment stuff

## TODO: do the enrichment analysis with reactome and kegg

###############################################################################
## Conitnuous BMI DEG stuff

# get normalised Creighton microarray gene expression data
crcontmat = cr_raw
dim(crcontmat) # 22283 103

# Correlate it with the BMI value of the samples:
bmicor = cor(t(crcontmat), crclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"

# Try it with standardised data:
crstdmat = crcontmat
crstdmat = t(apply(crstdmat, 1, function(x) (x-mean(x))/sd(x)))

bmicorstd = cor(t(crstdmat), crclin$bmi, method = "spearman")
colnames(bmicorstd) = "correlation"

## The values didn't have any difference in the results
summary(bmicor)
summary(bmicorstd)

## Do some MClust stuff:
mc <- Mclust(bmicor,2)
y = seq(-0.4,0.4,l=100)
plot(y, type = 'l', dnorm(y, mc$parameters$mean[1], sqrt(mc$parameters$variance$sigmasq[1])), col = 'blue', ylim = c(0,6))
lines(y, dnorm(y, mc$parameters$mean[2], sqrt(mc$parameters$variance$sigmasq[2])), col = 'red')

# Do some data scrambling
set.seed(1) # set the seed for reproducibility


###############################################################################
## UCEC DEG stuff

datmat = UCEC
datbmi = UCECbmi
datmat = voomICGC(t(datmat), datbmi)

group = datbmi$bmiStatus
group = ifelse(group=='obese', 'obese', 'non-obese')
design = model.matrix(~group)
tt = make_tt(datmat,design)
ucecdegs = pull_deg(tt, adj = F)
adjucecdegs = pull_deg(tt, adj = T)

ucecoblist = rownames(ucecdegs)[1:799]
ucecadjoblist = rownames(adjucecdegs)

# check in UCEC data first
datmat = t(UCEC)
datmat = datmat[ucecoblist,]
datmatadj = datmat[ucecadjoblist,]
datmat    = standardise_data(datmat)
datmatadj = standardise_data(datmatadj)

# see log for why I'm flipping this metagene, but not the adjusted metagene
ucecsvd = svd(datmat)
ucecmeta = ucecsvd$v[,1]
ucecmeta = rank(ucecmeta)/length(ucecmeta)
ucecmeta = 1-ucecmeta
ucecord = order(ucecmeta)

ucecadjsvd = svd(datmatadj)
ucecadjmeta = ucecadjsvd$v[,1]
ucecadjmeta = rank(ucecadjmeta)/length(ucecadjmeta)
ucecadjord = order(ucecadjmeta)

maintxt = "UCEC unadjusted metagene"
maintxt2 = "UCEC adjusted metagene"

pdf('pdf/ucecmeta.pdf')
metaplot3(datmat, ucecmeta, datbmi, name=maintxt)
metaplot3(datmatadj, ucecadjmeta, datbmi, name=maintxt2)
dev.off()

# UCEC metagenes are both showing "dose-dependent" response
# (normal < overweight < obese).

# make transformation matrix:
ucecrawtransmat = diag(1/ucecsvd$d) %*% t(ucecsvd$u)
ucecadjtransmat = diag(1/ucecadjsvd$d) %*% t(ucecadjsvd$u)

# Try it on ICGC data first, then adjust the oblist for Creighton data:
tmp = paste(cancertypes,'bmi',sep='')

pdf('pdf/ucecicgc.pdf')
check_data(datlist=cancertypes, bmilist=tmp, genelist=ucecoblist, transmat=ucecrawtransmat, log=T, flip=T, name='UCEC metagene')
dev.off()

pdf('pdf/ucecadjicgc.pdf')
check_data(datlist=cancertypes, bmilist=tmp, genelist=ucecadjoblist, transmat=ucecadjtransmat, log=T, flip=F, name='UCEC adjusted metagene')
dev.off()

# seems like all the cancer types (except BLCA overweight samples) correlate
# with the UCEC obesity-associated genes

# adjust the oblists for Creighton data and check in the Creighton data (I
# don't think it'll work though)

## for each metagenes, identify the genes that are in the ICGC data and use these genes for validataion
cr_symmat = cr_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(cr_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_symmat) = tmpgenes
cr_symmat = cr_symmat[-(which(is.na(rownames(cr_symmat)))),]
cr_symmat = collapseRows(cr_symmat, unique(rownames(cr_symmat)), unique(rownames(cr_symmat)))
cr_symmat = cr_symmat$datETcollapsed
dim(cr_symmat) #13031 genes

crucecoblist    = ucecoblist[which(ucecoblist %in% rownames(cr_symmat))]
crucecadjoblist = ucecadjoblist[which(ucecadjoblist %in% rownames(cr_symmat))]

# make transformation matrix with adjusted metagenes
datmat = t(UCEC)
datmat = datmat[crucecoblist,]
datmat    = standardise_data(datmat)
datmatsvd  = svd(datmat)

datmatadj = datmat[crucecadjoblist,]
datmatadj = standardise_data(datmatadj)
datmatadjsvd = svd(datmatadj)

crucecrawtransmat = diag(1/ucecsvd$d) %*% t(ucecsvd$u)
crucecadjtransmat = diag(1/ucecadjsvd$d) %*% t(ucecadjsvd$u)

datmat = cr_symmat[crucecoblist,]
datmat = standardise_data(datmat, log=F)
tmpsvd = t(ucecrawtransmat %*% datmat)

datmatadj = cr_symmat[crucecadjoblist,]
datmatadj = standardise_data(datmatadj, log=F)
tmpsvd2 = t(ucecadjtransmat %*% datmatadj)

pdf('pdf/ucecmetacr.pdf')

dev.off()













################################################################################# Methylation stuff





























